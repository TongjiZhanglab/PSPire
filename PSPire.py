#!/usr/bin/env python
# Date: 2023-07-17
# Author: Shuang Hou
# Contact: houshuang@tongji.edu.cn


"""
Description: 
PSPire is a machine learning model based on integrated residue-level and structure-level features to predict phase-separating proteins.

Usage: -u/-f/-p/-d are required and you can only specify one of them.
1. ${SOFTWAREPATH}/PSPire.py -u P09651
2. ${SOFTWAREPATH}/PSPire.py -u P09651 O00444
3. ${SOFTWAREPATH}/PSPire.py -f ${SOFTWAREPATH}/demo/PDB_files_list.txt  #please change the PDB file path before running
4. ${SOFTWAREPATH}/PSPire.py -f ${SOFTWAREPATH}/demo/uniprotID_list.txt
5. ${SOFTWAREPATH}/PSPire.py -p ${SOFTWAREPATH}/demo/AF-A0A2R8QUZ1-F1-model_v2.pdb
6. ${SOFTWAREPATH}/PSPire.py -d ${SOFTWAREPATH}/demo
"""

from distutils.log import warn
import logging
import argparse
import os
import sys
import re

import requests
import pandas as pd
import numpy as np

from Bio import PDB
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser

from lib.Utility import *
from lib.get_pdb_info import *
from lib.dssp_run import dssp_run
from lib.psaia_run import psaia_run
from lib.calculate_features import *
from lib.predict import *

# ------------------------------------
# parse arguments
# ------------------------------------
description = 'PSPire can predict phase-separating probability of proteins based on integrated residue-level and structure-level features.'
parser = argparse.ArgumentParser(description=description)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-u', '--uniprot', type=str, nargs='+',
                   help='UniProt IDs. Multiple IDs should be separated by space.')
group.add_argument('-f', '--file', type=argparse.FileType('r'),
                   help="List file with UniProt IDs or absolute path of protein pdb files. Each ID or pdb file name should take one line.")
group.add_argument('-p', '--pdbfile', type=str,
                   help="PDB file of a protein.")
group.add_argument('-d', '--directory', type=str,
                   help="Absolute directory path of pdb files. The script will automatically search files with pdb suffix under the specified directory.")

parser.add_argument('-o', '--output', default=sys.stdout,
                    help="Output file name. (Default: standard out)")
parser.add_argument('-n', '--name', type=str, default='PSPire',
                    help="Project name. PSPire would use this name to create temporary directory. (Default: PSPire)")
parser.add_argument('--ignore', action="store_true",
                    help="When this parameter is set, PSPire would ignore intrinsically disordered region(IDR)-related features for proteins with IDRs.")
parser.add_argument('-s', '--phos', type=str, default='',
                    help='''Absolute path of the phos feature file. If this parameter is specified, PSPire would use the model with the Phos feature. 
                            User can check the demo directory of PSPire software package for example format of the phos feature file. (Default: '')''')
parser.add_argument('--mobidb', action="store_true",
                    help='''By default, PSPire would assume the pdb files you provide have pLDDT score in the B-factor column calculated by AlphaFold, 
                            and use the score to get IDRs. When this parameter is set, PSPire would get IDRs by MobiDB-lite software.''')
parser.add_argument('-t', '--threshold', default=50, type=int,
                    help='Threshold of pLDDT score to get idr regions. (Default: 50)')
parser.add_argument('-c','--cutoff', default=25, type=int,
                    help='If the RSA percentage of a residue is greater than this cutoff, it will be assigned as exposed surface residue, otherwise as buried residue. (Default: 25)')
parser.add_argument('-j', '--jobs', default=10, type=int,
                    help='If mobidb parameter is on, PSPire would use the given number of cpus to run MobiDB-lite. (Default: 10)')
parser.add_argument('--resume', dest="resume", action="store_true",
                    help='''By default, PSPire would clean up the temporary files and start from the beginning. 
                    When resume is on, each re-run would use previous temporary files to resume from the step it crashed.''')
parser.add_argument('--dont_remove', dest="nremove", action="store_true",
                    help='By default, PSPire would clean up temporary files. When dont_remove is on, PSPire would keep temporary files.')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filemode="w",
                    filename=args.name + '_running.log')

# ------------------------------------
# global parameters
# ------------------------------------
requests.packages.urllib3.disable_warnings()
PWDIR = os.getcwd()
PACKAGEDIR = sys.path[0]
tmpDir = PWDIR + '/' + args.name
make_dir(tmpDir,args.resume)


# ------------------------------------
# check input arguments
# ------------------------------------
def download_pdb_file(uid, dir):
    '''Download pdb file from AlphaFold database to a new folder'''
    filename = 'AF-' + uid + '-F1-model_v2.pdb'
    r = requests.get(f'https://alphafold.ebi.ac.uk/files/{filename}',headers={'Connection':'close'},verify=False)
    if r.status_code != 200:
        warn(f'Could not find pdb file for {uid} in AlphaFold database. Skip this file.')
        return False
    else:
        content = r.content.decode('utf-8')
        with open(dir + '/' + filename, 'w') as f:
            f.write(content)
        return True


def check_pdb_format(file):
    '''check whether the input file is in valid pdb format'''
    if not os.path.isfile(file):
        return False
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('_', file)
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                if residue.has_id("CA"):
                    return True
                else:
                    return False


def get_lists():
    ori_ids = list()
    ori_files =list()
    uniprot_ids = list()
    inputfiles = list()
    if args.uniprot:
        uniprot_ids = args.uniprot
        ori_ids = args.uniprot
    elif args.file:
        lists = [i.strip() for i in args.file.readlines()]
        if re.match(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', lists[0]):
            ori_ids = lists
            uniprot_ids.append(lists[0])
            if len(lists) > 1:
                for i in lists[1:]:
                    if re.match(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', i):
                        uniprot_ids.append(i)
                    else:
                        warn(f'{i} is not a valid UniProt ID. Skip this ID.')
        elif check_pdb_format(lists[0]):
            ori_files = lists
            inputfiles.append(lists[0])
            if len(lists) > 1:
                for i in lists[1:]:
                    if check_pdb_format(i):
                        inputfiles.append(i)
                    else:
                        warn(f'{i} is not in valid pdb format. Skip this file.')
        else:
            raise_error(f'{args.file} is not valid. List file should contain UniProt IDs or absolute path of protein pdb files. Each ID or pdb file name should take one line.')
    elif args.pdbfile:
        if check_pdb_format(args.pdbfile):
            inputfiles = [args.pdbfile]
        else:
            raise_error(f'{args.pdbfile} is not in valid pdb format.')
    else:
        ori_files = [os.path.join(args.directory, f) for f in os.listdir(args.directory) if os.path.isfile(os.path.join(args.directory, f)) and f.endswith('.pdb')]
        # out = run_cmd('ls '+os.path.join(args.directory, '*pdb'), 'result')
        if len(ori_files) == 0:
            raise_error(f'There is no pdb file under the specified directory: {args.directory}')
        else:
            # tmpfiles = out.split('\n')
            # ori_files = tmpfiles
            for i in ori_files:
                if check_pdb_format(i):
                    inputfiles.append(i)
                else:
                    warn(f'{i} is not in valid pdb format. Skip this file.')

    # download or link pdb files
    files = list()
    dirName = tmpDir + '/pdbFiles'
    make_dir(dirName,args.resume)
    if inputfiles:
        for i in inputfiles:
            fileName = dirName + '/' + os.path.basename(i)
            if not os.path.isfile(fileName):
                run_cmd('cp '+i+' '+dirName)
            files.append(fileName)

    if uniprot_ids:
        for uid in uniprot_ids:
            fileName = dirName + '/AF-'+uid+'-F1-model_v2.pdb'
            if not os.path.isfile(fileName):
                if download_pdb_file(uid, dirName):
                    files.append(fileName)
            elif not check_pdb_format(fileName):
                warn(f'{fileName} already exists but is not a valid pdb file. Re-download the file.')
                if download_pdb_file(uid, dirName):
                    files.append(fileName)
            else:
                info(f'{fileName} already exists.')
                files.append(fileName)
    
    if len(files)==0:
        warn('There is no id to analyze.')
        sys.exit(1)

    if type(args.output) == str:
        outputFile = open(args.output, 'w')
    else:
        outputFile = args.output
    
    if os.path.basename(files[0]).startswith('AF'):
        uids = [os.path.basename(i).replace('AF-', '').replace('-F1-model_v2.pdb', '') for i in files]
    else:
        uids = [os.path.basename(i).replace('.pdb', '') for i in files]
    
    if len(ori_files)>0:
        if os.path.basename(ori_files[0]).startswith('AF'):
            ori_ids = [os.path.basename(i).replace('AF-', '').replace('-F1-model_v2.pdb', '') for i in ori_files]
        else:
            ori_ids = [os.path.basename(i).replace('.pdb', '') for i in ori_files]
        
    records = list()
    for num,pdbFile in enumerate(files):
        pdbparser = PDBParser()
        structure = pdbparser.get_structure('aa', pdbFile)
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
        seq = chains['A']
        out = SeqRecord(Seq(seq),id='sp|'+uids[num]+'|'+uids[num])
        records.append(out)
        
    fastaFile = tmpDir+'/allSeq.fasta'
    SeqIO.write(records, fastaFile, "fasta")

    return list(uids), list(files), outputFile, ori_ids, ori_files,fastaFile


# ------------------------------------
# main functions
# ------------------------------------
def mobidb_run(fastaFile,ofile,jobs):
    make_dir(tmpDir + '/mobidb_result',args.resume)
    if os.path.isfile(ofile):
        info(f'{ofile} already exists.')
    else:
        mobidbExe = PACKAGEDIR+'/software/MobiDB-lite-3.8.4/mobidb_lite.py'
        cmd = 'python '+mobidbExe+' '+fastaFile+' -o '+ofile+' -t ' + str(jobs)
        run_cmd(cmd)


def main():
    uids, files, outputFile, ori_ids, ori_files,fastaFile = get_lists()
    
    mobidb_ofile = tmpDir + '/mobidb_result/mobidb.out'
    dsspFile = tmpDir + '/dssp_result/dssp_result.pkl'
    pdbFile = tmpDir + '/pdbInfo.pkl'
    
    if os.path.isfile(tmpDir + '/final_features.pkl') and os.path.isfile(pdbFile):
        features = pd.read_pickle(tmpDir + '/final_features.pkl')
        pdbInfo = pd.read_pickle(pdbFile)
    else:
        if not tool_exists('singularity'):
            warn('''Singularity is not callable from command line.
                 PSPire would use DSSP to calculate relative solvent accessible surface area (RSA). 
                 If you want to use PSAIA to calculate RSA as the model published in our paper, 
                 please install singularity and add it to your PATH.''')
        
        simgFile = PACKAGEDIR + '/data/psaia.simg'
        if not os.path.isfile(simgFile):
            warn('''Container image file needed for PSAIA running does not exist. 
                 PSPire would use DSSP to calculate relative solvent accessible surface area (RSA). 
                 If you want to use PSAIA to calculate RSA as the models published in our paper, 
                 please build the container image and make sure the image file under the suggested folder.''')
        
        if tool_exists('singularity') and os.path.isfile(simgFile):
            psaia_job = MyThread(psaia_run, args=(files, tmpDir, PACKAGEDIR, args.resume))
            psaia_job.start()
            info('PSAIA start running')

        if args.mobidb:
            mobidb_job = MyThread(mobidb_run, args=(fastaFile,mobidb_ofile,args.jobs))
            mobidb_job.start()
            info('Mobidb start running')
        
        if not args.mobidb or not (tool_exists('singularity') and os.path.isfile(simgFile)):
            dssp_job = MyThread(dssp_run, args=(uids, files, tmpDir, args.resume,dsspFile))
            dssp_job.start()
            info('DSSP start running.')
        
        if args.mobidb:        
            mobidb_job.join()
            
        if not args.mobidb or not (tool_exists('singularity') and os.path.isfile(simgFile)):
            dssp_job.join()

        if os.path.isfile(pdbFile):
            pdbInfo = pd.read_pickle(pdbFile)
        else:
            info('PDB info start running.')
            pdbInfo, files, uids = get_pdb_info(mobidb_ofile,dsspFile, uids, files, args.threshold, tmpDir, args.mobidb)
        
        if tool_exists('singularity') and os.path.isfile(simgFile):
            psaia_job.join()

        features = calculate_features(tmpDir, args.cutoff, files, uids, pdbInfo, args.ignore,dsspFile)
        
    if len(features) == 0:
        raise_error('The features is empty. Can not predict.')
        
    if args.phos:
        if os.path.isfile(args.phos):
            phos_fea = pd.read_csv(args.phos)
            features = pd.merge(features,phos_fea,on='uniprot_id',how='left')
        else:
            raise_error(f'The phos file does not exist: {args.phos}')
    
    dataDir = PACKAGEDIR + '/data/'
    outData = predict(dataDir,features)
    
    idrIds = pdbInfo[pdbInfo['idr_length_filter_mobidb']!=0].index.tolist()
    outData['Include_IDRs'] = np.where(outData['uniprot_id'].isin(idrIds),'Yes','No')
    cols = ['uniprot_id','ssup_regions','pos_sticker_regions','neg_sticker_regions']
    outData = pd.merge(outData,features[cols],on='uniprot_id',how='left')
    outData.columns = ['Uniprot_ID','Score','Include_IDRs','SSUP_regions','Pos_Sticker_Regions','Neg_Sticker_Regions']
    otable = outData.to_csv(index=False)
    outputFile.write(otable)
    
    if not args.nremove:
        run_cmd('rm -rf '+tmpDir)
        
    nids = ','.join(list(set(ori_ids)-set(uids)))
    if nids:
        if args.directory or (args.file and len(ori_files)!=0):
            print(f'The pdb file format of the following ids are not valid: {nids}')
        if args.uniprot or (args.file and len(ori_files)==0):
            print(f'The following ids is not valid or have no AlphaFold structure: {nids}')
        

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        error("User interrupt me ^_^ \n")
        sys.exit(1)
