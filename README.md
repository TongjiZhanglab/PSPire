# PSPire

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Version: 1.0.0](https://img.shields.io/badge/Version-1.0.0-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)

PSPire is a machine learning model based on integrated residue-level and structure-level features to predict phase-separating proteins. It is written in Python3 and is available as a command line tool.

# Installation

## Installation using Conda

1. Create conda environment:

   ```
   conda create -n PSPire python=3.8
   conda activate PSPire
   ```

2. Install [DSSP](https://github.com/PDB-REDO/dssp). Conda installation is recommended:

   ```shell
   conda install -c salilab -y dssp=2.2.1
   ```

3. Install [Pymol](https://pymol.org) for sticker calculation:

   ```shell
   conda install -c conda-forge -y pymol-open-source
   ```

4. Install the following python package:

   ```shell
   conda install -y xgboost=1.6.2 
   conda install -y scikit-learn=1.1.2 
   conda install -y biopython numpy pandas requests
   ```

5. (Optional) The published models of PSPire used PSAIA software to calculate relative solvent accessible surface area (RSA). It is recommended to prepare the running environment for PSAIA. User can first install Singularity and then use Singularity to build a container. If Singularity or the container image is not found, PSPire would use DSSP to calculate RSA.

   + Install [Singularity](https://apptainer.org/admin-docs/master/installation.html#). As Singularity is written primarily in Go, you should install [Go](https://go.dev/doc/install) first. After installation, you can type the command below to check Singularity has been installed successfully.

     ```shell
     singularity -h
     ```

     > Note: you need to use the following command to source the Singularity bash completion file to make sure the usage of bash completion in new shells.
     >
     > ```shell
     > echo ". Singularity_Installation_Path/etc/bash_completion.d/singularity" >> ~/.bashrc  # you should replce Singularity_Installation_Path with your installation path
     > ```

   + Build Qt 4.8.6 libraray container image needed for PSAIA software running:

   ```shell
   cd /path/to/PSPire/data
   singularity -d build psaia.simg docker://msoares/qt4-dev
   ```

6. Add PSPire to `$PATH`. To enable global access to PSPire from any location on your system, it's recommended to add the PSPire's directory to your system's `$PATH` environment variable. 

   ```bash
   chmod o+x /path/to/PSPire/PSPire.py
   echo 'export PATH=/path/to/PSPire:$PATH' >> ~/.bashrc
   source ~/.bashrc
   ```

## Installation using docker

1. Pull the PSPire docker image:

   ```shell
   docker pull houshuang2020/pspire:latest
   ```

2. Replace the `lib/psaia_run.py` and `PSPire.py` files with the files under the docker_script folder.
3. Make sure the PSAIA script is executable:
   ```
   chmod o+x /path/to/PSPire/software/PSAIA/psa
   ```
# Usage on Code Ocean
PSPire have been submitted to the Code Ocean platform (https://codeocean.com/capsule/3528836/tree/v2) and the required environment to run the code have been setup. Users can edit the "run" and "test.sh" files on the platform to predict proteins of interest.

# Parameters

Users can use `PSPire.py -h` to get help information of PSPire. 

Required parameters: the following parameters are required and you can only specify one of them.

```python
-u UNIPROT [UNIPROT ...], --uniprot UNIPROT [UNIPROT ...]
                      UniProt IDs. Multiple IDs should be separated by space.
-f FILE, --file FILE  List file with UniProt IDs or absolute path of protein pdb files.
                        Each ID or pdb file name should take one line.
-p PDBFILE, --pdbfile PDBFILE
                      PDB file of a protein.
-d DIRECTORY, --directory DIRECTORY
                      Absolute directory path of pdb files. The script will automatically
                        search files with pdb suffix under the specified directory.
```

Optional parameters:

```python
-o OUTPUT, --output OUTPUT
                      Output file name. (Default: standard out)
-n NAME, --name NAME  Project name. PSPire would use this name to create temporary
                        directory. (Default: PSPire)
--ignore              When this parameter is set, PSPire would ignore intrinsically
                        disordered region(IDR)-related features for proteins with IDRs.
-s PHOS, --phos PHOS  Absolute path of the user-defined phosphorylation (Phos) feature file.
                        If this parameter is specified, PSPire would use the model with the
                        Phos feature. User can check the demo directory of PSPire software 
                        package for example format of the phos feature file. (Default: '')
--mobidb              By default, PSPire would assume the pdb files you provide have pLDDT 
                        score in the B-factor column calculated by AlphaFold, and use the 
                        score to get IDRs. When this parameter is set, PSPire would get IDRs 
                        by MobiDB-lite software.
-t THRESHOLD, --threshold THRESHOLD
                      Threshold of pLDDT score to get idr regions. (Default: 50)
-c CUTOFF, --cutoff CUTOFF
                      If the RSA percentage of a residue is greater than this cutoff, it 
                        will be assigned as exposed surface residue, otherwise as buried 
                        residue. (Default: 25)
-j JOBS, --jobs JOBS  If mobidb parameter is on, PSPire would use the given number of cpus 
                        to run MobiDB-lite. (Default: 10)
--resume              By default, PSPire would clean up the temporary files and start from 
                        the beginning. When resume is on, each re-run would use previous 
                        temporary files to resume from the step it crashed.
--dont_remove         By default, PSPire would clean up temporary files. When dont_remove is 
                        on, PSPire would keep temporary files.
```

## Usage

1. Specify uniprot ids:

   ```shell
   PSPire.py -u P09651
   PSPire.py -u P09651 O00444
   ```

2. Specify list file with UniProt IDs or absolute path of protein pdb files:

   ```shell
   PSPire.py -f ${SOFTWAREPATH}/demo/PDB_files_list.txt
   PSPire.py -f ${SOFTWAREPATH}/demo/uniprotID_list.txt
   ```

3. Specify PDB file of a protein:

   ```shell
   PSPire.py -p ${SOFTWAREPATH}/demo/AF-A0A2R8QUZ1-F1-model_v2.pdb
   ```

4. Specify absolute directory path of pdb files:

   ```shell
   PSPire.py -d ${SOFTWAREPATH}/demo
   ```

5. Ignore intrinsically disordered region(IDR)-related features for proteins with IDRs:

   ```shell
   PSPire.py -u P09651 --ignore
   ```

6. Specified user-defined phosphorylation feature:

   ```shell
   PSPire.py -u P09651 -s ${SOFTWAREPATH}/demo/phos_feature_example.csv
   ```

7. Use MobiDB-lite software to calculate IDRs:

   ```shell
   PSPire.py -u P09651 --mobidb
   ```

# Output

By default, the results will be sent to the standard output. User can use the "-o" parameter to specify output file name in CSV format. The "Score" column represents the PS propensity of the protein, indicating that a higher value corresponds to a higher likelihood of PS. The "Include_IDRs" column indicates whether the protein contains IDRs. The last three columns indicate the structured superficial regions, positive and negative sticker regions.

```python
Uniprot_ID,Score,Include_IDRs,SSUP_regions,Pos_Sticker_Regions,Neg_Sticker_Regions
P45973,0.9956347,Yes,"1-21,23-24,26-27,29,31-34,42-48,50-51,54-56,58,60-62,64-65,68-69,71-80,113-115,118-122,124-125,127,129-132,134,136,139,141,143-144,146-148,152,154-155,157-159,161-162,165-166,168-191","[(2, 3, 4, 5, 6, 7), (74, 75, 76, 77, 79), (29, 152)]","[(179, 180, 181), (12, 13, 14, 15, 16, 17), (18, 19, 20, 21, 23, 24, 42, 50, 56, 58)]"
```

## Pre-calculated scores

Pre-calculated PSPire predicted scores of the following model organism proteomes using protein structures provided by [AlphaFold](https://alphafold.ebi.ac.uk/download) were provided under the pre_calculated_scores folder. The model incorporating the Phos feature was specifically utilized for human proteins, while the model without the Phos feature was employed for proteins from other species. Additionally, the files also contain the predicted scores of proteins when the IDR-related features are ignored. As for human, the file also contains relative surface exposure and secondary structure state data of each residue.

- [Arabidopsis thaliana](https://compbio-zhanglab.org/release/PSPire_scores/Arabidopsis_thaliana_scores.csv)
- [Caenorhabditis elegans](https://compbio-zhanglab.org/release/PSPire_scores/Caenorhabditis_elegans_scores.csv)
- [Candida albicans](https://compbio-zhanglab.org/release/PSPire_scores/Candida_albicans_scores.csv)
- [Danio rerio](https://compbio-zhanglab.org/release/PSPire_scores/Danio_rerio_scores.csv)
- [Dictyostelium discoideum](https://compbio-zhanglab.org/release/PSPire_scores/Dictyostelium_discoideum_scores.csv)
- [Drosophila melanogaster](https://compbio-zhanglab.org/release/PSPire_scores/Drosophila_melanogaster_scores.csv)
- [Escherichia coli](https://compbio-zhanglab.org/release/PSPire_scores/Escherichia_coli_scores.csv)
- [Glycine max](https://compbio-zhanglab.org/release/PSPire_scores/Glycine_max_scores.csv)
- [Homo sapiens](https://compbio-zhanglab.org/release/PSPire_scores/Homo_sapiens_scores.csv)
- [Methanocaldococcus jannaschii](https://compbio-zhanglab.org/release/PSPire_scores/Methanocaldococcus_jannaschii_scores.csv)
- [Mus musculus](https://compbio-zhanglab.org/release/PSPire_scores/Mus_musculus_scores.csv)
- [Oryza sativa](https://compbio-zhanglab.org/release/PSPire_scores/Oryza_sativa_scores.csv)
- [Rattus norvegicus](https://compbio-zhanglab.org/release/PSPire_scores/Rattus_norvegicus_scores.csv)
- [Saccharomyces cerevisiae](https://compbio-zhanglab.org/release/PSPire_scores/Saccharomyces_cerevisiae_scores.csv)
- [Schizosaccharomyces pombe](https://compbio-zhanglab.org/release/PSPire_scores/Schizosaccharomyces_pombe_scores.csv)
- [Zea mays](https://compbio-zhanglab.org/release/PSPire_scores/Zea_mays_scores.csv)

## Citing PSPire
If you use the code or data in this repository, please cite:

```
@article{Hou2024,
  author = {Hou, Shuang and Hu, Jiaojiao and Yu, Zhaowei and Li, Dan and Liu, Cong and Zhang, Yong},
  title = {Machine learning predictor PSPire screens for phase-separating proteins lacking intrinsically disordered regions},
  journal = {Nature Communications},
  volume = {15},
  number = {1},
  pages = {2147},
  year = {2024},
  doi = {10.1038/s41467-024-46445-y},
  url = {https://doi.org/10.1038/s41467-024-46445-y},
}
```
