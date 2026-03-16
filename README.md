# ITSoneDB population pipeline

___
**ITSoneDB** is a database collecting Eukayotic ITS1 sequences and consistent taxonomic annotations. It is available at [http://itsonedb.cloud.ba.infn.it/](http://itsonedb.cloud.ba.infn.it/). 
___
## RATIONALE  
<div align=justify>The pipeline designed for ITSoneDB population integrates <em>ad-hoc</em> Python and BASH scripts and third-party tools (See the figure below).

1. In the initial step the ENA entries are locally downloaded and eukaryotic entries are extracted.  

2. From each entry specific information (i.e. accession number, version, description line and annotation under specific keys) are pulled out and stored in a TSV file and consistent sequence data are annotated in FASTA files. TSV and FASTA files are analyzed by two parallel procedures to extract or to infer the ITS1 location. The TSV files are parsed out to extract the annotation relative to ITS1 boundaries by means of a commonly used ITS1 synonyms dictionary.
   
3. In parallel, HMM profiles for 18S and 5.8S rRNA genes are mapped on FASTA files by means of hmmsearch (HMMER 3.1) (right diagram part).  The ITS1 boundaries information obtained by both procedures are merged in order to produce the files needed to populate the database.</div>

___
![Alt text](ITSoneDB_Eukaryotes.png)

___

## REQUIREMENTS
1. The **HMMER** is required (for conda installation info see [https://anaconda.org/bioconda/hmmer](https://anaconda.org/bioconda/hmmer)).
2. The Species Representative Entry Identification procedure require **VSEARCH** (for installation info see [https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)).
3. Python3 dependencies are **argparse**, **argcomplete**, **getopt**, **multiprocessing**, **numpy**, **biopython**, **datetime**.

## USAGE
Following the instruction to execute the scripts:
```
$./ITSoneDB_upgrade_pipeline.sh
        -s full path directory containing managed scripts
        -x full path auxiliary files directory
        -r full path previous releases directory
        -c cpus number
        
$./ETL_population_pipeline.sh
        -s full path directory containing managed scripts
        -n release number
        -p previous release number
        -r full path releases directory
        
```
Auxiliary files are 16S and 5.8S HMM models and txt file containing ITS1 synonyms dictionary.
 

# ITSoneDB_devel

This is the ITSoneDB population pipeline used from 1.138 database version.

## File System Description

ITSoneDB_devel repository contains several forlders:

- **_ITSoneDB_upgrade_pipeline_**: it contains the python scripts that are menaged by 
**ITSoneDB_upgrade_pipeline.sh**

- **_ITSoneDB_etl_popul_pipeline_** : it contains the python scripts that are menaged by 
**ETL_population_pipeline.sh**

- _librs_ : it contains libraries used by python scripts

- _additionals_ : additional python scripts not updated


## **ITSoneDB_upgrade_pipeline.sh** usage

This script menages the entire procedure to produce ITSoneDB annotations on ENA sequences.
It verifies whether there is a newer ENA version with respect to the latest in _releases_ directory. 
If true, the program downloads in local the ENA flat files. Then applyies the procedure to extract ITS1 location from ENA notations and de-novo annotate ITS1 from NHMMER mapping 18S and 5.8 HMM profiles.
Then, a final table is produced and represents the starting point of **ETL_population_pipeline.sh** script.

#### USAGE:
```
$./ITSoneDB_upgrade_pipeline.sh
        -s full path to ITSoneDB_upgrade_pipeline directory
        -x full path to aux directory (containing auxiliary files such as 18S and 5.8S HMM profiles and 
           ITS1 annotation synonyms)
        -r full path to releases directory (containing previous releases of ITSoneDB)
        -c cpus number
```

## **ETL_population_pipeline.sh** usage

This script menages the procedure to populate ETL directory after
ITSoneDB_upgrade_pipeline procedure completion. The ETL directory contains tables and files
which have to be charged on [ITSoneDB web site](http://itsonedb.cloud.ba.infn.it/)

#### USAGE:

```
$./ETL_population_pipeline.sh
        -s full path to ITSoneDB_etl_popul_pipeline directory
        -n new release number
        -p previous release number
        -r full path to releases directory (containing all ITSoneDB releases)
```
