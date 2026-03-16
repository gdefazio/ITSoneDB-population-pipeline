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
