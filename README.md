# ITSoneDB population pipeline

___
**ITSoneDB** is a database colleciting Eukayotic ITS1 sequencens and consistent taxonomic annotation.  
___  
<div align=justify>The designed pipeline integrates ad-hoc Python and BASH scripts and third-party tools (See the figure below).
In the initial step the ENA entries are locally downloaded and eukaryotic entries are extracted.  
From each entry specific information (i.e. accession number, version, description line and annotation under specific keys) are pulled out and stored in a TSV file and consistent sequence data are annotated in FASTA files. TSV and FASTA files are analyzed by two parallel procedures to extract or to infer the ITS1 location. The TSV files are parsed out to extract the annotation relative to ITS1 boundaries by means of a commonly used ITS1 synonyms dictionary.   
In parallel, HMM profiles for 18S and 5.8S rRNA genes are mapped on FASTA files by means of hmmsearch (HMMER 3.1) (right diagram part).  The ITS1 boundaries information obtained by both procedures are merged in order to produce the files needed to populate the database.</div>


![Alt text](ITSoneDB_Eukaryotes.tiff "Pipeline steps developed to generate ITSoneDB")