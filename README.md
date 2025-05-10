# rough notes on 10x data for tr ref

## source 
- PRJNA603748. This has data from all 4 loci for ~ 40 subjects.
- 'ricotta' genomic analysis

## personalised set from ricotta data
- extract refs for each locus
- combine into a single TR set:
- python/calc_support_and_set.py
- igblast dbs in ricotta/db
- ricotta/ricotta_TR.aux
- ricotta/ricotta_TR.ndm
  
## full-length VDJ sequences from PRJNA603748
- PRJNA603748/SRRs.txt - list of sequence accessions
- PRJNA603748/download_sras: download sequence sets, using index into SRRs.txt. Uses SRA Toolkit (available as am amaconda package)
- PRJNA603748/process_PRJNA603748.slurm: processes the sets using SRA toolkit and cellranger. Uses a slurm array wth indexes into SRRs.txt.
- python/aggregate_contigs.py - pulls full_length VDJ sequences from all samples into a single file

## annotation
- details in PRJNA603748/notes.txt
- results are in dropbox (ask me for the links)
- python/calc_igblast_support.py. This determines which reference set alleles are supported in igblast, looking for unambigous calls with 100% identity.


  




