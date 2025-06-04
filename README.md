# notes on data for tr ref

## reference sets
- current version in ./ricotta_june_2
- genomic allele tables are in a sub-directory

## 10x data
- see notes.txt in ./all_10x_projects

## MUSA-style file for the TR loci
- in directories under ./db
- initially created with imgt sequences using imgt_to_db.py
- updated with genomic sequences using allele_table_to_ref.py

## TODO - Coding
- Submit ricotta assemblies to the SRA and create a table mapping repertoires to SRRs so that we can list the 'best' SRR supporting each allele (OGRDB submission needs this)
- Update tcr_db.csv with usage information from annotated repertoires: include number of supporting repertoires, max supporting repertoire (SRR), max supporting count
- Produce a submission sheet for IUIS - talk to Corey about what he needs, use the columns in the allele tabke as a starting point. May be enough just to filter this to the 'best' example for each allele that meets the submission criteria (see code in allele_table_to_ref.py for the criteria)
- Produce submission sheets for OGRDB. See notes.txt in ./sample_reports for what's required.
- Create some quick overview reports,e.g. how many alleles do we have of each type, how does this compare with IMGT (missing/novel). How many supported by repertoires.

## TODO - Analysis
- Review all alleles with subsequence/supersequence notes in tcr_db.csv. Be on the lookout for coordinate errors in the reference, eg RSS wrongly demarcated.
- Review the genomically derived L-PART1 and 2 against the leaders expressed in 10x data, to identify and fix any coordinate problems in the reference.
- Review IMGT functional alleles (probably F, not (F) or [F]) that are not included in our set




  




