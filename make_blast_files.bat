make_igblast_ndm ricotta_TRAV_gapped.fasta VA ricotta_TRA.ndm
make_igblast_ndm ricotta_TRBV_gapped.fasta VB ricotta_TRB.ndm
make_igblast_ndm ricotta_TRDV_gapped.fasta VD ricotta_TRD.ndm
make_igblast_ndm ricotta_TRGV_gapped.fasta VG ricotta_TRG.ndm
cat ricotta_TRA.ndm >ricotta_TR.ndm
tail -n +2  ricotta_TRB.ndm >>ricotta_TR.ndm
tail -n +2  ricotta_TRD.ndm >>ricotta_TR.ndm
tail -n +2  ricotta_TRG.ndm >>ricotta_TR.ndm

annotate_j ricotta_TRAJ.fasta ricotta_TRA.aux
annotate_j ricotta_TRBJ.fasta ricotta_TRB.aux
annotate_j ricotta_TRDJ.fasta ricotta_TRD.aux
annotate_j ricotta_TRGJ.fasta ricotta_TRG.aux
cat ricotta_TRA.aux >ricotta_TR.aux
tail -n +2  ricotta_TRB.aux >>ricotta_TR.aux
tail -n +2  ricotta_TRD.aux >>ricotta_TR.aux
tail -n +2  ricotta_TRG.aux >>ricotta_TR.aux

