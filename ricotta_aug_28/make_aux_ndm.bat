cat ricotta_TRAJ.fasta ricotta_TRBJ.fasta ricotta_TRDJ.fasta ricotta_TRGJ.fasta > ricotta_TRJ.fasta
cat ricotta_TRAV.fasta ricotta_TRBV.fasta ricotta_TRDV.fasta ricotta_TRGV.fasta > ricotta_TRV.fasta
cat ricotta_TRBD.fasta ricotta_TRDD.fasta > ricotta_TRD.fasta

annotate_j ricotta_TRJ.fasta ricotta_TR.aux
annotate_j ricotta_TRAJ.fasta ricotta_TRA.aux
annotate_j ricotta_TRBJ.fasta ricotta_TRB.aux
annotate_j ricotta_TRDJ.fasta ricotta_TRD.aux
annotate_j ricotta_TRGJ.fasta ricotta_TRG.aux

make_igblast_ndm ricotta_TRAV.fasta VA ricotta_TRA.ndm
make_igblast_ndm ricotta_TRBV.fasta VA ricotta_TRB.ndm
make_igblast_ndm ricotta_TRDV.fasta VA ricotta_TRD.ndm
make_igblast_ndm ricotta_TRGV.fasta VA ricotta_TRG.ndm

cp ricotta_TRA.ndm ricotta_TR.ndm
tail -n +2 ricotta_TRB.ndm >> ricotta_TR.ndm
tail -n +2 ricotta_TRD.ndm >> ricotta_TR.ndm
tail -n +2 ricotta_TRG.ndm >> ricotta_TR.ndm
