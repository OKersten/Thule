

# This Fasta was produced as follows
# A realigned bam was produced mapping the entirety of Razorbill sequence data to the reference genome
# The sample was called RAZ_Anc
# min and max depth half and double ave. depth -> 33.12
# angsd -C 50 -doFasta 2 -doCounts 1 -i RAZ_Anc.Assembly_Puffin_NUMT.realigned.bam -out Razorbill.NU.MT -setMinDepth 15 -setMaxDepth 67 -minMapQ 30 -minQ 30 -ref Assembly_Puffin_NU.MT.fasta
# -> 31,809,951 sites covered enough
# gzip Razorbill.NU.MT.fa
