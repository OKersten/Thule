
Remove mitochondrial sequence from bam

```bash

cat Modern.samples | grep -v "RAZ_Ancestral" | grep -v "IOM001" | while read i ;
do sbatch BamToNU_bam.sh ${i} ; done

rm slurm*
rm IOM001*
```

BamToNU_bam.sh

```bash

module load SAMtools/1.9-GCC-8.2.0-2.31.1

samtools index ${1}.Assembly_Puffin_NUMT.realigned.bam
samtools view -@ 8 -b -h -o ${1}.Assembly_Puffin_NU.bam ${1}.Assembly_Puffin_NUMT.realigned.bam Scaffolds_chromosome_10 Scaffolds_chromosome_11 Scaffolds_chromosome_12 Scaffolds_chromosome_13 Scaffolds_chromosome_14 Scaffolds_chromosome_15 Scaffolds_chromosome_16 Scaffolds_chromosome_17 Scaffolds_chromosome_18 Scaffolds_chromosome_19 Scaffolds_chromosome_1 Scaffolds_chromosome_20 Scaffolds_chromosome_21 Scaffolds_chromosome_22 Scaffolds_chromosome_23 Scaffolds_chromosome_24 Scaffolds_chromosome_25 Scaffolds_chromosome_2 Scaffolds_chromosome_3 Scaffolds_chromosome_4 Scaffolds_chromosome_5 Scaffolds_chromosome_6 Scaffolds_chromosome_7 Scaffolds_chromosome_8 Scaffolds_chromosome_9 Scaffolds_chromosome_Z Scaffolds_unplaced
samtools sort -@ 8 ${1}.Assembly_Puffin_NU.bam > ${1}.Assembly_Puffin_NU.sorted.bam
samtools index ${1}.Assembly_Puffin_NU.sorted.bam
rm ${1}.Assembly_Puffin_NU.bam
rm ${1}.Assembly_Puffin_NUMT.realigned.bam
rm ${1}.Assembly_Puffin_NUMT.realigned.bam.bai
```

ANGSD - QualityCheck

```bash
# listing files for ANGSD to work on:
ls *.bam > bams

# get reference for ANGSD commands
cp /XX/XX/Assembly_Puffin_NU.MT.fasta .

sbatch Puffin_ANGSD_QC.sh

```

Puffin_ANGSD_QC.sh

```bash

mkdir tmp
cp -r *bam* tmp/
cp -r *fasta* tmp/
cd tmp/

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 25 -maxDepth 1078 -checkBamHeaders 1 -C 50 -baq 2 "

TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

module load SAMtools/1.9-GCC-8.2.0-2.31.1
samtools faidx Assembly_Puffin_NU.MT.fasta

module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="12"

angsd -b bams -GL 1 $FILTERS $TODO -P ${Threads} -ref Assembly_Puffin_NU.MT.fasta -r Scaffolds_chromosome_1: -out PuffQC

cp -r PuffQC* ../
cd ../
rm -r tmp/
```

THEN

```bash

gunzip PuffQC.counts.gz 
#subsample counts for every 10th SNP starting after header (reducing file size so R can handle - should be about 20M SNPs)
head -1 PuffQC.counts > PuffQC.header
awk 'NR % 10 == 2' PuffQC.counts > PuffQC.counts2  
cat PuffQC.header PuffQC.counts2 > PuffQC.counts3
wc -l PuffQC.counts3 #used for Rscript, insert as no. of rows for last graph - read.table
#edit and run plotQC.R # using cannibalized Matteo Fumagalli's script
Rscript ./plotQC.R PuffQC > qranks #warning messages are fine but need to adjust the number of individuals/columns and n's for last graph

rm PuffQC.counts2
rm PuffQC.counts3
gzip PuffQC.counts 

# percentages of sites with coverage >5x in each sample, from worst to best:
cat qranks
# look at PuffQC.pdf for more details

```