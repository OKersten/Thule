Let's GO

```bash


ls /XX/XX/*.bam | grep -v "IOM001" | grep -v "RAZ_Anc" > bams_good
cp /XX/XX/Assembly_Puffin_NU.MT.fasta .
cp /XX/XX/ChromosomeList .

cp /XX/XX/SFS/Results/sites2do .

mkdir Results

cat bams_good | awk -F '/' '{print $2}' | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
sbatch Hetero.sh ${sample} ;
done

```

Hetero.sh

```bash


mkdir temp.${1}_${SLURM_ARRAY_TASK_ID}
cp -r Assembly_Puffin_NU.MT.fasta temp.${1}_${SLURM_ARRAY_TASK_ID}/ #reference genome
cd temp.${1}_${SLURM_ARRAY_TASK_ID}/

Bamfile=$(ls ../../*.bam | grep ${1})

module purge 
module load SAMtools/1.9-GCC-8.2.0-2.31.1
samtools faidx Assembly_Puffin_NU.MT.fasta

module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="1"
GENOME_REF=Assembly_Puffin_NU.MT.fasta

ChromosomeEdit=$(echo "Scaffolds_chromosome_"${SLURM_ARRAY_TASK_ID}":")

grep -w Scaffolds_chromosome_${SLURM_ARRAY_TASK_ID} ../sites2do > sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do

sleep 120

angsd sites index sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do

sleep 120 

#run angsd
angsd -i ${Bamfile} -anc ${GENOME_REF} -ref ${GENOME_REF} -doSaf 1 -GL 1 -r ${ChromosomeEdit} -sites sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do -doCounts 1 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 2 -minMapQ 30 -minQ 30 -P ${Threads} -out ${1}_Chromosome${SLURM_ARRAY_TASK_ID}

mv ${1}* ../Results/

cd ..
rm -r temp.${1}_${SLURM_ARRAY_TASK_ID}/

```
Then

```bash

cd Results/

cat ../bams_good | awk -F '/' '{print $2}' | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
sbatch SFS.creation.sh ${sample} ;
done

```

SFS.creation.sh

```bash


module load angsd/0.931-GCC-8.2.0-2.31.1
Threads="16"

realSFS cat ${1}_Chromosome*.saf.idx -outnames ${1} ;
realSFS ${1}.saf.idx -P ${Threads} -fold 1 > ${1}.sfs ;

rm ${1}_Chromosome*.saf.gz
rm ${1}_Chromosome*.saf.idx
rm ${1}_Chromosome*.saf.pos.gz
rm ${1}_Chromosome*.arg 

```

```bash


module purge
module load R/3.6.0-intel-2019a

touch Heterozygosity.txt
cat ../bams_good | awk -F '/' '{print $2}' | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
Hetero=$(Rscript Hetero.R ${sample}.sfs) ;
echo -e ${sample}"\t"${Hetero} >> Heterozygosity.txt ;
done ;

```

Hetero.R

```r
infile <- commandArgs(trailingOnly = TRUE)
a<-scan(infile)
hetero <- a[2]/sum(a)
cat(hetero)
```
