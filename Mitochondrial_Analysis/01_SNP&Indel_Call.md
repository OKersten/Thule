
Get all bam files from Paleomix run and extract mitochondrial mapping from bam files

```bash

cat Modern.samples | grep -v "RAZ_Ancestral" | while read i ;
do echo ${i} ;
sbatch IndexSortMT.sh ${i} ;
done ;

```

IndexSortMT.sh

```bash


module load SAMtools/1.9-GCC-8.2.0-2.31.1

samtools index ${1}.Assembly_Puffin_NUMT.realigned.bam ;
samtools view -b -h -o ${1}.Assembly_Puffin_MT.bam ${1}.Assembly_Puffin_NUMT.realigned.bam scaffold1766 ;
samtools sort ${1}.Assembly_Puffin_MT.bam > ${1}.Assembly_Puffin_MT.sorted.bam ;
samtools index ${1}.Assembly_Puffin_MT.sorted.bam ;
rm ${1}.Assembly_Puffin_MT.bam
rm ${1}.Assembly_Puffin_NUMT.realigned.bam*
rm ${1}.Assembly_Puffin_NUMT.realigned.bai*

```

Call SNPs with GATK4 - Haplotype Caller

```bash

mkdir GATK_Data
mkdir GATK_Log

for f in $(cat Modern.samples | grep -v "RAZ_Ancestral" | grep -v "IOM001") ; 
do sbatch mtDNA_SNPs_pt1.sh ${f} ; 
done

```

mtDNA_SNPs_pt1.sh

```bash

module load GATK/4.1.4.0-GCCcore-8.3.0-Java-1.8

gatk --java-options "-Xmx8g" HaplotypeCaller \
-VS STRICT \
-R Assembly_Puffin_NU.MT.fasta \
-I ${1}.Assembly_Puffin_MT.sorted.bam \
-ploidy 1 \
-ERC GVCF \
-L scaffold1766 \
-O ${1}.g.vcf.gz 2> Haplotype_caller.${1}.out

mv Haplotype_caller.${1}.out GATK_Log/ ;
mv ${1}.g.vcf.gz GATK_Data/ ;
mv ${1}.g.vcf.gz.tbi GATK_Data/ ;

```

```bash

cd GATK_Data/

mkdir tmp

ls *g.vcf.gz > gvcf.list

sbatch mtDNA_SNPs_pt2.sh 
```

mtDNA_SNPs_pt2.sh

```bash

module load GATK/4.1.4.0-GCCcore-8.3.0-Java-1.8

gatk --java-options "-Xmx25g" CombineGVCFs \
-R Assembly_Puffin_NU.MT.fasta \
-L scaffold1766 \
-V gvcf.list \
-O cohort.mtDNA.g.vcf.gz \
--tmp-dir=tmp 2> Puffin_MT_Combining_Log.out

gatk --java-options "-Xmx25g" GenotypeGVCFs \
-R Assembly_Puffin_NU.MT.fasta \
-L scaffold1766 \
-V cohort.mtDNA.g.vcf.gz \
-O Puffin_MT_SnpsAndIndels.raw.vcf.gz \
--tmp-dir=tmp 2> Puffin_MT_Genotyping_Log.out

mv Puffin_MT_*_Log.out ../GATK_Log/

```
