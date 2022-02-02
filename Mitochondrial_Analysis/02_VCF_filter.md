
Prepping Raw VCF

```bash

#Total Variants (without Outgroup)
gunzip Puffin_MT_SnpsAndIndels.raw.vcf.gz
module purge
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
module load BCFtools/1.9-intel-2018b
vcftools --vcf Puffin_MT_SnpsAndIndels.raw.vcf --remove-indels --out Puffin_MT_SnpsAndIndels.raw

```

Analysis of VCF

```bash
# More filtering
bcftools filter --SnpGap 10 -e 'QD < 2.0 || MQ < 40 || FS > 60.0 || SOR > 3 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' -o MT_SnpsAndIndels.filtered1.vcf Puffin_MT_SnpsAndIndels.raw.vcf

# More filtering on Depth and genotypeQuality
bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<15' -o MT_SnpsAndIndels.filtered2.vcf MT_SnpsAndIndels.filtered1.vcf

# Filtering on monomorphic SNPs and remove indels - only keep biallelic
bcftools filter -e 'AC==0 || AC==AN' MT_SnpsAndIndels.filtered2.vcf | bcftools view -m2 -M2 -v snps -o MT_Snps.filtered1.vcf
bcftools view -H MT_Snps.filtered1.vcf | wc -l

# check samples with missing genotypes and sites with many missing genotypes -> remove unusual samples(25%)/sites(20%, 0.8)
vcftools --vcf MT_Snps.filtered1.vcf --missing-indv --out MT_Filt_SNPs_stats
cat MT_Filt_SNPs_stats.imiss | sort -rk5 | head

vcftools --vcf MT_Snps.filtered1.vcf --missing-site --out MT_Filt_SNPs_stats
cat MT_Filt_SNPs_stats.lmiss | sort -rk6 | head 

# remove sites with missing data
vcftools --vcf MT_Snps.filtered1.vcf --max-missing 1.0 --recode --recode-INFO-all --out MT_Snps.filtered2.vcf
mv MT_Snps.filtered2.vcf.recode.vcf MT_Snps.filtered2.vcf

# no. alt allele in only 1 sample (singletons)
vcftools --vcf MT_Snps.filtered2.vcf --singletons --out MT_Filt_SNPs_stats
cat MT_Filt_SNPs_stats.singletons | tail -n+2 | grep -w "S" | wc -l

cat MT_Filt_SNPs_stats.singletons | tail -n+2 | awk '{print $5}' | sort | uniq -c | sort -rnk1 | head

# sample wise stats
cat Modern.samples | grep -vw "RAZ_Ancestral" | grep -v "IOM001" > samples.list
bcftools stats -S samples.list MT_Snps.filtered2.vcf > NoThu.MT.SNPs.stats.persample.txt

# Final
cp MT_Snps.filtered2.vcf PuffinMT_noOut_SNPs.filtered.vcf

bgzip -c -i PuffinMT_noOut_SNPs.filtered.vcf > PuffinMT_noOut_SNPs.filtered.vcf.gz
bcftools index PuffinMT_noOut_SNPs.filtered.vcf.gz

cat samples.list | while read line ;
do sample=$(echo ${line} | awk -F '.' '{print $1}') ;
echo ${sample} ;
bcftools consensus -s ${sample} -H 1 -M N -f Assembly_Puffin_mitogenome.final.fasta PuffinMT_noOut_SNPs.filtered.vcf.gz > ${sample}.CompleteMitoRef.SNPs.fasta ;
done ;
touch CompleteMitoRef.SNPs.AllPuffin.fasta
cat samples.list | while read line ;
do sample=$(echo ${line} | awk -F '.' '{print $1}') ;
cat ${sample}.CompleteMitoRef.SNPs.fasta | grep -v "^>" - | awk -v VAR=${sample} 'BEGIN { ORS=""; print ">"VAR"\n" } { print $0"\n"}' >> CompleteMitoRef.SNPs.AllPuffin.fasta ;
done
rm *.CompleteMitoRef.SNPs.fasta

```