Align mitochondrial fastas with Muscle, then call Snps with snp-sites

```
# Get Outgroup species mitogenome from NCBI
# Accession: CM018102.1 - Alcidae - Alca torda - Razorbill
# Accession: NC_007978.1 - Alcidae - Synthliboramphus antiquus - Ancient murrelet
# Accession: NC_029328.1 - Alcidae - Synthliboramphus wumizusume - Japanese murrelet
# Accession: NC_045517.1 - Alcidae - Aethia cristatella - Crested auklet


# make sure outgroup fasta start at the same nucl. as Puffin mitogenome (only ancient murrelet)
/cluster/home/oliverke/Programs/seqkit restart -i 3081 Synthliboramphus.antiquus.fasta > Synthliboramphus.antiquus.restart1.fasta

# align with Puffin mitogenome
module load MUSCLE/3.8.31-foss-2018b

muscle -profile -in1 Assembly_Puffin_mitogenome.final.fasta -in2 Alca.torda.fasta -out Puffin_razorbill.fasta
muscle -in Puffin_razorbill.fasta -out Puffin_razorbill.refined.fasta -refine

muscle -profile -in1 Puffin_razorbill.refined.fasta -in2 Synthliboramphus.antiquus.restart1.fasta -out Puffin_razorbill_murrelet1.fasta
muscle -in Puffin_razorbill_murrelet1.fasta -out Puffin_razorbill_murrelet1.refined.fasta -refine

muscle -profile -in1 Puffin_razorbill_murrelet1.refined.fasta -in2 Synthliboramphus.wumizusume.fasta -out Puffin_razorbill_murrelet2.fasta
muscle -in Puffin_razorbill_murrelet2.fasta -out Puffin_razorbill_murrelet2.refined.fasta -refine

muscle -profile -in1 Puffin_razorbill_murrelet2.refined.fasta -in2 Aethia.cristatella.fasta -out Puffin_razorbill_murrelet2_auklet.fasta
muscle -in Puffin_razorbill_murrelet2_auklet.fasta -out Puffin_razorbill_murrelet2_auklet.refined.fasta -refine


# check manuallyin jalview - remove indels on Puffin reference, save as Puffin_razorbill_murrelet2_auklet.refined.manual.fasta
# extract snp sites
snp-sites -v -p -o Puffin_razorbill.murrelet.auklet.snps Puffin_razorbill_murrelet2_auklet.refined.manual.fasta

```

Merge SNPs of Puffin with SNPs Puffin-Razorbill-Murrelet

```bash

module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
module load BCFtools/1.9-intel-2018b

vcftools --remove-indv scaffold1766/1-17084 --recode --vcf Puffin_razorbill.murrelet.auklet.snps.edit.vcf --out Puffin_razorbill.murrelet.auklet.snps.edit2.vcf
bgzip -c -i Puffin_razorbill.murrelet.auklet.snps.edit2.vcf.recode.vcf > Puffin_razorbill.murrelet.auklet.snps.edit.vcf.gz
bcftools index Puffin_razorbill.murrelet.auklet.snps.edit.vcf.gz
echo "1 scaffold1766" > chr_name_conv.txt
bcftools annotate --rename-chrs chr_name_conv.txt Puffin_razorbill.murrelet.auklet.snps.edit.vcf.gz | bgzip > Puffin_razorbill.murrelet.auklet.snps.edit2.vcf.gz
bcftools index Puffin_razorbill.murrelet.auklet.snps.edit2.vcf.gz

bgzip -c -i PuffinMT_noOut_SNPs.filtered.vcf > PuffinMT_noOut_SNPs.filtered.vcf.gz
bcftools index PuffinMT_noOut_SNPs.filtered.vcf.gz
bcftools merge --missing-to-ref PuffinMT_noOut_SNPs.filtered.vcf.gz Puffin_razorbill.murrelet.auklet.snps.edit2.vcf.gz > PuffinMT_RazorMurreAukl_SNPs.join.vcf
grep -v "#" PuffinMT_RazorMurreAukl_SNPs.join.vcf | wc -l


```