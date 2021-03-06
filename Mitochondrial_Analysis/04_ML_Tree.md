

Get Fasta's from VCF (with outgroup)

```bash

module load BCFtools/1.9-intel-2018b

bgzip -c -i PuffinMT_RazorMurreAukl_SNPs.join.vcf > PuffinMT_RazorMurreAukl_SNPs.join.vcf.gz
bcftools index PuffinMT_RazorMurreAukl_SNPs.join.vcf.gz


#Get sample fasta from  multisample vcf
cat Modern.samples | grep -v "RAZ_Ancestral" - | grep -v "IOM001" - | while read line ;
do sample=$(echo ${line} | awk -F '.' '{print $1}') ;
echo ${sample} ;
bcftools consensus -s ${sample} -H 1 -M N -f Assembly_Puffin_mitogenome.final.fasta PuffinMT_RazorMurreAukl_SNPs.join.vcf.gz > ${sample}.CompleteMitoRef.SNPs.fasta ;
done ;

for sample in Alca_torda/1-17052 Synthliboramphus_antiquus/1-16634 Synthliboramphus_wumizusume/1-16630 Aethia/1-16626 ;
do ID=$(echo ${sample} | awk -F '/' '{print $1}') ;
echo ${sample} ;
echo ${ID} ;
bcftools consensus -s ${sample} -H 1 -M N -f Assembly_Puffin_mitogenome.final.fasta PuffinMT_RazorMurreAukl_SNPs.join.vcf.gz > ${ID}.CompleteMitoRef.SNPs.fasta ;
done ;

mv Aethia.CompleteMitoRef.SNPs.fasta Aethia_cristatella.CompleteMitoRef.SNPs.fasta
```

ML Tree (All SNPs but no Ns)

```bash

# remove N’s and merge them all with new labels


module load BEDTools/2.27.1-intel-2018b
cat Modern.samples | grep -v "RAZ_Ancestral" | grep -v "IOM001" > Sample.list
echo "Alca_torda" >> Sample.list
echo "Synthliboramphus_antiquus" >> Sample.list
echo "Synthliboramphus_wumizusume" >> Sample.list
echo "Aethia_cristatella" >> Sample.list

echo -e "scaffold1766\t0\t15555" > Mitogenome.NoNs.bed
echo -e "scaffold1766\t15698\t16849" >> Mitogenome.NoNs.bed
echo -e "scaffold1766\t16859\t17084" >> Mitogenome.NoNs.bed

touch All.CompleteMitoRef.SNPs.withOut.NoNs.fasta
cat Sample.list | while read line ;
do echo ${line} ;
bedtools getfasta -fi ${line}.CompleteMitoRef.SNPs.fasta -bed Mitogenome.NoNs.bed | grep -v "^>" - | awk -v VAR=${line} 'BEGIN { ORS=""; print ">"VAR"\n" } { print $0"\n"}' >> All.CompleteMitoRef.SNPs.withOut.NoNs.fasta ;
done

# replace all (potentially) *'s with -
sed 's/*/-/g' All.CompleteMitoRef.SNPs.withOut.NoNs.fasta > All.CompleteMitoRef.SNPs.withOut.NoNs.edit.fasta

# Make tree
module load IQ-TREE/1.6.12-foss-2018b
echo "Start IQ tree"
iqtree -s All.CompleteMitoRef.SNPs.withOut.NoNs.edit.fasta -m MFP -alrt 1000 -bb 1000 -AICc -bnni -o Alca_torda,Synthliboramphus_antiquus,Synthliboramphus_wumizusume,Aethia_cristatella -nt AUTO -ntmax 8
echo "IQ Tree done"

```

Gene ML Tree with a annotation bed-file of the Puffin mitogenome (scaffold1766-2.bed)

```bash

#########################################
######### Coding Gene Tree ##############
#########################################
#1 partition 1st codon positon of ProtCodGene
#1 partition 2nd codon position of ProtCodGene
#1 Partition 3rd codon position of ProtCodGene
#1 Partition of tRNA genes
#1 Partition of rRNA genes
#1 Partition of Control Region
#1 Partition of intergenic

module load BEDTools/2.27.1-intel-2018b
module load seqtk/1.3-foss-2018b

#adjust bed for whatever region your interested 
grep -v "rrn" scaffold1766-2.bed | grep -v "trn" | head -n-3 > scaffold1766_coding.bed 
grep "rrn" scaffold1766-2.bed > scaffold1766_rRNA.bed 
grep "trn" scaffold1766-2.bed > scaffold1766_tRNA.bed
tail -4 scaffold1766-2.bed > scaffold1766_CR.bed 
#manually change to have one region span end of tRNA (15548) to end of genome
nano scaffold1766_CR.bed 
#manually change all codons by -3 (end) except Cox3 and NAD3_1 and NAD6 (change to start 3 later)
nano scaffold1766_coding.bed 

cat scaffold1766_*.bed > scaffold1766_intergenic.bed


#CDS
cat Sample.list | while read sample ;
do echo ${sample} ;
bedtools getfasta -fi ${sample}.CompleteMitoRef.SNPs.fasta -bed scaffold1766_coding.bed -name -s -fo ${sample}.CompleteRef.coding.fasta ;
done ;

#tRNA
cat Sample.list | while read sample ;
do echo ${sample} ;
bedtools getfasta -fi ${sample}.CompleteMitoRef.SNPs.fasta -bed scaffold1766_tRNA.bed -name -s -fo ${sample}.CompleteRef.trna.fasta ;
done ;

#rRNA
cat Sample.list | while read sample ;
do echo ${sample} ;
bedtools getfasta -fi ${sample}.CompleteMitoRef.SNPs.fasta -bed scaffold1766_rRNA.bed -name -s -fo ${sample}.CompleteRef.rrna.fasta ;
done ;

#CR
cat Sample.list | while read sample ;
do echo ${sample} ;
bedtools getfasta -fi ${sample}.CompleteMitoRef.SNPs.fasta -bed scaffold1766_CR.bed  -name -s -fo ${sample}.CompleteRef.cr.fasta ;
sed 's/N//g' ${sample}.CompleteRef.cr.fasta > ${sample}.CompleteRef.cr.edit.fasta ;
rm ${sample}.CompleteRef.cr.fasta ;
done ;

#intergenic
cat Sample.list | while read sample ;
do echo ${sample} ;
bedtools maskfasta -fi ${sample}.CompleteMitoRef.SNPs.fasta -bed scaffold1766_intergenic.bed -fo ${sample}.CompleteRef.inter.fasta -mc X ;
seqtk seq -l 0 ${sample}.CompleteRef.inter.fasta > ${sample}.CompleteRef.inter.edit.fasta ;
sed 's/X//g' ${sample}.CompleteRef.inter.edit.fasta > ${sample}.CompleteRef.inter.edit2.fasta ;
done

##Concat multifasta into 1 fasta (per sample) sequence and all fastas together
#CDS
cat Sample.list | while read sample ;
do echo ${sample} ;
grep -v "^>" ${sample}.CompleteRef.coding.fasta | awk -v VAR=${sample} 'BEGIN { ORS=""; print ">"VAR"\n" } { print }' > ${sample}.CompleteRef.CDS.fasta ;
seqtk seq -l 80 ${sample}.CompleteRef.CDS.fasta > ${sample}.CompleteRef.CDS.edit.fasta ;
done ;

cat *CompleteRef.CDS.edit.fasta | seqtk seq -l 0 - > All.Aligned.CDS.fasta
tail -4 Sample.list > outlier
grep ">" All.Aligned.CDS.fasta | grep -v -f outlier - | awk -F '>' '{print $2}' > nooutlier
seqtk subseq All.Aligned.CDS.fasta outlier | seqtk seq -l 80 - > All.Aligned.CDS.outlier.fasta
seqtk subseq All.Aligned.CDS.fasta nooutlier | seqtk seq -l 80 - > All.Aligned.CDS.nooutlier.fasta
cat All.Aligned.CDS.outlier.fasta All.Aligned.CDS.nooutlier.fasta | seqtk seq -l 0 - > All.Aligned.CDS.edit.fasta

#tRNA
cat Sample.list | while read sample ;
do echo ${sample} ;
grep -v "^>" ${sample}.CompleteRef.trna.fasta | awk -v VAR=${sample} 'BEGIN { ORS=""; print ">"VAR"\n" } { print }' > ${sample}.CompleteRef.trna.join.fasta ;
seqtk seq -l 80 ${sample}.CompleteRef.trna.join.fasta > ${sample}.CompleteRef.trna.join.edit.fasta ;
done ;

cat *CompleteRef.trna.join.edit.fasta | seqtk seq -l 0 - > All.Aligned.tRNA.fasta
seqtk subseq All.Aligned.tRNA.fasta outlier | seqtk seq -l 80 - > All.Aligned.trna.outlier.fasta
seqtk subseq All.Aligned.tRNA.fasta nooutlier | seqtk seq -l 80 - > All.Aligned.trna.nooutlier.fasta
cat All.Aligned.trna.outlier.fasta All.Aligned.trna.nooutlier.fasta | seqtk seq -l 0 - > All.Aligned.tRNA.edit.fasta

#rRNA
cat Sample.list | while read sample ;
do echo ${sample} ;
grep -v "^>" ${sample}.CompleteRef.rrna.fasta | awk -v VAR=${sample} 'BEGIN { ORS=""; print ">"VAR"\n" } { print }' > ${sample}.CompleteRef.rrna.join.fasta ;
seqtk seq -l 80 ${sample}.CompleteRef.rrna.join.fasta > ${sample}.CompleteRef.rrna.join.edit.fasta ;
done ;

cat *CompleteRef.rrna.join.edit.fasta | seqtk seq -l 0 - > All.Aligned.rRNA.fasta
seqtk subseq All.Aligned.rRNA.fasta outlier | seqtk seq -l 80 - > All.Aligned.rrna.outlier.fasta
seqtk subseq All.Aligned.rRNA.fasta nooutlier | seqtk seq -l 80 - > All.Aligned.rrna.nooutlier.fasta
cat All.Aligned.rrna.outlier.fasta All.Aligned.rrna.nooutlier.fasta | seqtk seq -l 0 - > All.Aligned.rRNA.edit.fasta

#CR
cat Sample.list | while read sample ;
do echo ${sample} ;
grep -v "^>" ${sample}.CompleteRef.cr.edit.fasta | awk -v VAR=${sample} 'BEGIN { ORS=""; print ">"VAR"\n" } { print }' > ${sample}.CompleteRef.ContrReg.fasta ;
seqtk seq -l 80 ${sample}.CompleteRef.ContrReg.fasta > ${sample}.CompleteRef.ContrReg.edit.fasta ;
done ;

cat *CompleteRef.ContrReg.edit.fasta | seqtk seq -l 0 - > All.Aligned.CR.fasta
seqtk subseq All.Aligned.CR.fasta outlier | seqtk seq -l 80 - > All.Aligned.cr.outlier.fasta
seqtk subseq All.Aligned.CR.fasta nooutlier | seqtk seq -l 80 - > All.Aligned.cr.nooutlier.fasta
cat All.Aligned.cr.outlier.fasta All.Aligned.cr.nooutlier.fasta | seqtk seq -l 0 - > All.Aligned.CR.edit.fasta

#Intergenic
cat Sample.list | while read sample ;
do echo ${sample} ;
grep -v "^>" ${sample}.CompleteRef.inter.edit2.fasta | awk -v VAR=${sample} 'BEGIN { ORS=""; print ">"VAR"\n" } { print }' > ${sample}.CompleteRef.inter.edit3.fasta ;
seqtk seq -l 80 ${sample}.CompleteRef.inter.edit3.fasta > ${sample}.CompleteRef.inter.edit4.fasta ;
done ;

cat *CompleteRef.inter.edit4.fasta | seqtk seq -l 0 - > All.Aligned.inter.fasta
seqtk subseq All.Aligned.inter.fasta outlier | seqtk seq -l 80 - > All.Aligned.inter.outlier.fasta
seqtk subseq All.Aligned.inter.fasta nooutlier | seqtk seq -l 80 - > All.Aligned.inter.nooutlier.fasta
cat All.Aligned.inter.outlier.fasta All.Aligned.inter.nooutlier.fasta | seqtk seq -l 0 - > All.Aligned.inter.edit.fasta

#replace all *'s with -
sed 's/*/-/g' All.Aligned.CDS.edit.fasta > All.Aligned.CDS.edit2.fasta
sed 's/*/-/g' All.Aligned.tRNA.edit.fasta > All.Aligned.tRNA.edit2.fasta
sed 's/*/-/g' All.Aligned.rRNA.edit.fasta > All.Aligned.rRNA.edit2.fasta
sed 's/*/-/g' All.Aligned.CR.edit.fasta > All.Aligned.CR.edit2.fasta
sed 's/*/-/g' All.Aligned.inter.edit.fasta > All.Aligned.inter.edit2.fasta

#from each of them
#remove SPI015 and SPI002 (identical to BJO001)
#remove HOR008 (identical to BRE004)
#remove PAP006 (identical to PAP003)
#remove ROS001 and PAP004 and ROS006 (identical to BRE001)

cat Sample.list | grep -v SPI015 | grep -v SPI002 | grep -v HOR008 | grep -v PAP006 | grep -v ROS006 | grep -v ROS001 | grep -v PAP004 > Good.samples
seqtk subseq All.Aligned.CDS.edit2.fasta Good.samples > All.Aligned.CDS.edit3.fasta
seqtk subseq All.Aligned.tRNA.edit2.fasta Good.samples > All.Aligned.tRNA.edit3.fasta
seqtk subseq All.Aligned.rRNA.edit2.fasta Good.samples > All.Aligned.rRNA.edit3.fasta
seqtk subseq All.Aligned.CR.edit2.fasta Good.samples > All.Aligned.CR.edit3.fasta
seqtk subseq All.Aligned.inter.edit2.fasta Good.samples > All.Aligned.inter.edit3.fasta

#NEXUS
#check CDS sequence length
module purge
module load SAMtools/1.9-intel-2018b
samtools faidx All.Aligned.CDS.edit3.fasta
cut -f1-2 All.Aligned.CDS.edit3.fasta.fai | head
#then
nano Puffin.mt.part.nex

#nexus
begin sets;
    charset part1 = All.Aligned.CDS.edit3.fasta: 1-11374\3;
    charset part2 = All.Aligned.CDS.edit3.fasta: 2-11375\3;
    charset part3 = All.Aligned.CDS.edit3.fasta: 3-11376\3;
    charset part4 = All.Aligned.rRNA.edit3.fasta: *;
    charset part5 = All.Aligned.tRNA.edit3.fasta: *;
    charset part6 = All.Aligned.CR.edit3.fasta: *;
    charset part7 = All.Aligned.inter.edit3.fasta: *;
end;

########### IQ tree ############

module purge
module load IQ-TREE/1.6.12-foss-2018b
echo "Start IQ tree"
iqtree -spp Puffin.mt.part.nex -m MFP+MERGE -bnni -AICc -alrt 1000 -bb 1000 -bsam GENESITE -o Alca_torda,Synthliboramphus_antiquus,Synthliboramphus_wumizusume,Aethia_cristatella -nt AUTO -ntmax 8

echo "IQ Tree done"

# alter treefile to include excluded taxa with branch length 0
make new copy -> Puffin.mt.part.nex.edit.treefile
# SPI015 and SPI002 identical to BJO001 -> replace BJO001 with ((SPI002:0.0000000000,SPI015:0.0000000000)100/100:0.0000000000,BJO001:0.0000000000)100/100
# HOR008 identical to BRE004 -> replace BRE004 with (HOR008:0.0000000000,BRE004:0.0000000000)100/100
# PAP006 identical to PAP003 -> replace PAP003 with (PAP006:0.0000000000,PAP003:0.0000000000)100/100
# ROS001 and PAP004 and ROS006 identical to BRE001 -> replace BRE001 with (((ROS001:0.0000000000,ROS006:0.0000000000)100/100:0.0000000000,PAP004:0.0000000000)100/100:0.0000000000,BRE001:0.0000000000)100/100

```

Plotting

```bash

#For plotting, remove the aLHRT test
sed 's/[0-9][0-9]\.[0-9]\///g' Puffin.mt.part.nex.edit.treefile | sed 's/[0-9][0-9][0-9]\///g' | sed 's/[0-9][0-9]\///g' | sed 's/[0-9]\///g' > Puffin.mt.part.nex.edit2.treefile

#use R afterwards

```