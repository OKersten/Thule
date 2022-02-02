
```bash

cp /XX/XX/Razorbill.NU.MT.fa.gz .
cp /XX/XX/ChromosomeList .
cp /XX/XX/Assembly_Puffin_NU.MT.fasta .

#make bam list
ls /XX/XX/*.bam | grep -v "IOM001" | grep -v "RAZ_Anc" | grep -v "THU" > bams_good
ls /XX/XX/*.bam | grep "THU008" >> bams_good
ls /XX/XX/*.bam | grep "THU010" >> bams_good
ls /XX/XX/*.bam | grep "THU016" >> bams_good
ls /XX/XX/*.bam | grep "THU001" >> bams_good
ls /XX/XX/*.bam | grep "THU002" >> bams_good
ls /XX/XX/*.bam | grep "THU009" >> bams_good

#make pop sizes
echo -e "6\n6\n6\n6\n6\n6\n6\n5\n6\n6\n6\n6\n3\n3\n1" > sizeFile.size
#make pop names
awk -F '/' '{print $2}' bams_good | awk -F '.' '{print $1}' | awk -F '0' '{print $1}' | sort | uniq > popNames.name
echo "THU_Large" >> popNames.name
echo "THU_Small" >> popNames.name
echo "RAZ" >> popNames.name

#edit chromsome file
cat ChromosomeList | grep -v "unplaced" | grep -v "Z" > ChromosomeList.edit

#call script
cat ChromosomeList.edit | while read line ;
do echo ${line} ;
mkdir temp.${line} ;
sbatch ABBABABA.sh ${line} ;
done 
```

ABBABABA.sh

```bash

cd temp.${1}

ls ../../*.bam | grep -v "IOM001" | grep -v "RAZ_Anc" | grep -v "THU" > bams_good
ls ../../*.bam | grep "THU008" >> bams_good
ls ../../*.bam | grep "THU010" >> bams_good
ls ../../*.bam | grep "THU016" >> bams_good
ls ../../*.bam | grep "THU001" >> bams_good
ls ../../*.bam | grep "THU002" >> bams_good
ls ../../*.bam | grep "THU009" >> bams_good
cp ../Razorbill.NU.MT.fa.gz .
cp ../sizeFile.size .
cp ../Assembly_Puffin_NU.MT.fasta .

module purge 
module load SAMtools/1.9-GCC-8.2.0-2.31.1
gunzip Razorbill.NU.MT.fa.gz
samtools faidx Razorbill.NU.MT.fa
samtools faidx Assembly_Puffin_NU.MT.fasta

module purge
module load angsd/0.931-GCC-8.2.0-2.31.1

#Specify chromosome
Chromosome=$(echo "Scaffolds_"${1}":")

#run angsd
angsd -doAbbababa2 1 -doCounts 1 -bam bams_good -sizeFile sizeFile.size -out ${1} -ref Assembly_Puffin_NU.MT.fasta -anc Razorbill.NU.MT.fa -r ${Chromosome} -C 50 -baq 2 -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 77 -p 1

cp -r *.abbababa2 ../
cp -r *.arg ../

cd ../

rm -r temp.${1}
```

Then

```bash

cat chromosome_1.abbababa2 > AllChromosomes.Multipop.abbababa2
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 ;
do cat chromosome_${i}.abbababa2 | tail -n +2 >> AllChromosomes.Multipop.abbababa2 ;
done

#estimate Z score etc
module purge
module load R/3.6.0-intel-2019a
Rscript estAvgError.R angsdFile="AllChromosomes.Multipop" out="Result.AllChromosomes.Multipop" sizeFile=sizeFile.size nameFile=popNames.name


Rscript Dstats.R
```

