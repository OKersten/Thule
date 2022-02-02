
```bash

cp /XX/XX/Razorbill.NU.MT.fa.gz .
cp /XX/XX/ChromosomeList .
cp /XX/XX/Assembly_Puffin_NU.MT.fasta .
cp /XX/XX/SFS/sites2do .
wc -l sites2do

touch bams_good ;
for i in SPI BJO GAN GUL BRE FAR GRI HOR PAP ROS WES IOM THU ;
do echo ${i} ;
ls /XX/XX/*.bam | grep -v "IOM001" | grep ${i} >> bams_good ;
done
ls /XX/XX/*.bam | grep "THU008" >> bams_good ;
ls /XX/XX/*.bam | grep "THU010" >> bams_good ;
ls /XX/XX/*.bam | grep "THU016" >> bams_good ;
ls /XX/XX/*.bam | grep "THU001" >> bams_good ;
ls /XX/XX/*.bam | grep "THU002" >> bams_good ;
ls /XX/XX/*.bam | grep "THU009" >> bams_good ;

#make pop sizes
echo -e "6\n6\n12\n42\n5\n6\n3\n3\n1" > sizeFile.size

#make pop names
touch popNames.name
for i in SPI BJO CAN MAI IOM THU ;
do echo ${i} >> popNames.name ;
done ;
echo "THU_Large" >> popNames.name
echo "THU_Small" >> popNames.name
echo "RAZ" >> popNames.name

#edit chromsome file
cat ChromosomeList | grep -v "unplaced" | grep -v "Z" > ChromosomeList.edit

###windows per chromosome
#ROH
tail -n+2 /XX/XX/AllSamples_ROH_Length > RoHLength
module purge
module load R/4.1.0-foss-2021a
R
df <- read.table("RoHLength")
median(df$V4) 


module load BEDTools/2.30.0-GCC-10.2.0
module load SAMtools/1.11-GCC-10.2.0
samtools faidx Assembly_Puffin_NU.MT.fasta 
cat Assembly_Puffin_NU.MT.fasta.fai | head -25 > Assembly_Puffin_NU.MT.bed
bedtools makewindows -g Assembly_Puffin_NU.MT.bed -w 150000 > Assembly_Puffin.150K.windows.bed

cat Assembly_Puffin.150K.windows.bed | awk '{print $1":"$2"-"$3-1}' > Assembly_Puffin.150K.windows.bed.edit

mkdir Results

cat ChromosomeList.edit | while read line ;
do echo ${line} ;
sbatch ABBABABA.sh ${line} ;
done 

```

ABBABABA.sh

```bash


mkdir temp.${1}_${SLURM_ARRAY_TASK_ID}

cd temp.${1}_${SLURM_ARRAY_TASK_ID}

touch bams_good ;
for i in SPI BJO GAN GUL BRE FAR GRI HOR PAP ROS WES IOM THU ;
do echo ${i} ;
ls /XX/XX/*.bam | grep -v "IOM001" | grep ${i} >> bams_good ;
done
ls /XX/XX/*.bam | grep "THU008" >> bams_good ;
ls /XX/XX/*.bam | grep "THU010" >> bams_good ;
ls /XX/XX/*.bam | grep "THU016" >> bams_good ;
ls /XX/XX/*.bam | grep "THU001" >> bams_good ;
ls /XX/XX/*.bam | grep "THU002" >> bams_good ;
ls /XX/XX/*.bam | grep "THU009" >> bams_good ;

cp /XX/XX/Razorbill.NU.MT.fa.gz .
cp /XX/XX/sizeFile.size .
cp /XX/XX/Assembly_Puffin_NU.MT.fasta .
cp /XX/XX/Assembly_Puffin.150K.windows.bed.edit .
cp /XX/XX/sites2do .

module purge 
module load SAMtools/1.9-GCC-8.2.0-2.31.1
gunzip Razorbill.NU.MT.fa.gz
samtools faidx Razorbill.NU.MT.fa
samtools faidx Assembly_Puffin_NU.MT.fasta

module purge
module load angsd/0.931-GCC-8.2.0-2.31.1

region=$(echo "${1}:")
grep ${region} Assembly_Puffin.150K.windows.bed.edit > region2do

grep -w Scaffolds_${1} sites2do > sites_${1}_2do

#run angsd
echo "Array number ${SLURM_ARRAY_TASK_ID}"
awk -v val1=${SLURM_ARRAY_TASK_ID} 'NR % 4 == val1' region2do | while read line ;
do echo "Performing ABBABABA2 on ${line}" ;
line_edit=$(echo ${line} | sed 's/:/_/g') ;
Start=$(echo ${line} | awk -F ':' '{print $2}' | awk -F '-' '{print $1}') ;
End=$(echo ${line} | awk -F ':' '{print $2}' | awk -F '-' '{print $2}') ;
cat sites_${1}_2do | awk -v val1=${Start} -v val2=${End} '$2 >= val1 && $2 <= val2 {print $1"\t"$2}' > sites_${line_edit}_2do ;
sleep 60 ;
angsd sites index sites_${line_edit}_2do ;
angsd -doAbbababa2 1 -doCounts 1 -bam bams_good -blockSize 1500 -sample 1 -sizeFile sizeFile.size -out Abbababa.150K.RoH.${line_edit} -sites sites_${line_edit}_2do -ref Assembly_Puffin_NU.MT.fasta -anc Razorbill.NU.MT.fa -r ${line} -useLast 0 -C 50 -baq 2 -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -p 1 ;
mv Abbababa.150K.RoH.${line_edit}* ../Results/ ;
done ;
cd ../

rm -r temp.${1}_${SLURM_ARRAY_TASK_ID}

```

Then

```bash

#estimate Z score etc
mkdir Zscore
nano estAvgError.R # from https://github.com/ANGSD/angsd/blob/master/R/estAvgError.R

cat ChromosomeList.edit | while read line ;
do echo ${line} ;
mkdir temp.${line} ;
sbatch estAveErrorABBABABA.sh ${line} ;
done 

```

estAveErrorABBABABA.sh

```bash


cd temp.${1}

cp /XX/XX/Results/Abbababa.150K.RoH.Scaffolds_${1}_* .
cp /XX/XX/sizeFile.size .
cp /XX/XX/popNames.name .
cp /XX/XX/Assembly_Puffin.150K.windows.bed.edit .
cp /XX/XX/estAvgError.R .

module purge
module load R/3.6.0-intel-2019a

region=$(echo "${1}:")
grep ${region} Assembly_Puffin.150K.windows.bed.edit > region2do

#run angsd
cat region2do | while read line ;
do echo "Performing error estim. on ${line}" ;
line_edit=$(echo ${line} | sed 's/:/_/g') ;
Count=$(wc -l Abbababa.150K.RoH.${line_edit}.abbababa2 | awk '{print $1}')
if [ ${Count} -gt 169 ]
then
	Rscript estAvgError.R angsdFile="Abbababa.150K.RoH.${line_edit}" out="EstError.Abbababa.150K.RoH.${line_edit}" sizeFile=sizeFile.size nameFile=popNames.name ;
	mv EstError.Abbababa.150K.RoH.${line_edit}* ../Zscore/ 
else 
	echo "Not enough blocks in ${line_edit}" 
fi ;
done ;
cd ../
rm -r temp.${1}

```

Finally

```bash

rm slurm*
cd Zscore
touch EstError.Abbababa.150K.RoH.All.Observed.txt
echo -e Chromosome"\t"Region > partOne
head -1 EstError.Abbababa.150K.RoH.Scaffolds_chromosome_9_7050000-7199999.Observed.txt | paste partOne - >> EstError.Abbababa.150K.RoH.All.Observed.txt
rm partOne
for i in {1..25} ;
do echo $i ;
ls EstError.Abbababa.150K.RoH.Scaffolds_chromosome_${i}_*.Observed.txt | while read file ;
do echo ${file} ;
Region=$(echo ${file} | awk -F '.' '{print $5}' | awk -F '_' '{print $4}') ;
echo ${Region} ;
cat ${file} | tail -n+2 | while read line ;
do echo -e ${i}"\t"${Region}"\t"${line} >> EstError.Abbababa.150K.RoH.All.Observed.txt ;
done ;
done ;
done ;


```
