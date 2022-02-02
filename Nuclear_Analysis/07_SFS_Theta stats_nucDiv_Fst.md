
```bash

cp /XX/XX/Assembly_Puffin_NU.MT.fasta .
cp /XX/XX/ChromosomeList .

grep -v "Z" ChromosomeList | grep -v "unplaced" > ChromosomeList.edit

# make neccessary pop files
ls ../*.bam | grep -v "IOM001" | grep -v "RAZ" > bams_good

mkdir temp
cd temp
ls ../../*.bam | grep -v "IOM001" | grep -v "RAZ" > bams_good
mv bams_good Global.bamlist
cp Global.bamlist ../
cd ../
rm -r temp

#running a global ANGSD run for all sites to get sites info
# make result directory
mkdir /cluster/work/users/oliverke/Thule/ANGSD/SFS/Results
#run Puffin_ANGSD1_SFS.sh
cat ChromosomeList.edit | while read i ;
do echo ${i} ;
sbatch Puffin_ANGSD_SFS1.sh ${i} ;
done

```

Puffin_ANGSD_SFS1.sh

```bash

mkdir temp.${1}
cp Assembly_Puffin_NU.MT.fasta temp.${1}/
cp Global.bamlist temp.${1}/

cd temp.${1}/

module load SAMtools/1.9-GCC-8.2.0-2.31.1
samtools faidx Assembly_Puffin_NU.MT.fasta

module purge
module load angsd/0.931-GCC-8.2.0-2.31.1

Chromosome=$(echo ${1} | awk -F ':' '{print $1}' | awk -F '_' '{print $2"_"$3}')

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -dosnpstat 1 -C 50 -baq 2 -checkBamHeaders 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setmaxDepth 670 -setminDepth 370"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 3 -doPost 2"

Threads="16"

angsd -b Global.bamlist -GL 1 $FILTERS $TODO -P ${Threads} -ref Assembly_Puffin_NU.MT.fasta -r ${1} -out PuffinAngsd_SFS_${Chromosome}

mv PuffinAngsd_SFS_${Chromosome}.snpStat.gz ../Results/

cd ../
rm -r temp.${1}/
```

```bash

# merging all sites
zcat PuffinAngsd_SFS_chromosome_1.snpStat.gz > sfsSites.snpStat

sbatch Concat.sh #took 1 hour
```

Concat.sh

```bash

for i in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 2 3 4 5 6 7 8 9 ;
do zcat PuffinAngsd_SFS_chromosome_${i}.snpStat.gz | tail -n +2 >> sfsSites.snpStat ;
done ;

```

```bash

tail -n+2 sfsSites.snpStat | wc -l 
gzip sfsSites.snpStat

# filtering out sites where heterozygote counts comprise more than 50% of all counts (likely lumped paralogs)
zcat sfsSites.snpStat.gz | awk '($3+$4+$5+$6)>0' | awk '($12+$13+$14+$15)/($3+$4+$5+$6)<0.5' | cut -f 1,2 > sites2do
wc -l sites2do

cp sites2do ../

cd ..

# make neccessary pop files
for i in BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU ;
do grep ${i} bams_good > ${i}.bamlist ;
done

#Thu_Big
touch THU_Large.bamlist
grep "THU008" bams_good >> THU_Large.bamlist
grep "THU010" bams_good >> THU_Large.bamlist
grep "THU016" bams_good >> THU_Large.bamlist

#Thu_Small
touch THU_Small.bamlist
grep "THU001" bams_good >> THU_Small.bamlist
grep "THU002" bams_good >> THU_Small.bamlist
grep "THU009" bams_good >> THU_Small.bamlist

touch CAN.bamlist
for i in GAN GUL ;
do grep ${i} bams_good >> CAN.bamlist ;
done

touch MAI.bamlist
for i in BRE FAR GRI HOR PAP ROS WES ;
do grep ${i} bams_good >> MAI.bamlist ;
done

mkdir Results2
cp Results/sites2do Results2/

cd Results2/

ls /XX/XX/*.bam | grep -v "IOM001" | grep -v "RAZ" > bams_good

#colonies
for i in BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU;
do grep ${i} bams_good > ${i}.bamlist ;
done
#Thu_Big
touch THU_Large.bamlist
grep "THU008" bams_good >> THU_Large.bamlist
grep "THU010" bams_good >> THU_Large.bamlist
grep "THU016" bams_good >> THU_Large.bamlist
#Thu_Small
touch THU_Small.bamlist
grep "THU001" bams_good >> THU_Small.bamlist
grep "THU002" bams_good >> THU_Small.bamlist
grep "THU009" bams_good >> THU_Small.bamlist

touch CAN.bamlist
for i in GAN GUL ;
do grep ${i} bams_good >> CAN.bamlist ;
done

touch MAI.bamlist
for i in BRE FAR GRI HOR PAP ROS WES ;
do grep ${i} bams_good >> MAI.bamlist ;
done

cd ..

# Estimate SFS from sites2do (pop-wise)
# do the same with alt script - more juice for chr 1 2 3
### Do this in batches - So
ls *.bamlist | head -5 | while read line ;
do echo ${line} ;
sbatch SFS_Angsd.sh ${line} ;
done
ls *.bamlist | head -5 | while read line ;
do echo ${line} ;
sbatch SFS_AngsdAlt.sh ${line} ;
done

ls *.bamlist | head -10 | tail -5 | while read line ;
do echo ${line} ;
sbatch SFS_Angsd.sh ${line} ;
done
ls *.bamlist | head -10 | tail -5 |while read line ;
do echo ${line} ;
sbatch SFS_AngsdAlt.sh ${line} ;
done

ls *.bamlist | tail -8 | while read line ;
do echo ${line} ;
sbatch SFS_Angsd.sh ${line} ;
done
ls *.bamlist | tail -8 | while read line ;
do echo ${line} ;
sbatch SFS_AngsdAlt.sh ${line} ;
done
```

SFS_Angsd.sh ( SFS_AngsdAlt.sh with 20G RAM)

```bash


mkdir temp.${1}_${SLURM_ARRAY_TASK_ID}

cp -r Assembly_Puffin_NU.MT.fasta temp.${1}_${SLURM_ARRAY_TASK_ID}/ #reference genome

cp ChromosomeList.edit temp.${1}_${SLURM_ARRAY_TASK_ID}/ #list of chromosomes
cp Results2/${1} temp.${1}_${SLURM_ARRAY_TASK_ID}/
cd temp.${1}_${SLURM_ARRAY_TASK_ID}/

# Run command:

module load SAMtools/1.9-GCC-8.2.0-2.31.1
samtools faidx Assembly_Puffin_NU.MT.fasta

module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="1"
GENOME_ANC=Assembly_Puffin_NU.MT.fasta
GENOME_REF=Assembly_Puffin_NU.MT.fasta
Population=$(echo ${1} | awk -F '.' '{print $1}')

TODO="-doSaf 1"

Chromosome=$(cat ChromosomeList.edit | awk -F ':' '{print $1}' | grep -w Scaffolds_chromosome_${SLURM_ARRAY_TASK_ID} )
ChromosomeEdit=$(echo ${Chromosome}":")
grep -w Scaffolds_chromosome_${SLURM_ARRAY_TASK_ID} ../sites2do > sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do

sleep 120

angsd sites index sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do

sleep 120 

angsd -r ${ChromosomeEdit} -sites sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do -b ${1} -GL 1 -P ${Threads} -anc ${GENOME_ANC} -ref ${GENOME_REF} $TODO -out ${Population}_Chromosome${SLURM_ARRAY_TASK_ID}

mv sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do ${Population}.sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do
cp -r ${Population}* ../Results2/

cd ..
rm -r temp.${1}_${SLURM_ARRAY_TASK_ID}/

```

```bash

cd Results2/

for i in BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU THU_Large THU_Small CAN MAI;
do sbatch SFS_Calc1.sh ${i};
done

```

SFS_Calc1.sh

```bash

module load angsd/0.931-GCC-8.2.0-2.31.1
Threads="16"
Population=$(echo ${1})

realSFS cat ${Population}_Chromosome*.saf.idx -outnames ${Population} 
realSFS ${Population}.saf.idx -P ${Threads} -fold 1 > ${Population}.sfs 

rm ${Population}_Chromosome*.saf.gz
rm ${Population}_Chromosome*.saf.idx
rm ${Population}_Chromosome*.saf.pos.gz
rm ${Population}_Chromosome*.arg 
rm ${Population}.sites_Chromosome*_2do

```

THEN Theta

```bash
cd Results2

ls /XX/XX/*.bam | grep -v "IOM001" | grep -v "RAZ_Anc" > bams_good

for i in BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU ;
do grep ${i} bams_good > ${i}.theta.bamlist ;
done

#Thu_Big
touch THU_Large.theta.bamlist
grep "THU008" bams_good >> THU_Large.theta.bamlist
grep "THU010" bams_good >> THU_Large.theta.bamlist
grep "THU016" bams_good >> THU_Large.theta.bamlist
#Thu_Small
touch THU_Small.theta.bamlist
grep "THU001" bams_good >> THU_Small.theta.bamlist
grep "THU002" bams_good >> THU_Small.theta.bamlist
grep "THU009" bams_good >> THU_Small.theta.bamlist

for i in BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU THU_Large THU_Small ;
do sbatch SFS_Theta.sh ${i} ; #time of 2-5hrs
done

```

SFS_Theta.sh 

```bash

Population=$(echo ${1})

mkdir temp.${Population}

cp ../Results2/${Population}.sfs temp.${Population}/
cp ../Results2/${Population}.saf.idx temp.${Population}/
cp ../Results2/${Population}.saf.gz temp.${Population}/
cp ../Results2/${Population}.saf.pos.gz temp.${Population}/

cd temp.${Population}/

# Run command:
module purge
module load angsd/0.931-GCC-8.2.0-2.31.1
Threads="16"

realSFS saf2theta ${Population}.saf.idx -outname ${Population} -sfs ${Population}.sfs -fold 1 -P ${Threads}

# Sitewise Log Theta and nucDiversity globally/population wise
#thetaStat print ${Population}.thetas.idx > ${Population}.LogscaleTheta.persite.txt

#calculate Tajimas D globally/population wise per chromosome
thetaStat do_stat ${Population}.thetas.idx

#calculate Tajimas D globally/population wise per chromosome - SLIDING WINDOW
#thetaStat do_stat ${Population}.thetas.idx -win 100000 -step 100000 -outnames ${Population}.thetasWindow

rm ${Population}.sfs
rm ${Population}.saf.idx
rm ${Population}.saf.gz
rm ${Population}.saf.pos.gz

cp * ../

cd ../
rm -r temp.${Population}

```



```bash

Then R script for plotting

```

Fst

```bash

cd Results_Fst/

POPS="BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU THU_Small THU_Large MAI CAN "

LIST=$POPS
dequeue_from_list() {
  shift;
  LIST=$@
}
for pop1 in $POPS ; 
do dequeue_from_list $LIST ;
for pop2 in $LIST ; 
do echo ${pop1}" "${pop2} ;
sbatch SFS_Fst.sh ${pop1} ${pop2} ;
done ;
done ;

```

SFS_Fst.sh

```bash

mkdir temp.${1}.${2}/
cp ../Results2/${1}.saf.idx temp.${1}.${2}/
cp ../Results2/${1}.saf.gz temp.${1}.${2}/
cp ../Results2/${1}.saf.pos.gz temp.${1}.${2}/
cp ../Results2/${2}.saf.idx temp.${1}.${2}/
cp ../Results2/${2}.saf.gz temp.${1}.${2}/
cp ../Results2/${2}.saf.pos.gz temp.${1}.${2}/
cd temp.${1}.${2}/

module purge
module load angsd/0.931-GCC-8.2.0-2.31.1
Threads="1"

#2dsfs - can also be used as prior (e.g. Fst)
realSFS ${1}.saf.idx ${2}.saf.idx -P ${Threads} -fold 1 > ${1}.${2}.ml

# get Fst
realSFS fst index ${1}.saf.idx ${2}.saf.idx -sfs ${1}.${2}.ml -whichFst 1 -fstout ${1}.${2} -P ${Threads} -fold 1

# get the global estimate
realSFS fst stats ${1}.${2}.fst.idx -P ${Threads} -fold 1 > ${1}.${2}.Fst

rm ${1}.saf.idx
rm ${1}.saf.gz
rm ${1}.saf.pos.gz
rm ${2}.saf.idx
rm ${2}.saf.gz
rm ${2}.saf.pos.gz

cp * ../
cd ../
rm -r temp.${1}.${2}/
```

THEN

```bash

#Summary 

touch Colony.Fst.txt
POPS="BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU THU_Small THU_Large MAI CAN"

LIST=$POPS
dequeue_from_list() {
  shift;
  LIST=$@
}
for pop1 in $POPS ; 
do dequeue_from_list $LIST ;
for pop2 in $LIST ; 
do echo ${pop1}" "${pop2} ;
Fst=$(awk '{print $2}' ${pop1}.${pop2}.Fst)
echo -e ${pop1}"\t"${pop2}"\t"${Fst} >> Colony.Fst.txt
done ;
done ; 

```
