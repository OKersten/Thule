```bash


zcat /XX/XX/PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz | tail -n+2 | cut -f 1 | awk -F '_' '{print $1"_"$2"_"$3"\t"$4}' > sites2do
wc -l sites2do

cp /XX/XX/ChromosomeList .
grep -v "Z" ChromosomeList | grep -v "unplaced" > ChromosomeList.edit

mkdir test
cd test
ls /XX/XX/*.Assembly_Puffin_NU.sorted.bam | grep -v "IOM001" > bams_good.edit
cp bams_good.edit  ../
cd ../
rm -r test

mkdir Results

sbatch Outgroup_ANGSD_EEMS.sh

```

Outgroup_ANGSD_EEMS.sh

```bash


mkdir temp_${SLURM_ARRAY_TASK_ID}

cp -r /XX/XX/Assembly_Puffin_NU.MT.fasta temp_${SLURM_ARRAY_TASK_ID}/ #reference genome
cp bams_good.edit temp_${SLURM_ARRAY_TASK_ID}/  
cp sites2do temp_${SLURM_ARRAY_TASK_ID}/ 

cd temp_${SLURM_ARRAY_TASK_ID}/

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -C 50 -baq 2 -checkBamHeaders 1 -doHWE 1 -skipTriallelic 1"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -doIBS 2"

ChromosomeEdit=$(echo "Scaffolds_chromosome_"${SLURM_ARRAY_TASK_ID}":")
grep -w Scaffolds_chromosome_${SLURM_ARRAY_TASK_ID} sites2do > sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do

# Run command:
module load SAMtools/1.9-GCC-8.2.0-2.31.1
samtools faidx Assembly_Puffin_NU.MT.fasta

module purge
module load angsd/0.931-GCC-8.2.0-2.31.1

sleep 120

angsd sites index sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do

sleep 120 

angsd -b bams_good.edit -GL 1 $FILTERS $TODO -ref Assembly_Puffin_NU.MT.fasta -r ${ChromosomeEdit} -sites sites_Chromosome${SLURM_ARRAY_TASK_ID}_2do -out PuffinAngsd_IBS_withOutgroup_Chromosome${SLURM_ARRAY_TASK_ID}

cp -r PuffinAngsd_IBS_withOutgroup_Chromosome${SLURM_ARRAY_TASK_ID}* ../Results/

cd ../

rm -r temp_${SLURM_ARRAY_TASK_ID}/

```

#Samples

Run 100 replicates (different seed) → best likelihood as tree

Run 100 bootstraps (different seed) → bootstrap support

#Pops

Run 100 replicates (different seed) for each migration → best likelihood to display as tree and SE

check how likelihood and variation explained change → graph m vs like or m vs. var expl. dotplo

Choose best likelihood tree for best m and do 100 bootstrap  

```bash

#No Z or unplaced
zcat PuffinAngsd_IBS_withOutgroup_Chromosome1.ibs.gz > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ibs 
for i in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 2 3 4 5 6 7 8 9 ;
do zcat PuffinAngsd_IBS_withOutgroup_Chromosome${i}.ibs.gz | tail -n +2 >> PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ibs ;
done ;
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ibs
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ibs.gz | tail -n +2 | wc -l

mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ibs.gz ../
cd ../
rm -r Results/
rm slurm*

cat bams_good.edit| awk -F '/' '{print $4}' | awk -F '.' '{print $1}' > samples
cat samples | grep -v "THU" | awk -F '0' '{print $1}' | sort | uniq > pops
echo "THU_Large" >> pops
echo "THU_Small" >> pops

#Samples
Samples=$(cat samples | tr "\n" " ")
echo -e chr" "pos" "major" "minor" "${Samples} > header
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ibs.gz | tail -n+2 | cat header - > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.edit.ibs
rm header
cat samples | while read i ;
do echo ${i} ;
ColumnNo=$(head -1 PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.edit.ibs | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
awk -v c1=$ColumnNo -F ' ' '{print $c1}' PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.edit.ibs | tail -n+2 | sed 's/0/A/g ; s/-1/B/g ; s/1/C/g; s/A/1,0/g ; s/B/0,0/g ; s/C/0,1/g' > ${i}.column ;
echo ${i} | cat - ${i}.column > ${i}.column.edit ;
rm ${i}.column ;
done
paste -d " " *.column.edit > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples ;
rm *.column.edit
rm PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.edit.ibs
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples

#Pops
cat pops | grep -v "THU" | while read line ;
do echo ${line} ;
touch pop_refallele ;
touch pop_altallele ;
grep ${line} samples > samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;
done ;

line=THU_Large;
echo ${line};
touch pop_refallele ;
touch pop_altallele ;
touch samples.pop ;
grep "THU008" samples >> samples.pop ;
grep "THU010" samples >> samples.pop ;
grep "THU016" samples >> samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;

line=THU_Small ;
echo ${line} ;
touch pop_refallele ;
touch pop_altallele ;
touch samples.pop ;
grep "THU001" samples >> samples.pop ;
grep "THU002" samples >> samples.pop ;
grep "THU009" samples >> samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;

paste -d " " PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_* >> PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops

rm PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_*
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops

#Change RAZ to random ref or alt allele when missing data
#samples
gunzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz
ColumnNo=$(head -1 PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples | tr " " "\n" | grep -w -n "RAZ_Anc" | awk -F ':' '{print $1}') ;
awk -v c1=$ColumnNo -F ' ' '{print $c1}' PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples |awk '{$0=gensub(/0,0/, 2+int(rand()*10000), "g", $0)}1' | awk '($1 <5002 && $1 >1) {print "1,0"; next}1' | awk '($1 >5001 && $1 <10002 ) {print "0,1"; next}1' > RAZ.column ; 
cut -d ' ' -f 1-53,55-78 PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples | paste -d ' ' - RAZ.column > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.edit
rm RAZ.column ;
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.edit
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples

#pops
gunzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz
ColumnNo=$(head -1 PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops | tr " " "\n" | grep -w -n "RAZ_Anc" | awk -F ':' '{print $1}') ;
awk -v c1=$ColumnNo -F ' ' '{print $c1}' PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops |awk '{$0=gensub(/0,0/, 2+int(rand()*10000), "g", $0)}1' | awk '($1 <5002 && $1 >1) {print "1,0"; next}1' | awk '($1 >5001 && $1 <10002 ) {print "0,1"; next}1' > RAZ.column ; 
cut -d ' ' -f 1-9,11-15 PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops | paste -d ' ' - RAZ.column > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.edit
rm RAZ.column ;
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.edit
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops

```

ADDING CAN AND MAI to IOM BJO SPI

```bash

line=CAN ;
echo ${line} ;
touch pop_refallele ;
touch pop_altallele ;
touch samples.pop ;
echo "GUL" >> samples.pop ;
echo "GAN" >> samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;

line=MAI ;
echo ${line} ;
touch pop_refallele ;
touch pop_altallele ;
touch samples.pop ;
echo "ROS" >> samples.pop ;
echo "HOR" >> samples.pop ;
echo "WES" >> samples.pop ;
echo "BRE" >> samples.pop ;
echo "PAP" >> samples.pop ;
echo "GRI" >> samples.pop ;
echo "FAR" >> samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;

line=THU_Large;
echo ${line};
touch pop_refallele ;
touch pop_altallele ;
touch samples.pop ;
grep "THU008" samples >> samples.pop ;
grep "THU010" samples >> samples.pop ;
grep "THU016" samples >> samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;

line=THU_Small ;
echo ${line} ;
touch pop_refallele ;
touch pop_altallele ;
touch samples.pop ;
grep "THU001" samples >> samples.pop ;
grep "THU002" samples >> samples.pop ;
grep "THU009" samples >> samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;

for line in IOM SPI BJO ;
do echo ${line} ;
touch pop_refallele ;
touch pop_altallele ;
grep ${line} samples > samples.pop ;
cat samples.pop | while read i ;
do echo ${i} ;
ColumnNo=$(zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | head -1 | tr " " "\n" | grep -w -n ${i} | awk -F ':' '{print $1}') ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $1}' > refallele${i} ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.gz | awk -v c1=$ColumnNo -F ' ' '{print $c1}' | tail -n +2 | awk -F ',' '{print $2}' > altallele${i} ;
done ;
paste refallele* >> pop_refallele ;
paste altallele* >> pop_altallele ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_refallele > refallele_sum ;
awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' pop_altallele > altallele_sum ;
paste -d "," refallele_sum altallele_sum > popsum ;
echo ${line} | cat - popsum > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_${line} ;
rm refallele* ;
rm altallele* ;
rm samples.pop ;
rm pop_refallele ;
rm pop_altallele ;
rm popsum ;
done ;

paste -d " " PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_* >> PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters

rm PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix_*
gzip PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters
##########################################
```

```bash
##########################################

mkdir sample_like_boot
for N in $(seq 1 100) ;
do echo ${N};
sbatch Treemix_Sample_Like.sh ${N};
done ;

mkdir pop_like_boot
for N in $(seq 0 10) ;
do echo ${N};
for M in {1..100} ;
do echo ${M} ;
sbatch Treemix_Pop_Mig_Like.sh ${N} ${M};
done ;
done ;

#starting them all at the same time might cause erros with SAGA ....
rm slurm*
ls | grep temp > toRepeat.txt
rm -r temp*
cat toRepeat.txt | while read line ;
do N=$(echo ${line} | awk -F '.' '{print $3}')
M=$(echo ${line} | awk -F '.' '{print $4}')
echo ${N}"_"${M}
sbatch Treemix_Pop_Mig_Like.sh ${N} ${M};
done ;

```

Treemix_Sample_Like.sh 

```bash

mkdir temp.sample.${1}

cp PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.edit.gz temp.sample.${1}/

cd temp.sample.${1}

module load TreeMix/1.13-intel-2018b

#making a random seed
R1=$(echo $RANDOM % 100 | bc)
R2=$(echo $RANDOM % 100 | bc)
S=$(bc <<<"scale=0; ${R1} * ${R2} * ${1}")

## Run treemix
treemix -i PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.edit.gz -global -noss -root RAZ_Anc -o PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.${1} -seed ${S}

mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.${1}* ../sample_like_boot/

cd ../
rm -r temp.sample.${1}

```

Treemix_Pop_Mig_Like.sh

```bash


mkdir temp.pop.${1}.${2}

cp PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.edit.gz temp.pop.${1}.${2}/

cd temp.pop.${1}.${2}/

sleep 120 

module load TreeMix/1.13-intel-2018b

#making a random seed
R1=$(echo $RANDOM % 100 | bc)
R2=$(echo $RANDOM % 100 | bc)
S=$(bc <<<"scale=0; ${R1} * ${R2} * ${1}")

## Run treemix
treemix -i PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.edit.gz -m ${1} -global -root RAZ_Anc -o PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.${2}.m${1} -seed ${S}

mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.${2}.m${1}* ../pop_like_boot/

cd ../
rm -r temp.pop.${1}.${2}

```

```bash

#samples
cd sample_like_boot
Like=$(cat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.*.llik | awk -F ':' '{print $2}' | sed 's/ //g' | sort -r -n | uniq | head -1) 
Replicate=$(grep -l -w ${Like} PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.*.llik | head -1 | awk -F '.' '{print $6}')
mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.${Replicate}.* ../
cd ..
rm sample_like_boot/*
mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.${Replicate}.* sample_like_boot/
cd sample_like_boot/
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.${Replicate}.treeout.gz | head -1 > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.main_tree ;

#pop
#check best m 
#Plot.Variation.R


cd pop_like_boot
for N in $(seq 0 10) ;
do echo ${N};
Like=$(cat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.*.m${N}.llik | awk -F ':' '{print $2}' | sed 's/ //g' | sort -r -n | uniq | head -1) ;
Replicate=$(grep -l -w ${Like} PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.*.m${N}.llik | head -1 | awk -F '.' '{print $5}') ;
mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.${Replicate}.m${N}.* ../ ;
cd .. ;
rm pop_like_boot/PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.*.m${N}.* ;
mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.${Replicate}.m${N}.* pop_like_boot/ ;
cd pop_like_boot/ ;
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.${Replicate}.m${N}.treeout.gz | head -1 > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m${N}.main_tree ;
done

#plot best m's Panel on local
for N in {0..10} ;
do echo ${N} ;
for i in vertices.gz treeout.gz modelcov.gz llik edges.gz covse.gz cov.gz ;
do mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.*.m${N}.${i} PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m${N}.${i} ;
done ;
done ; 

#use R with script plotting_funcs.R
R
prefix="PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m"
library(RColorBrewer)
library(R.utils)
source("plotting_funcs.R")
#trees
pdf(file = "TreemixPanel_All.pdf")
par(mfrow=c(4,3))
for(edge in 0:10){
plot_tree(cex=0.5,paste0(prefix,edge))
title(paste(edge,"migration(s)"))
}
dev.off()

pdf(file = "TreemixPanel_m0_1.pdf")
par(mfrow=c(2,1))
for(edge in 0:1){
plot_tree(cex=0.5,paste0(prefix,edge))
title(paste(edge,"migration(s)"))
}
dev.off()

pdf(file = "TreemixPanel_m2_5.pdf")
par(mfrow=c(3,2))
for(edge in 2:5){
plot_tree(cex=0.5,paste0(prefix,edge))
title(paste(edge,"migration(s)"))
}
dev.off()

pdf(file = "TreemixPanel_m6_10.pdf")
par(mfrow=c(3,2))
for(edge in 6:10){
plot_tree(cex=0.5,paste0(prefix,edge))
title(paste(edge,"migration(s)"))
}
dev.off()

#residuals
pdf(file = "TreemixPanel_Residuals_All.pdf")
par(mfrow=c(4,3))
for(edge in 0:10){
plot_resid(stem=paste0(prefix,edge),pop_order="pops")
}
dev.off()

pdf(file = "TreemixPanel_Residuals_m0_1.pdf")
par(mfrow=c(2,1))
for(edge in 0:1){
plot_resid(stem=paste0(prefix,edge),pop_order="pops")
}
dev.off()

pdf(file = "TreemixPanel_Residuals_m2_5.pdf")
par(mfrow=c(3,2))
for(edge in 2:5){
plot_resid(stem=paste0(prefix,edge),pop_order="pops")
}
dev.off()

pdf(file = "TreemixPanel_Residuals_m6_10.pdf")
par(mfrow=c(3,2))
for(edge in 6:10){
plot_resid(stem=paste0(prefix,edge),pop_order="pops")
}
dev.off()

#variation explained:
sink(file = "VariationExplained.txt", append = FALSE)
for(edge in 0:10){
var <- get_f(paste0(prefix,edge))
print(paste0("edge=",edge,"   VariationExplained=",var))
}

# BOOTSTRAP
for N in $(seq 1 100) ;
do echo ${N};
sbatch Treemix_Sample_Boot.sh ${N};
done ;

for m in 0 1 2 ;
do echo $m
for N in $(seq 1 100) ;
do echo ${N};
sbatch Treemix_Pop_Boot.sh ${N} ${m};
done ;
done ;

```

Treemix_Sample_Boot.sh

```bash

mkdir temp.sample.${1}

cp PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.edit.gz temp.sample.${1}/

cd temp.sample.${1}

module load TreeMix/1.13-intel-2018b

#making a random seed
R1=$(echo $RANDOM % 100 | bc)
R2=$(echo $RANDOM % 100 | bc)
S=$(bc <<<"scale=0; ${R1} * ${R2} * ${1}")

## Run treemix
treemix -i PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.edit.gz -global -noss -k 500 -root RAZ_Anc -bootstrap -o PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.boot.${1} -seed ${S}

mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.boot.${1}* ../sample_like_boot/

cd ../
rm -r temp.sample.${1}
```

Treemix_Pop_Boot.sh

```bash


mkdir temp.pop.${1}.${2}

cp PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.edit.gz temp.pop.${1}.${2}/

cd temp.pop.${1}.${2}/

## Load modules
module load TreeMix/1.13-intel-2018b

#making a random seed
R1=$(echo $RANDOM % 100 | bc)
R2=$(echo $RANDOM % 100 | bc)
S=$(bc <<<"scale=0; ${R1} * ${R2} * ${1}")

## Run treemix
treemix -i PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.edit.gz -m ${2} -global -k 500 -root RAZ_Anc -bootstrap -o PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.boot.${1}.m${2} -seed ${S}

mv PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.boot.${1}.m${2}* ../pop_like_boot/

cd ../
rm -r temp.pop.${1}.${2}

```

THEN

```bash

#pops

for i in 0 1 2 ;
do touch PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m${i}.boot_tree ;
for N in {1..100} ;
do zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.boot.${N}.m${i}.treeout.gz | head -1  >> PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m${i}.boot_tree ;
done ;
rm *boot*m${i}.*.gz ;
rm *boot*m${i}.llik ;
module purge ;
module load IQ-TREE/1.6.12-foss-2018b ;
iqtree -o RAZ_Anc -sup PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m${i}.main_tree PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.m${i}.boot_tree ;
done ;

#check mirgrations
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.*.m1.treeout.gz
zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.*.m2.treeout.gz

#samples
touch PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.boot_tree ;
for N in {1..100} ;
do zcat PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.boot.${N}.treeout.gz | head -1  >> PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.boot_tree ;
done ;
rm *boot*.gz
rm *boot*.llik
module purge ;
module load IQ-TREE/1.6.12-foss-2018b ;
iqtree -o RAZ_Anc -sup PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.main_tree PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.samples.m0.boot_tree ;


THEN

TreeInR.pop.R
```

#Migration

```bash
#MIGRATION
Vert=$(ls pop_like_boot/ | grep -v "boot" | grep "vertices")
Edge=$(ls pop_like_boot/ | grep -v "boot" | grep "edges")
for N in $(seq 1 10) ;
do echo ${N};
for M in {1..100} ;
do echo ${M} ;
sbatch Treemix_Pop_Mig.sh ${N} ${Vert} ${Edge} ${M};
done ;
done ;

```

Treemix_Pop_Mig.sh

```bash

mkdir temp.pop.${1}.${4}

cp pop_like_boot/${2} temp.pop.${1}.${4}/
cp pop_like_boot/${3} temp.pop.${1}.${4}/
cp ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.edit.gz temp.pop.${1}.${4}/

cd temp.pop.${1}.${4}/

sleep 120 

module load TreeMix/1.13-intel-2018b

#making a random seed
R1=$(echo $RANDOM % 100 | bc)
R2=$(echo $RANDOM % 100 | bc)
S=$(bc <<<"scale=0; ${R1} * ${R2} * ${1}")

## Run treemix
treemix -i ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.edit.gz -m ${1} -se -root RAZ -g ${2} ${3} -o ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m${1}.${4} -seed ${S}

mv ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m${1}.${4}* ../pop_like_boot/

cd ../
rm -r temp.pop.${1}.${4}
```

Plotting

```bash
zcat ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m0.1.covse.gz | head -1 | tr ' ' '\n' | head -13 > pops

#get highest likelihood m

#pops
for i in cov covse edges modelcov treeout vertices ; do mv ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m0.1.${i}.gz ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m0.${i}.gz ; done
mv ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m0.1.llik ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m0.llik

#put plotting_funcs.R in same folder 
touch VariationExplained.txt
# https://speciationgenomics.github.io/Treemix/

#use R with script plotting_funcs.R
R
prefix="ModernPuffinAngsd_NoZ_NoUnpl.IBS.ldpruned.withOut.treemix.pops.m"
library(RColorBrewer)
library(R.utils)
source("plotting_funcs.R")
#trees
par(mfrow=c(2,2))
for(edge in 0:3){
plot_tree(cex=0.5,paste0(prefix,edge))
title(paste(edge,"edges"))
}
#residuals
par(mfrow=c(1,1))
for(edge in 0){
plot_resid(stem=paste0(prefix,edge),pop_order="pops")
}
#variation explained:
sink(file = "VariationExplained.txt", append = FALSE)
for(edge in 0:12){
var <- get_f(paste0(prefix,edge))
print(paste0("edge=",edge,"   VariationExplained=",var))
}
sink(file=NULL)
# ctrl + d to escape

#then
sed 's/VariationExplained=/,/g' VariationExplained.txt | sed 's/"/,/g' | sed 's/=/,/g' > VariationExplained.Edit.txt
Plot.Variation.R

```

F3 stats

```bash
sbatch Treemix_Pop_F3.sh
sbatch Treemix_Pop_F3_edit.sh

```

Treemix_Pop_F3.sh

```bash


mkdir temp.f3

cp PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters.gz temp.f3/

cd temp.f3/

module load TreeMix/1.13-intel-2018b

## Run treemix
threepop -i PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters.gz -k 500 > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters.F3stat

rm PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters.gz
cp * ../f3/
cd ../
rm -r temp.f3/

```

Treemix_Pop_F3_edit.sh

```bash


mkdir temp.f3_edit

cp PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz temp.f3_edit/

cd temp.f3_edit/


module load TreeMix/1.13-intel-2018b

## Run treemix
threepop -i PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz -k 500 > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.F3stat

rm PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.pops.gz
cp * ../f3/
cd ../
rm -r temp.f3_edit/
```

```bash

#remove unneccessary lines and add header
echo "Admixed,A,B,f3,stdErr,Zscore" > header
tail -n+5 PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters.F3stat | grep -v "total_nsnp" | grep -v "Estimating f_3" | sed 's/;/,/g' | sed 's/ /,/g' | cat header - > PuffinAngsd_IBS_withOutgroup_NoZ_NoUn.ldpruned.treemix.Clusters.F3stat.edit
rm header


```
