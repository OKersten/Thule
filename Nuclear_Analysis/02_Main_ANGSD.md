
```bash

echo -e "Scaffolds_chromosome_10:\nScaffolds_chromosome_11:\nScaffolds_chromosome_12:\nScaffolds_chromosome_13:\nScaffolds_chromosome_14:\nScaffolds_chromosome_15:\nScaffolds_chromosome_16:\nScaffolds_chromosome_17:\nScaffolds_chromosome_18:\nScaffolds_chromosome_19:\nScaffolds_chromosome_1:\nScaffolds_chromosome_20:\nScaffolds_chromosome_21:\nScaffolds_chromosome_22:\nScaffolds_chromosome_23:\nScaffolds_chromosome_24:\nScaffolds_chromosome_25:\nScaffolds_chromosome_2:\nScaffolds_chromosome_3:\nScaffolds_chromosome_4:\nScaffolds_chromosome_5:\nScaffolds_chromosome_6:\nScaffolds_chromosome_7:\nScaffolds_chromosome_8:\nScaffolds_chromosome_9:\nScaffolds_chromosome_Z:\nScaffolds_unplaced:" > ChromosomeList

#run Puffin_ANGSD_Abel.sh
cat ChromosomeList | while read i ;
do echo ${i} ;
sbatch Puffin_ANGSD_main.sh ${i} ;
done
```

Puffin_ANGSD_main.sh

```bash

mkdir temp.${1}
cp -r Assembly_Puffin_NU.MT.fasta temp.${1}/ #reference genome
cp bams_good temp.${1}/ #list of bamfiles
cp -r *.Assembly_Puffin_NU.sorted.bam temp.${1}/
cp -r *.Assembly_Puffin_NU.sorted.bam.bai temp.${1}/

cd temp.${1}

Chromosome=$(echo ${1} | awk -F ':' '{print $1}' | awk -F '_' '{print $2"_"$3}')
#script is called per chromosome and file names will be named with this EnvVariable


FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -dosnpstat 1 -C 50 -baq 2 -checkBamHeaders 1 -doHWE 1 -HWE_pval 1e-2 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -snp_pval 1e-6 -minMaf 0.05 -setMaxDepth 635 -setMinDepth 365"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 8 -doPost 1 -doGlf 2"

module load SAMtools/1.9-GCC-8.2.0-2.31.1
samtools faidx Assembly_Puffin_NU.MT.fasta

module purge
module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="16"

angsd -b bams_good -GL 1 $FILTERS $TODO -P ${Threads} -ref Assembly_Puffin_NU.MT.fasta -r ${1} -out PuffinAngsd_${Chromosome}

cp -r PuffinAngsd_${Chromosome}* /cluster/work/users/oliverke/Thule/ANGSD/GenoLikelihoods/

cd ..

rm -r temp.${1}

```

SUMMARY

```bash

zcat PuffinAngsd_chromosome_1.mafs.gz > PuffinAngsd_All.mafs

for i in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 2 3 4 5 6 7 8 9 Z ;
do zcat PuffinAngsd_chromosome_${i}.mafs.gz | tail -n +2 >> PuffinAngsd_All.mafs ;
done ;
zcat PuffinAngsd_unplaced_.mafs.gz | tail -n +2 >> PuffinAngsd_All.mafs ;
gzip PuffinAngsd_All.mafs
zcat PuffinAngsd_All.mafs.gz | tail -n +2 | wc -l

```