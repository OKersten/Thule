ngsLD

```bash
cp /XX/XX/ChromosomeList .
sed -i 's/Scaffolds_//g' ChromosomeList
sed -i 's/://g' ChromosomeList
sed -i 's/unplaced/unplaced_/g' ChromosomeList

cat ChromosomeList | while read i ;
do echo ${i} ;
sbatch Puffin_ngsLD_decay.sh ${i} ;
done

```

Puffin_ngsLD_decay.sh

```bash

mkdir temp_${1}
cd temp_${1}

cp -r /XX/XX/PuffinAngsd_${1}.mafs.gz .
cp -r /XX/XX/PuffinAngsd_${1}.beagle.gz .

zcat PuffinAngsd_${1}.mafs.gz | cut -f 1,2 | tail -n +2 > PuffinAngsd_${1}.pos
NS=$(cat PuffinAngsd_${1}.pos | wc -l)
echo "Number of sites = "$NS

#ngsLD
module load ngsLD/191108-GCC-8.2.0-2.31.1

ngsLD --n_threads 16 --geno PuffinAngsd_${1}.beagle.gz --pos PuffinAngsd_${1}.pos --probs --n_ind 77 --n_sites $NS --max_kb_dist 50 | gzip --best > PuffinAngsd_${1}.ld.gz

#DecayPlot
#subsample output (5% of lines) for Decay Plot
zcat PuffinAngsd_${1}.ld.gz | awk 'rand()<0.05' | gzip --best > PuffinAngsd_${1}.ld_sampled.gz

cp -r PuffinAngsd_${1}.pos ../LD/
cp -r PuffinAngsd_${1}.ld.gz ../LD/
cp -r PuffinAngsd_${1}.ld_sampled.gz ../LD/

cd ../

rm -r temp_${1}/
```

Pruning

```bash
cat ChromosomeList | while read i ;
do echo ${i} ;
sbatch Puffin_LDpruning.sh ${i} ;
done

```

Puffin_LDpruning.sh

```bash

mkdir temp_${1}
cd temp_${1}

cp -r /XX/XX/PuffinAngsd_${1}.beagle.gz .
cp -r /XX/XX/PuffinAngsd_${1}.ld.gz .

module load OrthoMCL/2.0.9-intel-2018b-Perl-5.28.0

gunzip PuffinAngsd_${1}.ld.gz
gunzip PuffinAngsd_${1}.beagle.gz
tail -n +2 PuffinAngsd_${1}.beagle > PuffinAngsd_${1}.beagle.edit

mcl <(cut -f1,2,7 PuffinAngsd_${1}.ld | awk '$3>0.2') --abc -I 2.0 -o PuffinAngsd_${1}.ld.mcl

#non-cluster snps
awk '{print $1"\t"$2"\t"$7}' PuffinAngsd_${1}.ld | awk '$3>0.2' | awk '{print $1}' > tmp1
awk '{print $1"\t"$2"\t"$7}' PuffinAngsd_${1}.ld | awk '$3>0.2' | awk '{print $2}' > tmp2
cat tmp1 tmp2 | sort | uniq > cluster.snps
sed -i 's/:/_/g' cluster.snps
rm tmp1
rm tmp2
LANG=C fgrep -v -w -f cluster.snps PuffinAngsd_${1}.beagle.edit > tmp1

#rep snps of cluster
touch Rep.snps
TotalLines=$(cat PuffinAngsd_${1}.ld.mcl | wc -l)
echo "Number of total lines = "$TotalLines
for line in $(seq 1 1 $TotalLines) ;
do ColumnNumbers=$(head -${line} PuffinAngsd_${1}.ld.mcl | tail -1 | awk '{print NF}' | sort -nu | tail -n 1) ;
AlmostMiddleColumn=$(bc <<<"scale=0; $ColumnNumbers / 2" ) ;
let "MiddleColumn=$AlmostMiddleColumn+1" ;
Scaffold=$(head -${line} PuffinAngsd_${1}.ld.mcl | tail -1 | awk -v VAR=$MiddleColumn '{print $VAR}') ;
echo $Scaffold >> Rep.snps ;
done ;

cat Rep.snps | sort | uniq > Rep.snps.uniq
sed -i 's/:/_/g' Rep.snps.uniq

LANG=C fgrep -w -f Rep.snps.uniq PuffinAngsd_${1}.beagle.edit > tmp2

#merging
cat tmp1 tmp2 | sort -t '_' -n -k4,4 > PuffinAngsd_${1}.beagle.tmp
head -1 PuffinAngsd_${1}.beagle | cat - PuffinAngsd_${1}.beagle.tmp > PuffinAngsd_${1}.beagle.ldpruned
gzip PuffinAngsd_${1}.beagle.ldpruned

cp -r PuffinAngsd_${1}.beagle.ldpruned.gz ../LD/
cd ../

rm -r temp_${1}/

```

Combine chromosomes

```bash

#ALL
zcat PuffinAngsd_chromosome_1.beagle.ldpruned.gz > PuffinAngsd_All.beagle.ldpruned

for i in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 2 3 4 5 6 7 8 9 Z ;
do zcat PuffinAngsd_chromosome_${i}.beagle.ldpruned.gz | tail -n +2 >> PuffinAngsd_All.beagle.ldpruned ;
done ;
zcat PuffinAngsd_unplaced_.beagle.ldpruned.gz | tail -n +2 >> PuffinAngsd_All.beagle.ldpruned ;
gzip PuffinAngsd_All.beagle.ldpruned
zcat PuffinAngsd_All.beagle.ldpruned.gz | tail -n +2 | wc -l


#No Z or unplaced
zcat PuffinAngsd_chromosome_1.beagle.ldpruned.gz > PuffinAngsd_NotAll.beagle.ldpruned
for i in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 2 3 4 5 6 7 8 9;
do zcat PuffinAngsd_chromosome_${i}.beagle.ldpruned.gz | tail -n +2 >> PuffinAngsd_NotAll.beagle.ldpruned ;
done ;
gzip PuffinAngsd_NotAll.beagle.ldpruned
mv PuffinAngsd_NotAll.beagle.ldpruned.gz PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz
zcat PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz | tail -n +2 | wc -l

```