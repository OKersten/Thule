
```bash

cp /XX/XX/PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz  .

# NgsAdmix for K from 1 to 10
for K in $(seq 1 10) ;
do echo ${K};
for i in $(seq 1 50) ;
do sbatch NgsAdmix.sh ${K} ${i};
done ;
done ;


```

NgsAdmix.sh

```bash

mkdir temp.${1}.${2}
cd temp.${1}.${2}
cp ../PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz .

#making a random seed
R=$(echo $RANDOM % 100 | bc)
S=$(bc <<<"scale=0; ${R} * ${1} * ${2}")

#ngsAdmix
module load NGSadmix/32-GCC-7.3.0-2.30

Threads=16

NGSadmix -likes PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz -K ${1} -P ${Threads} -o PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned_ngsAdmix_k${1}_n${2} -maxiter 5000 -seed ${S}

rm -r PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz
cp -r * ../
cd ../
rm -r temp.${1}.${2}
```

Then

```bash

#log likelihoods for best K
(for log in `ls *.log`; do grep -Po 'like=\K[^ ]+' $log; done) > logfile

#Zip file for Clumpak run
for i in 1 2 3 4 5 6 7 8 9 10 ;
do zip K${i}_PuffinAngsd_NoZ_NoUnpl.zip PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned_ngsAdmix_k${i}_*.qopt ;
done ;
zip PuffinAngsd_NoZ_NoUnpl.zip K*_PuffinAngsd_NoZ_NoUnpl.zip


#use NgsAdmix2Clumpak.R

#sort file
sort -n -k1,1 logfile_formatted.txt | awk -F ' ' '{print $1"\t"$2}' > logfile_formatted_clumpak.txt

#then clumpak ; http://clumpak.tau.ac.il/bestK.html
# with logfile and zip file
# save output zipfolders as EvannoOut and ClumpakOut

#check best K's (Evanno) and select mean Qmatrix file based Clumpakk

awk '{print $6"\t"$7}' Clumpak/1621093048/K=2/MajorCluster/CLUMPP.files/ClumppIndFile.output > K2_Clumpak_MajorQmatrix.txt
awk '{print $6"\t"$7"\t"$8}' Clumpak/1621093048/K=3/MajorCluster/CLUMPP.files/ClumppIndFile.output > K3_Clumpak_MajorQmatrix.txt
awk '{print $6"\t"$7"\t"$8"\t"$9}' Clumpak/1621093048/K=4/MajorCluster/CLUMPP.files/ClumppIndFile.output > K4_Clumpak_MajorQmatrix.txt

# then Admixture_GenoLike.R for Plot

```

