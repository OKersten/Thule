Inbreeding according to Kersten et al. 2021

```bash

cat bams_good | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
cp /XX/XX/IndHeterozygosity/Results/${sample}.* . ;
done

rm *.sfs

module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="16"

sample=BJO001

realSFS ${sample}.saf.idx ${sample}.saf.idx -P ${Threads} -fold 1 > test.ml
realSFS fst index ${sample}.saf.idx ${sample}.saf.idx -sfs test.ml -whichFst 1 -fstout test -fold 1
realSFS fst stats2 test.fst.idx -win 100000 -step 50000 -fold 1 > slidingwindow #this will determine RoH size
rm test.fst.gz
touch Regions.txt
cat slidingwindow | tail -n+2 | awk '{print $2}' | uniq | while read chromo ;
do echo ${chromo} ;
echo -e ${chromo}":1-100000" >> Regions.txt ;
grep -w ${chromo} slidingwindow | while read line ;
do Region=$(echo ${line} | awk '{print $1}' | awk -F '(' '{print $4}' | sed 's/)//g' | sed 's/,/-/g') ;
echo -e ${chromo}":"${Region} >> Regions.txt ;
done ;
done ;
######
wc -l Regions.txt

#split into 22974 (42*547 )+ 33

cat bams_good | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
sbatch realSFS_ROH.sh ${sample} ;
done ;


```

realSFS_ROH.sh

```bash


module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="1"

touch ${1}_${SLURM_ARRAY_TASK_ID}.sfs
cp Regions.txt Regions_${1}_${SLURM_ARRAY_TASK_ID}.txt

NumberOfLines=$(echo "547") 
BottomLine=$(bc <<<"scale=0; ${SLURM_ARRAY_TASK_ID} * 547 ") 

cat Regions_${1}_${SLURM_ARRAY_TASK_ID}.txt | head -${BottomLine} | tail -${NumberOfLines} | while read line ;
do realSFS ${1}.saf.idx -r ${line} -P ${Threads} -fold 1 >> ${1}_${SLURM_ARRAY_TASK_ID}.sfs ;
done ;

rm Regions_${1}_${SLURM_ARRAY_TASK_ID}.txt

```

realSFS_ROH_end.sh

```bash


module load angsd/0.931-GCC-8.2.0-2.31.1

Threads="1"

touch ${1}_43.sfs
cp Regions.txt Regions_${1}_43.txt

cat Regions_${1}_43.txt | tail -33 | while read line ;
do realSFS ${1}.saf.idx -r ${line} -P ${Threads} -fold 1 >> ${1}_43.sfs ;
done ;

rm Regions_${1}_43.txt

```

THEN

```bash

cat bams_good | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
touch ${sample}.sfs ;
for i in {1..43} ;
do cat ${sample}_${i}.sfs >> ${sample}.sfs ;
rm ${sample}_${i}.sfs ;
done ;
done ;

nano HeteroRoh.R
### HeteroRoh.R ###
infile <- commandArgs(trailingOnly = TRUE)
sample <- substr(infile, 1, 6)
df <- as.data.frame(read.table(file = infile, stringsAsFactors = FALSE))
df$V4 <- rowSums(df)
df$V5 <- df$V2 / df$V4
df$V6 <- sample
sample_table <- as.data.frame(cbind(df$V6, df$V5))
write.table(sample_table, file = paste(sample,"heterozygosity.txt", sep = "_"), col.names = F, row.names = F, quote = F, sep = "\t")
###############

module purge
module load R/3.6.0-intel-2019a
cat bams_good | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
paste Regions.txt ${sample}.sfs > ${sample}.regionSFS ;
Rscript HeteroRoh.R ${sample}.sfs ;
paste ${sample}.regionSFS ${sample}_heterozygosity.txt > ${sample}.RoH
done ;

cat *heterozygosity.txt > HeterozygosityCombined.txt

# plot data to assess cutoff


#then
nano ROH_Length_Count.R #change to RoH sizes

##### ROH_Length_Count.R
library(tidyr)
library(tidyverse)
infile <- commandArgs(trailingOnly = TRUE)
sample <- substr(infile, 1, 6)
df <- as.data.frame(read.table(file = infile, stringsAsFactors = FALSE))
df_edit <- cbind(df$V5, data.frame(do.call('rbind', strsplit(as.character(df$V1),':',fixed=TRUE))), df$V6)
colnames(df_edit) <- c("Sample", "Chromosome", "Region", "Hetero")
for (i in 1:25) {
df_edit2 <- df_edit %>% filter(Chromosome == paste("Scaffolds_chromosome",i, sep = "_"))
df_edit2$ROH <- ifelse(df_edit2$Hetero < 0.001441109, 1, 0) # enter cutoff here
value <- rle(df_edit2$ROH)$values
lengths <- as.numeric(rle(df_edit2$ROH)$lengths)
df2 <- cbind.data.frame(value, lengths)
df2$chromosome <- paste("Chromosome",i, sep = "_")
df3 <- df2 %>% filter(value == 1) %>% filter(lengths > 1)
df3$ROH_Length <- df3$lengths * 0.1 - ((df3$lengths - 1) * 0.05 )
assign(paste("Chromo",i,"df3", sep = "_"), df3)
}
Dataframes <- lapply(paste("Chromo_",1:25,"_df3", sep = ""), get)
All_chromoMerge_df3 <- bind_rows(Dataframes)
All_chromoMerge_df3$value <- sample
colnames(All_chromoMerge_df3) <- c("Sample", "Number", "Chromosome","ROH_Length")
write.table(All_chromoMerge_df3, file = paste(sample,"_ROH_Lengths.txt", sep = ""), row.names = F, quote = F, sep = "\t")
sample_inbreeding <- unname(colSums(All_chromoMerge_df3[,"ROH_Length", drop = FALSE]) / (nrow(df_edit) * 0.1 - ((nrow(df_edit) - 1) * 0.05)))
Inbreeding_table <- as.data.frame(cbind(sample, sample_inbreeding))
write.table(Inbreeding_table, file = paste(sample,"_F_ROH.txt", sep = ""), row.names = F, quote = F, sep = "\t")
###############

module purge
module load R/3.6.0-intel-2019a

echo -e "Sample\tNumber\tChromosome\tROH_Length" > AllSamples_ROH_Length
echo -e "Sample\tF_ROH" > AllSamples_F_ROH
cat bams_good | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
Rscript ROH_Length_Count.R ${sample}.RoH ;
cat ${sample}_ROH_Lengths.txt | tail -n+2 >> AllSamples_ROH_Length ;
cat ${sample}_F_ROH.txt | tail -n+2 >> AllSamples_F_ROH ;
done ;


#########
check ROHs in detail
#########
nano RoH_Distribution.R
##### RoH_Distribution.R
library(tidyr)
library(tidyverse)
infile <- commandArgs(trailingOnly = TRUE)
sample <- substr(infile, 1, 6)
df <- as.data.frame(read.table(file = infile, stringsAsFactors = FALSE))
df_edit <- cbind(df$V5, data.frame(do.call('rbind', strsplit(as.character(df$V1),':',fixed=TRUE))), df$V6)
colnames(df_edit) <- c("Sample", "Chromosome", "Region", "Hetero")

for (i in 1:25) {
  df_edit2 <- df_edit %>% filter(Chromosome == paste("Scaffolds_chromosome",i, sep = "_"))
  df_edit2$ROH <- ifelse(df_edit2$Hetero < 0.001441109, 1, 0) # enter cutoff here
  df_edit3 <- df_edit2 %>% filter(ROH == "1") %>% separate("Region", c("Start", "End"), sep = "-") %>% 
    arrange("Start") %>% 
    group_by(g = cumsum(cummax(lag(End, default = first(Start))) < Start)) %>% 
    summarise(Start = first(Start), End = max(End))
  df_edit3$Chromosome <- c(paste("Chromosome",i, sep = "_"))
  df_edit3$Start <- as.numeric(df_edit3$Start)
  df_edit3$End <- as.numeric(df_edit3$End)
  df_edit3$RoH_Length <- (df_edit3$End - df_edit3$Start) / 1000000
  df_edit3$Sample <- sample
  df_edit4 <- subset(df_edit3, select = -c(g))
  df_edit5 <- df_edit4[, c("Sample", "Chromosome", "Start", "End", "RoH_Length")]
  assign(paste("Chromo",i,"df2", sep = "_"), df_edit5)
}

Dataframes <- lapply(paste("Chromo_",01:25,"_df2", sep = ""), get)
All_chromoMerge_df3 <- bind_rows(Dataframes)
write.table(All_chromoMerge_df3, file = paste(sample,"_ROH_Distribution.txt", sep = ""), row.names = F, quote = F, sep = "\t")
#########

echo -e "Sample\tChromosome\tStart\tEnd\tRoH_Length" > AllSamples_ROH_Distribution
module purge
module load R/4.0.3-foss-2020b
cat bams_good | awk -F '.' '{print $1}' | while read sample ;
do echo ${sample} ;
Rscript RoH_Distribution.R ${sample}.RoH ;
cat ${sample}_ROH_Distribution.txt | tail -n+2 >> AllSamples_ROH_Distribution ;
done ;

```
