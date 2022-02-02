
```bash
cp /XX/XX/PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz .

module load PCAngsd/200115-foss-2019a-Python-2.7.15

pcangsd.py -beagle PuffinAngsd_NoZ_NoUnpl.beagle.ldpruned.gz -snp_weights -selection -sites_save -o Puffin.NoZ_NoUnpl.LDpruned.pcangsd -threads 16 > PCAngsd.log


```

Then plot with R/ggplot2

