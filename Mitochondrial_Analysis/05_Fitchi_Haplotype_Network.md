Run Fitchi with gene tree topology, but AllSNP/NoNs haplotypes

```bash
module load Python/3.7.2-GCCcore-8.2.0
python ElConcatenero.py -m N -if fasta -of nexus -in All.CompleteMitoRef.SNPs.withOut.NoNs.edit.fasta -o All.CompleteMitoRef.SNPs.withOut.NoNs.edit

sed -i 's/mixed ()/DNA/g' All.CompleteMitoRef.SNPs.withOut.NoNs.edit.nex
sed -i 's/interleave=no//g' All.CompleteMitoRef.SNPs.withOut.NoNs.edit.nex

#Run RemoveBootstrap.R
#library(ape) # Load the ape package
#tree1 <- read.tree("Puffin.mt.part.nex.edit.treefile") # A normal Newick with bootstraps from RAxML, for example
#tree1$node.label <- NULL # Erase the bootstrap values from the phylo object
#write.tree(tree1, file = "Puffin.mt.part.nex.edit.nobootstrap.treefile") # Save it


cp All.CompleteMitoRef.SNPs.withOut.NoNs.edit.nex All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.nex

# add tree at end of nexus from treefile in textwrangler
begin trees;

tree ml = [&R] >TREE<

end;

#change sample names to include pops
cp All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.nex All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex

for i in BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU ;
do sed -i '' "s/${i}/${i}_/g" All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex ;
done ;

#change the outlier names to something short
sed -i '' "s/Aethia_cristatella/Aet_cri/g" All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex
sed -i '' "s/Alca_torda/Alc_tor/g" All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex
sed -i '' "s/Synthliboramphus_antiquus/Syn_ant/g" All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex
sed -i '' "s/Synthliboramphus_wumizusume/Syn_wum/g" All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex

conda activate fitchi
python3 fitchi.py All.CompleteMitoRef.SNPs.withOut.NoNs.haplotypes.edit.nex All.SNPs.pops.withOut.fitchi.html --haploid -m 0.002 -p BJO BRE FAR GRI HOR IOM PAP ROS SPI WES GAN GUL THU

python3 /Users/oliverke/BioinformaticsProgs/Fitchi/fitchi_extract.py -e svg All.SNPs.pops.withOut.fitchi.html > All.SNPs.pops.withOut.fitchi.svg

conda deactivate

```