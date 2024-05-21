# 2024_Calypsoinae
2024_Calypsoinae_phylogenomics

# Commands for running HybPiper and building ASTRAL trees

## The Angiosperm353 kit was used, sequenced on a NextSeq2000 p1 flow cell (~100M 2x150 bp reads)

#################

```bash
cd /data/cbarrett/2024_Calypsoinae


## make a reads dir, and edit file names

sudo mkdir reads
sudo mv *.fastq.gz reads
cd reads
rename 's/_001.fastq.gz/.fastq.gz/g' *.fastq.gz
rename 's/_L001//g' *.fastq.gz

# make a clean reads dir and run FASTP (trim adapters, poly-X, low quality bases)

sudo mkdir ../fastp_polyg

for f1 in *_R1.fastq.gz
	do
        f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
        fastp -i $f1 -I $f2 -w 16 --trim_poly_g --trim_poly_x -l 75 --cut_right -o "../fastp_polyg/fastp-$f1" -O "../fastp_polyg/fastp-$f2"
	done


# make a namefile and edit with sed

# output from ls command, saved to new text file
ls *R1.fastq.gz > ../namelist.txt

# sed replace "_R1.fastq" in namefile
sed 's/_R1.fastq.gz//g' namelist.txt > namelist2.txt

# check if everything is okay
cat namelist2.txt

# If needed, delete sample names from namelist if they have too few reads or failed in sequencing
```

#################

## Run Hybpiper

```bash
## make an output directory and cd to it

mkdir hybpiper_2024_01_08
sudo chmod 777 -R hybpiper_2024_01_08
cd hybpiper_2024_01_08

## run the pipeline -- all the file locations are specified via relative paths below in the while loop

conda activate hybpiper

# [OLD Command, for hybpiper v 1]
while read name 
do python /usr/local/src/HybPiper/reads_first.py \
  -b ../Angiosperms353_targetSequences.fasta \
  -r ../fastp_polyg/"$name"_R*.fastq.gz \
  --prefix $name \
  --bwa
done < ../namelist2.txt

# Proper command for hybpiper2

while read name; 
 
do hybpiper assemble -t_dna ../Angiosperms353_targetSequences.fasta -r ../fastp_polyg/$name*.fastq.gz --prefix $name --bwa; 
 
done < ../namelist2.txt

```

## alignments with mafft (auto)

```bash

for i in *.FNA; do
mafft --adjustdirection --thread 32 ${i} > ${i%.*}_mafft.fasta;
done

## May need to remove '_R_' from beginning of fasta headers, inserted by mafft, to indicate reverse complements

sed -i 's/_R_//g' *_mafft.fasta

```

## Trim alignments, but first convert 'n' to '-'

```bash

parallel seqkit -is replace -p "n" -r "-" {} -o {.}_seqkit.fasta ::: *.fasta 

parallel trimal -in {} -out {.}_trim.fasta -htmlout {.}_trim.html -automated1 ::: *_seqkit.fasta

mkdir ../trimmed_aligned_paralogs
sudo mv *_seqkit_trim.fasta ../trimmed_aligned_paralogs
cd ../trimmed_aligned_paralogs

```

## Build trees

```bash

iqtree2 -S trimmed_aligned_paralogs --prefix iqtree_paralogs -m GTR+G+I -T 32

```


## Remove '.main' '.0' '.1' annotation in paralogs for ASTER/A-Pro

```bash

sed 's/.main//g' < iqtree_paralogs.treefile | sed 's/\.\d:/:/g' > iqtree_paralogs2.treefile

```

## Used the ASTER GUI for windows, but can run A-pro from commandline (Astral-Pro, standard search, local PP)


## Reroot all trees for PhyParts

```bash

nano outgroup.list  ## fastp-Brassavola-glauca_S40

python /data/cbarrett/hybpiper_test/reroot_trees.py iqtree4.treefile outgroup.list > rerooted.treefile  

```

## Run phyparts for quartets

```bash

java -jar /usr/local/bin/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rerooted.treefile -m inferred_species_tree4_rooted.tre -o output_file_phyparts

## Plot pie charts on tree

python phypartspiecharts.py inferred_species_tree4_rooted.tre output_file_phyparts 288 --svg_name piechart_tree.svg

# conduct gene and site concordance factor analyses
iqtree2 -t inferred_species_tree4_rooted.tre -s calypsoinae_all.nex --gcf rerooted2_fixed.tre --scf 10000 --prefix concord_astraltree -T 28

```

# R commands to plot and compare trees

## Comparing tree with one "best" paralog selected vs. including all paralogs

```{r}

# Read trees
tre<-read.tree("Calypsoinae_astral_wo_and_w_paralogs.tre")

# Root the second tree (was not rooted)
tre2 <-root(tre,outgroup = "fastp-Brassavola-glauca_S40",resolve.root = TRUE)
is.rooted(tre2)
TRUE TRUE

# Drop the weird Dactylostalix species
tre3 <- drop.tip(tre2,"fastp-Dactylostalix-ringensoides_S41")

# separate out trees from multiphylo to phylo objects, ladderize the second one which wasn't ladderized
noparalogs <- tre3[[1]]
paralogs <- tre3[[2]]
para2<-ladderize(paralogs, right = FALSE)

# Calculate Robinson-Foulds and SPR distance
RF.dist(noparalogs,para2)
[1] 2

SPR.dist(noparalogs,para2)

spr 
  1 # Basically, one SPR rearrangement is needed to get from one tree to the next (placement of Cor bulbosa)

```









