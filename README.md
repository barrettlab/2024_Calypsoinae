# 2024_Calypsoinae
2024_Calypsoinae_phylogenomics
## The Angiosperm353 kit was used, sequenced on a NextSeq2000 p1 flow cell (~100M 2x150 bp reads)

#################

## Code used in Testing the monophyly of Corallorhiza and Oreorchis (Epidendroideae, Epidendreae, Calypsoinae) with nuclear sequence capture 

## Craig F. Barrett1*, John V. Freudenstein2, Samuel V. Skibicki1, Brandon T. Sinn3,4, Tomohisa Yukawa5, Kenji Suesugu6,7

## Installing hybpiper with Conda

```bash
conda create -n hybpiper bioconda::hybpiper # This creates a new environment and install hybpiper into that
conda activate hybpiper  # This should work, but if it doesn't, try:

conda config --set channel_priority flexible
conda activate hybpiper
conda install bioconda::hybpiper
```

## Quality and adapter trimming with fastp

```bash
cd /data/cbarrett/2024_Calypsoinae


## make a reads dir, and edit file names

sudo mkdir reads
sudo mv *.fastq.gz reads
cd reads

# Ubuntu (myco server, Thorny Flat)
rename 's/_001.fastq.gz/.fastq.gz/g' *.fastq.gz
rename 's/_L001//g' *.fastq.gz

# RedHat (newer server)
rename "_001.fastq.gz" ".fastq.gz"  *.fastq.gz
rename "_L001" ""  *.fastq.gz

## make a clean reads dir and run FASTP (trim adapters, poly-X, low quality bases)

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
sed 's/_R1.fastq.gz//g' namelist.txt > namelist3.txt

# check if everything is okay
cat namelist3.txt

# If needed, delete sample names from namelist if they have too few reads or failed in sequencing
```

## Hybpiper assembly, stats, sequence generation, and paralog retriever

```bash
## make an output directory and cd to it

mkdir hybpiper_2024_01_08
sudo chmod 777 -R hybpiper_2024_01_08
cd hybpiper_2024_01_08

## run the pipeline -- all the file locations are specified via relative paths below in the while loop
conda activate hybpiper

while read name; 
do hybpiper assemble -t_dna ../Angiosperms353_targetSequences.fasta -r ../fastp_polyg/$name*.fastq.gz --prefix $name --bwa; 
done < ../namelist3.txt

hybpiper stats -t_dna ../Angiosperms353_targetSequences.fasta gene ../namelist3.txt

## Get heatmap pdf
hybpiper recovery_heatmap seq_lengths.tsv --heatmap_filetype pdf --heatmap_dpi 1000

## Retrieve sequences, with "best" paralogs
hybpiper retrieve_sequences -t_dna ../Angiosperms353_targetSequences.fasta dna --sample_names ../namelist3.txt

## Retrieve sequences with all paralogs. 2 outputs dirs, paralogs_all and paralogs_no_chimeras (we want the second one)
hybpiper paralog_retriever ../namelist3.txt -t_dna ../Angiosperms353_targetSequences.fasta dna

conda deactivate

```

## A little housekeeping

```bash
mkdir best_paralog
sudo mv *.FNA best_paralog
cd best_paralog

```


## Plotting original vs newtargets recovery

```{r}

hyb <- read.csv("hybpiper_stats.csv",header=T)
library(ggplot2)
ggplot(hyb, aes(x=Name)) + 
    geom_point(aes(y = GenesWithSeqs), color = "red", size=3) + 
    geom_point(aes(y = GenesWithSeqs_newtargets), color="blue", size=3) + coord_flip()

```

## Alignment with mafft

```bash
for i in *.FNA; do
mafft --adjustdirection --thread 32 ${i} > ${i%.*}_mafft.fasta;
done
```

## TrimAL

```bash
mkdir seqkit 
mkdir trimal

# replace lower case 'n's with '-'
parallel seqkit -is replace -p "n" -r "-" {} -o seqkit/{.}_seqkit.fasta ::: *.fasta 

cd seqkit 

parallel trimal -in {} -out ../trimal/{.}_trim.fasta -htmlout ../trimal/{.}_trim.html -automated1 ::: *_seqkit.fasta

cd ../trimal
mkdir html
sudo mv *.html html

```

## Adjustments to data containing paralogs

```bash

cd paralogs_no_chimeras
mkdir seqkit_paralogs 
mkdir trimal_paralogs

## alignments with mafft (auto)

for i in *.fasta; do
mafft --adjustdirection --thread 32 ${i} > ${i%.*}_mafft.fasta;
done

parallel seqkit -is replace -p "n" -r "-" {} -o seqkit_paralogs/{.}_seqkit.fasta ::: *_mafft.fasta 

cd seqkit_paralogs

parallel trimal -in {} -out ../trimal_paralogs/{.}_trim.fasta -htmlout ../trimal_paralogs/{.}_trim.html -automated1 ::: *_seqkit.fasta

cd ../trimal_paralogs
mkdir html
sudo mv *.html html

```


## Running hybpiper on RNA-seq and DNA-seq-cap data simultaneously to identify potential A353 gene losses

```bash
## Created a text file called new_corallorhiza_namelist.txt, which contained the following lines:
## The first four are RNA-seq datasets, and the rest are A353 seq-cap (DNA-seq)

C_trifida
C_wisteriana
C_maculata
C_striata
Corallorhiza-trifida_S2
Corallorhiza-wisteriana-Eastern-US_S27
Corallorhiza-wisteriana-Western-US_S28
Corallorhiza-bentleyi_S3
Corallorhiza-bulbosa_S19
Corallorhiza-involuta_S1
Corallorhiza-macrantha_S21
Corallorhiza-mertensiana_S30
Corallorhiza-odontorhiza-Mexico_S47
Corallorhiza-odontorhiza-var-odontorhiza2_S38
Corallorhiza-odontorhiza-var-pringlei_S25
Corallorhiza-maculata-var-maculata2_S36
Corallorhiza-maculata-var-mexicana_S17
Corallorhiza-maculata-var-occidentalis_S15
Corallorhiza-striata-CA-Coast-Ranges_S11
Corallorhiza-striata-Sierra-Nevada_S9
Corallorhiza-striata-var-striata_S5
Corallorhiza-striata-var-vreelandii_S7


Hybpiper commands:

while read name;  do hybpiper assemble -t_dna mega353.fasta -r $name*.fastq --prefix $name --bwa; done < new_corallorhiza_namelist.txt

hybpiper stats -t_dna mega353.fasta gene new_corallorhiza_namelist.txt

hybpiper recovery_heatmap seq_lengths.tsv

## RNA-seq and GO analyses
ShinyGo: (http://bioinformatics.sdstate.edu/go/)

Panther GO overrepresentation testing: (https://pantherdb.org/tools/uploadFiles.jsp?#)

```


## IQtree best and all paralogs

```bash
## Do some additional filtering (need to automate this)
## Pulled all 'trimal' files into Geneious, sort by # sequences
## Remove all alignments with < 10 sequences, then sort by length
## remove all alignments with < 100 positions, export, pull into new dir called 'trimal_geneious'

cd best_paralog/trimal

## Need to remove '_R_' from fasta headers, stuck there from mafft *reversed sequences)

sed -i 's/_R_//g' *_trim.fasta

## remove paralog tags (.0, .1, .2, etc.) so all paralog leaves have identical names, as required by Astral-pro
sed 's/.main//g' < treeshrink_paralogs_output.tre | sed 's/\.\d:/:/g' > treeshrink_iqtree_paralogs2.treefile
cd ..

iqtree2 -S trimal_geneious --prefix iqtree_bestparalog -m GTR+G+I -T 32
```
## PhyKit to identify consistently long branches among gene trees

```bash

conda activate /usr/local/src/conda_envs/binf

# Input file containing gene trees in Newick format
gene_tree_file="output2.treefile"

# Output matrix file
output_matrix="lb_scores_matrix.txt"

# Temporary directory for storing intermediate tree files
temp_dir="./temp_trees"
mkdir -p "$temp_dir"

# Create an associative array to store the taxa and scores for each gene tree
declare -A lb_scores

# Initialize a list to track taxa
declare -a taxa_list

# Counter for gene trees
tree_counter=0

# Read each gene tree from the file
while read -r gene_tree; do
    # Skip empty lines
    if [ -z "$gene_tree" ]; then
        continue
    fi

    # Increment the tree counter
    tree_counter=$((tree_counter + 1))

    # Create a temporary file for the current gene tree
    temp_tree_file="${temp_dir}/temp_tree_${tree_counter}.treefile"
    echo "$gene_tree" > "$temp_tree_file"

    # Run phykit lb_score and capture the output
    lb_output=$(phykit lb_score "$temp_tree_file" --verbose)

    # Parse the output and store the scores in the associative array
    while read -r line; do
        if [[ "$line" =~ ^[^#]+ ]]; then
            taxon=$(echo "$line" | awk '{print $1}')
            score=$(echo "$line" | awk '{print $2}')
            
            # Store taxa in taxa_list if not already present
            if ! [[ " ${taxa_list[@]} " =~ " ${taxon} " ]]; then
                taxa_list+=("$taxon")
            fi
            
            # Store score in the associative array, key is taxon + gene tree number
            lb_scores["$taxon,$tree_counter"]="$score"
        fi
    done <<< "$lb_output"

done < "$gene_tree_file"

# Generate the matrix and write to the output file
{
    # Print the header row (Taxa + GeneTree columns)
    echo -n "Taxa"
    for i in $(seq 1 $tree_counter); do
        echo -n -e "\tGeneTree_$i"
    done
    echo

    # Print the matrix rows (Taxa + Scores for each gene tree)
    for taxon in "${taxa_list[@]}"; do
        echo -n "$taxon"
        for i in $(seq 1 $tree_counter); do
            score=${lb_scores["$taxon,$i"]}
            if [ -z "$score" ]; then
                echo -n -e "\tNA"  # Missing score
            else
                echo -n -e "\t$score"
            fi
        done
        echo
    done
} > "$output_matrix"

# Clean up temporary files
rm -rf "$temp_dir"

echo "Long branch score matrix saved to $output_matrix"

```

```{r}
# Load required libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
    install.packages("viridis")
}
library(ggplot2)
library(viridis)  # For the viridis color palette
library(dplyr)    # For data manipulation

# Read the lb_scores.csv file
lb_scores <- read.csv("lb_scores.csv", row.names = 1)

# Check the data
head(lb_scores)

# Calculate mean and standard error for each taxon
mean_scores <- rowMeans(lb_scores, na.rm = TRUE)
stderr_scores <- apply(lb_scores, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Create a data frame for plotting
plot_data <- data.frame(
  Taxon = rownames(lb_scores),
  Mean = mean_scores,
  SE = stderr_scores
)

# Remove "fastp-" from taxon names
plot_data$Taxon <- gsub("^fastp-", "", plot_data$Taxon)

# Extract genus (first part before the first "-")
plot_data$Genus <- sub("-.*$", "", plot_data$Taxon)

# Reorder the Taxon factor by mean LB score within each genus
plot_data <- plot_data %>%
  arrange(Genus, Mean) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))

# Create a color palette using viridis scaled by the mean LB score
plot_data$Color <- viridis(100)[cut(plot_data$Mean, breaks = 100)]
# used chatgpt to help with the plotting code
# Plot using ggplot2
ggplot(plot_data, aes(y = Taxon, x = Mean, fill = Mean)) +
  geom_col() +
  geom_errorbarh(aes(xmin = Mean - SE, xmax = Mean + SE), height = 0.2, color = "black") +
  scale_fill_viridis_c() +  # Apply viridis color palette scaled by mean LB scores
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, size = 8),  # Adjust y labels
    axis.title.y = element_blank()  # Remove y-axis title
  ) +
  labs(
    x = "Mean LB Score",
    title = "Mean LB Scores per Taxon with Error Bars",
    fill = "Mean LB Score"
  )
```

## TreeShrink

```bash
python3 /usr/local/src/TreeShrink-1.3.8b/run_treeshrink.py -t iqtree_bestparalog.treefile > treeshrink_bestparalog.log

iqtree2 -S Trimal_paralogs_geneious --prefix iqtree_allparalogs -m GTR+G+I -B 1000 -T 32

# Run treeshrink and remove '.main' '.0' '.1' annotations in paralog trees for ASTER/A-Pro

python3 /usr/local/src/TreeShrink-1.3.8b/run_treeshrink.py -t iqtree_allparalogs.treefile > treeshrink_allparalogs.log
sed 's/.main//g' < treeshrink_allparalogs_output.treefile | sed 's/\.\d:/:/g' > treeshrink_iqtree_allparalogs2.treefile

```

## IQtree concat (GTRG, GHOST models)
```bash
iqtree2 -s Concatenated_alignments_bestparalog.fasta --prefix GTRGI -m GTR+G+I -T 30
iqtree2 -s Concatenated_alignments_bestparalog.fasta --prefix GTRGI -m MFP -T 30
iqtree2 -s Concatenated_alignments_bestparalog.fasta --prefix GTRGI -m GTR+FO*H2 -T 30
iqtree2 -s Concatenated_alignments_bestparalog.fasta --prefix GTRGI -m GTR+FO*H4 -T 30
iqtree2 -s Concatenated_alignments_bestparalog.fasta --prefix GTRGI -m GTR+FO*H6 -T 30
iqtree2 -s Concatenated_alignments_bestparalog.fasta --prefix GTRGI -m GTR+FO*H8 -T 30

# Choose lowest BIC score + highest BIC weight

```


## Astral and Astral Pro
```bash
conda activate aster

astral4 -t 20 -o bestparalog_astral.tre -i bestparalog_treeshrink_output.treefile 2>LOG_FILE
```

## Plastid analyses & tree: GTRGI, Ghost heterotachy model (Crotty et al., 2020) for 2-8 rate classes
```bash
iqtree2 -s plastid_concat.fasta --prefix GTRGI -m GTR+G+I -T 30
iqtree2 -s plastid_concat.fasta --prefix GTRGI -m MFP -T 30
iqtree2 -s plastid_concat.fasta --prefix GTRGI -m GTR+FO*H2 -T 30
iqtree2 -s plastid_concat.fasta --prefix GTRGI -m GTR+FO*H4 -T 30
iqtree2 -s plastid_concat.fasta --prefix GTRGI -m GTR+FO*H6 -T 30
iqtree2 -s plastid_concat.fasta --prefix GTRGI -m GTR+FO*H8 -T 30

# Choose lowest BIC score + highest BIC weight

```

## gCF, sCF, plotting

```bash
iqtree -t bestparalog_astral.tre --gcf bestparalog_treeshrink_output.treefile -s Concatenated_alignments_bestparalog.fasta --scf 10000 --prefix concord_bestparalog
```

```{r}
## Visualize by plotting gCF vs sCF (gene and site concordance factors)

library(viridis)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(GGally)
library(entropy)

# read the data
d = read.delim("concord.cf.stat", header = T, comment.char=‘#')
               

# plot the values
ggplot(d, aes(x = gCF, y = sCF)) + 
    geom_point(aes(colour = Label)) + 
    scale_colour_viridis(direction = -1) + 
    xlim(0, 100) +
    ylim(0, 100) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_text_repel()
	
plot(ladderize(tre))
nodelabels(frame="none")  ## To reference the node numbers


## More exploration of concordance factors (http://www.robertlanfear.com/blog/files/concordance_factors.html)
```

## Topology tests

## PhyTop and plotting with treeio

## Phylonet, Phylonetworks

```bash
## Preparing gene tree file and running Phylonet in R

See the tutorial and citations therein (https://phylogenomics.rice.edu/html/phylonetTutorial.html)

## 0. Load libraries in R
library(phytools)

## 1. in R, read in gene trees as a multiphylo object
tre <- read.tree("treeshrink_bestparalog.treefile")

## 2. check to see that it is a multiphylo object
class(tre)
# "multiphylo"

## 3. specify tips to be pruned from tree
taxa_to_drop<_c("fastp_Coelia_triptera_S10","fastp_Yoania_prainii_S14","fastp_Yoania_japonica_S12","fastp_Calypso_bulbosa_var_americana_S4","fastp_Calypso_bulbosa_Asia_S37","fastp_Changnienia_amoena_S44","fastp_Tipularia_japonica_S45","fastp_Corallorhiza_maculata_var_maculata2_S36","fastp_Corallorhiza_mertensiana_S30","fastp_Corallorhiza_maculata_var_mexicana_S17","fastp_Corallorhiza_macrantha_S21","fastp_Corallorhiza_bulbosa_S19","fastp_Corallorhiza_odontorhiza_var_odontorhiza2_S38","fastp_Corallorhiza_odontorhiza_var_pringlei_S25","fastp_Corallorhiza_odontorhiza_Mexico_S47","fastp_Corallorhiza_wisteriana_Eastern_US_S27","fastp_Corallorhiza_striata_var_striata_S5","fastp_Corallorhiza_striata_Sierra_Nevada_S9","fastp_Corallorhiza_striata_CA_Coast_Ranges_S11","fastp_Corallorhiza_involuta_S1","fastp_Oreorchis_indica2_S46","fastp_Cremastra_variabilis_S24","fastp_Cremastra_aphylla_S18","fastp_Cremastra_saprophytica_S22","fastp_Cremastra_unguiculata_S43","fastp_Aplectrum_hyemale_S8","fastp_Govenia_superba_S6","fastp_Govenia_capitata_S42","fastp_Dactylostalix_ringens_S39","fastp_Ephippianthus_schmidtii_S26","fastp_Ephippianthus_sawadanus_S31","fastp_Brassavola_glauca_S40")

## 4. drop all tips except the ones you want to keep
pruned.tree<-drop.tip.multiPhylo(tre,taxa_to_drop)

## 5. plot the first gene tree to see if it worked
plot(pruned.tree[[1]])

## 6. write the tree to file
write.tree(pruned.tree,file="treeshrink_Calyps_subset.nex")

## 7. You need to convert the tree to nexus format, either manually or via phytools:

writeNexus(pruned.tree,file="nexustest.nex")

## 8. You'll need to edit the nexus file to look like this:

```bash

#NEXUS
 
BEGIN TREES;


Tree gt001=(((fastp_Oreorchis_coreana_S32:0.0176422487,fastp_Cremastra_appendiculata_S20:0.0085695487):1e-06,(((fastp_Oreorchis_indica1_S35:2e-06,fastp_Oreorchis_fargesii_S34:1e-06):1e-06,fastp_Oreorchis_bilamellata_S48:1e-06):2.0377e-06,(fastp_Oreorchis_erythrochrysea_S33:2e-06,fastp_Oreorchis_patens_S29:1e-06):1e-06):2.0306e-06):1e-06,(((fastp_Corallorhiza_bentleyi_S3:2.5214e-06,(fastp_Corallorhiza_striata_var_vreelandii_S7:3.0243e-06,fastp_Corallorhiza_trifida_S2:1e-06):1e-06):0.008596816,fastp_Corallorhiza_wisteriana_Western_US_S28:0.0085942908):2.0473e-06,fastp_Corallorhiza_maculata_var_occidentalis_S15:2e-06):0.0086608932);
Tree gt002=(((((fastp_Corallorhiza_bentleyi_S3:0.0037893074,((fastp_Corallorhiza_wisteriana_Western_US_S28:0.0096400708,fastp_Corallorhiza_maculata_var_occidentalis_S15:0.0018930636):0.0018905778,(fastp_Corallorhiza_striata_var_vreelandii_S7:0.0018965563,fastp_Corallorhiza_trifida_S2:0.0076179904):0.0018905325):1.001e-06):0.001890477,fastp_Oreorchis_coreana_S32:0.0038012431):1.001e-06,((fastp_Oreorchis_erythrochrysea_S33:0.0038074333,fastp_Oreorchis_indica1_S35:2.002e-06):0.0018910937,fastp_Oreorchis_patens_S29:0.0037836267):1.001e-06):1.001e-06,(fastp_Oreorchis_bilamellata_S48:1.001e-06,fastp_Oreorchis_fargesii_S34:0.0018865137):0.0018889098):0.0037844083,fastp_Cremastra_appendiculata_S20:0.0037834846);
Tree gt003=((((fastp_Oreorchis_bilamellata_S48:1.1108e-06,fastp_Oreorchis_coreana_S32:1.1108e-06):0.0041928447,fastp_Oreorchis_fargesii_S34:1.1108e-06):0.0190836901,(fastp_Corallorhiza_bentleyi_S3:0.0305309322,(fastp_Corallorhiza_trifida_S2:0.0126344377,(fastp_Oreorchis_erythrochrysea_S33:0.0045712793,fastp_Oreorchis_indica1_S35:0.0083840509):0.0041861576):1.1108e-06):0.0049402476):0.0167464358,fastp_Cremastra_appendiculata_S20:0.0219480731);
Tree gt004=(fastp_Cremastra_appendiculata_S20:0.0030859065,((((fastp_Corallorhiza_bentleyi_S3:1.0004e-06,fastp_Corallorhiza_striata_var_vreelandii_S7:0.0060026864):0.0029992081,fastp_Corallorhiza_trifida_S2:0.0069111331):0.0019470118,(((fastp_Oreorchis_bilamellata_S48:1.0004e-06,fastp_Oreorchis_fargesii_S34:0.0034964139):0.0035081263,fastp_Oreorchis_coreana_S32:0.0117322654):1.0004e-06,fastp_Oreorchis_patens_S29:1.0004e-06):0.0038676149):1.0004e-06,fastp_Oreorchis_indica1_S35:0.003441931):0.002733112);
Tree gt005=(((((fastp_Cremastra_appendiculata_S20:0.0146515673,((fastp_Corallorhiza_bentleyi_S3:0.00368338,(fastp_Oreorchis_bilamellata_S48:0.0057207567,((fastp_Oreorchis_coreana_S32:1.0322e-06,fastp_Oreorchis_fargesii_S34:1.0322e-06):1.0322e-06,fastp_Oreorchis_patens_S29:1.0322e-06):1.0322e-06):0.0074283596):1.0322e-06,(fastp_Corallorhiza_trifida_S2:1.0322e-06,fastp_Oreorchis_erythrochrysea_S33:1.0322e-06):2.0644e-06):1.0322e-06):0.0035918253,fastp_Oreorchis_indica1_S35:1.0322e-06):3.4491e-06,fastp_Corallorhiza_maculata_var_occidentalis_S15:1.0322e-06):1.0322e-06,fastp_Corallorhiza_wisteriana_Western_US_S28:1.0322e-06):3.7869e-06,fastp_Corallorhiza_striata_var_vreelandii_S7:1.0322e-06);

...the rest of the gene trees...

END;
 
BEGIN PHYLONET;
 
InferNetwork_MPL (all) 0 -pl 30;
 
END;



## 9. In the above, you need the PHYLONET block, with this general format

```bash

END;
 
BEGIN PHYLONET;

Algorithm (taxa_to_analyze) #hybridizations -pl <#threads>

END;



## 10. Now, you are ready to run PHYLONET. These analyses for each H-value (# hybridizations) take 10-30 minutes.
## 11. Running the block above specifies zero hybridizations, and essentially finds the "species tree"

java -jar /usr/local/src/PhyloNet.jar treeshrink_Calyps_subset.nex

## 12. When finished, you'll get a bunch of output. Take the first of the final five trees and the likelihood score to calculate the AIC

## 13. Repeat the analysis for however many numbers of H (0,1,2,3,4,5,...)

## 14. Save the likelihood scores as a column in excel (or do this in R).

## 15. The number of parameters (k) = the number of total internal + terminal branches in the tree PLUS the number of pre-specified hybridization events (k + H). Use these for AIC calcs, where AIC = 2*k - 2*lnL.

## 16. Calulate AIC scores in excel with " =((2*C3)-(2*LN(B3)))" where the likelihoodis in cell B3 and k is in B4.

## 17. Calculate the delta AIC (AIC for value of H - minAIC) in excel.

## 18. Calculate AIC weights (wAIC) as "=EXP(-0.5*E3)" where the deltaAIC is in cell E3, and drag down.

# Voila! 

## The H-values (# of hybridizations) with the lowest AIC and highest wAIC is the optimal. It could very well be zero.
```

# II. Analysis with PhyloNetworks
```bash
Following (https://crsl4.github.io/PhyloNetworks.jl/dev/man/snaq_plot/)

## 1. Start julia and add phylonetworks
julia
using PhyloNetworks

## 2. Read in gene trees, view tree #3

genetrees = readMultiTopology("calyps_phylonetworks.tre");
genetrees[3]

## 3. Load phyloplots and plot tree #3

using PhyloPlots
plot(genetrees[3]); # tree for 3rd gene

## 4. Calculate quartet Concordance Factors
q,t = countquartetsintrees(genetrees);

## 5. Write CFs to a table to save and view
```bash
using CSV
df = writeTableCF(q,t)   # data frame with observed CFs: gene frequencies
CSV.write("tableCF.csv", df); # to save the data frame to a file
raxmlCF = readTableCF("tableCF.csv") # read in the file and produces a "DataCF" object
less("tableCF.csv")

raxmlCF = readTrees2CF(genetrees, whichQ="rand", numQ=200, CFfile="tableCF10.txt")
## 6. Get a starting tree -- in this case, the astral "species tree"
astraltree = readTopology("astral.tre")

## 7. For multithreading/parallel jobs
using Distributed
addprocs(30)
@everywhere using PhyloNetworks

## 8. Run the first network analysis with H=0, essentially find the "species tree"
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234)

## 9. Got booted out of julia for some reason, needed to start over!
using PhyloPlots
q,t = countquartetsintrees(genetrees); # read in trees, calculate quartet CFs
df = writeTableCF(q,t)

using CSV
CSV.write("tableCF.csv", df);
raxmlCF = readTableCF("tableCF.csv")
raxmlCF = readTrees2CF(genetrees, whichQ="rand", numQ=200, CFfile="tableCF10.txt")

astraltree = readTopology("astral.tre")
plot(astraltree, showedgelength=true);

## 10 Run the first network analysis with H=0, essentially find the "species tree"
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234)

## 11. Copy the log-likelihood from the screen output after run finishes, or get from "net0.out" file

## 12. Run the second network analysis with H=1, allowing 1 hybridization event, using 'net0' as starting tree

net1 = snaq!(net0, raxmlCF, hmax=1, filename="net1", seed=2345) # this runs for ~30 min

## 13. Check out the output

plot(net1, showgamma=true);
less("net1.err") # would provide info about errors, if any
less("net1.out") # main output file with the estimated network from each run
less("net1.networks") # extra info
net1

## 14. Save the network file

writeTopology(net1, round=true, digits=2)

## 15. Run the third network analysis with H=2, allowing 2 hybridization events, using 'net0' as starting tree

net2 = snaq!(net0,raxmlCF, hmax=2, filename="net2", seed=3456)
plot(net2, showgamma=true);

## 16. Run the rest of the analyses analysis with H=3-5, allowing 3-5 hybridization events, using 'net0' as starting tree
### Save the likelihoods for AIC calcs


net3 = snaq!(net0,raxmlCF, hmax=3, filename="net3", seed=4567)
net4 = snaq!(net0,raxmlCF, hmax=4, filename="net4", seed=1437)
net5 = snaq!(net0,raxmlCF, hmax=5, filename="net5", seed=8701)
```


## Divtime with LSD in IQtree2
```bash
iqtree2 -s newtargets_concat_calyps.fasta --date lsd_calib2.txt -te concat_h6.treefile -m GTR -o fastp-Brassavola-glauca_S40 -T 32 --date-tip 0 --date-root -60:-40 -u 1 --date-ci 1000 --redo

#Open/edit with FigTree (http://tree.bio.ed.ac.uk/software/figtree/)

```

## Stochastic character mapping
```{r}
## Ancestral state reconstruction of Corallorhiza/Oreorchis floral characters

setwd("H:/0001_research/0001_manuscripts/2022_Calypsoinae/2024_03_22_newtargets")

library(phytools)
library(geiger)

# Charater state definitions
# ltype 	 labellum type: 	0 = thin, flat, white/yellow/cream/orange; 1= thick, boat-shaped, red/purple/yellow, striped
# lobe  	 labellum lobes:	0 = absent; 1 = present
# callus	 fused lamellae:	0 = absent; 1 = present
# mentum	 nectar spur:		0 = absent; 1 = present

# Read in tree. Chose to use the concatenated tree under GHOST H6 heterotachy model as the most accurate representation of branch lengths
tre <- read.tree("concat_h6_boot.contree")

# Drop tips not of interest, set tip labels from tree, save as csv, then add characters/states to csv and re-import

tre2<-drop.tip(tre,c("fastp-Brassavola-glauca_S40", "fastp-Aplectrum-hyemale_S8", "fastp-Cremastra-aphylla_S18", "fastp-Cremastra-saprophytica_S22", "fastp-Cremastra-appendiculata_S20", "fastp-Cremastra-variabilis_S24", "fastp-Cremastra-unguiculata_S43", "fastp-Govenia-capitata_S42", "fastp-Govenia-superba_S6", "fastp-Dactylostalix-ringens_S39", "fastp-Ephippianthus-sawadanus_S31", "fastp-Ephippianthus-schmidtii_S26", "fastp-Calypso-bulbosa-Asia_S37", "fastp-Calypso-bulbosa-var-americana_S4", "fastp-Changnienia-amoena_S44", "fastp-Tipularia-japonica_S45", "fastp-Coelia-triptera_S10", "fastp-Yoania-japonica_S12", "fastp-Yoania-prainii_S14"))

labels <- tre2$tip.label
write.csv(labels,file="tiplabels.csv")

# Have a look at the tree
plot(ladderize(tre2))

# Add a new tip to the tree

tre3<-bind.tip(tre2, "Corallorhiza_sinensis", edge.length=0.05, where=28, position=0.5)


# Use penalized likelihood to make the tree ultrametric
tre3<-chronos(tre2)

# Have a look
plot(ladderize(tre3))

# preserve branching order as ladderized
tre4<-ladderize(tre3)

# Have a look
plot(tre4)

# Read in the character data, with column 1 = to the tip labels, and 2-5 as the characters/states
tips <-read.csv("tiplabels.csv",header=T,row.names = 1)

# Define characters as vectors
ltype<-as.matrix(tips)[,1]
lobe<-as.matrix(tips)[,2]
callus<-as.matrix(tips)[,3]
mentum<-as.matrix(tips)[,4]
troph<-as.matrix(tips)[,6]

# Use an MCMC implementation of stochastic character mapping with phytools

models <- list("ER","ARD","SYM","loss")

# where ordered is...

loss<-matrix(c(0,0,0,1,0,0,1,1,0),nrow=3)

mtrees.ltype<-make.simmap(tre9,ltype,model=models,nsim=1000)
mtrees.lobe<-make.simmap(tre9,lobe,model=models,nsim=1000)
mtrees.callus<-make.simmap(tre9,callus,model=models,nsim=1000)
mtrees.mentum<-make.simmap(tre9,mentum,model=models,nsim=1000)
mtrees.troph<-make.simmap(tre9,troph,model=models,nsim=1000)

pd.ltype<-summary(mtrees.ltype,plot=FALSE)
pd.lobe<-summary(mtrees.lobe,plot=FALSE)
pd.callus<-summary(mtrees.callus,plot=FALSE)
pd.mentum<-summary(mtrees.mentum,plot=FALSE)
pd.troph<-summary(mtrees.troph,plot=FALSE)

# Set up a 3 x 2 plot
par(mfrow = c(3, 2))

#plotting

cols = c("lightblue", "yellow", "red")

plot(pd.ltype,fsize=0.8,ftype="i")
title("Labellum type (thin 0, thick 1)")
# add.simmap.legend(colors=cols)

plot(pd.lobe,fsize=0.8,ftype="i")
title("Labellum lobes (absent 0, present 1)")
# add.simmap.legend(colors=cols)

plot(pd.callus,fsize=0.8,ftype="i")
title("Callus (absent 0, present 1)")
# add.simmap.legend(colors=cols)

plot(pd.mentum,fsize=0.8,ftype="i")
title("Mentum (absent 0, present )1")
# add.simmap.legend(colors=cols)

plot(pd.troph,fsize=0.8,ftype="i")
title("Trophic mode (leafy 0, PM 1, HM 1")
## Add legend manually at the end
add.simmap.legend(colors=cols,prompt=FALSE)

#################

## Updates 5-state analysis

# character definitions
1 - Leafy, presumably autotrophic (based on leaves present)
2 - Leafy, partially mycoheterotrophic (based on isotopic evidence)
3 - Leafless, presumably partially mycoheterotrophic (based on retention of all plastid photosynthetic genes but a lack of isotopic evidence)
4 - Leafless, partially mycoheterotrophic (based on isotopic evidence)
5 - Leafless, holomycotrophic (based on losses of several plastid photosynthetic genes; isotopic evidence for C. striata and C. maculata but not C. sinensis)

# transition matrix
loss3<-matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow = 5, ncol = 5)
loss3
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    1    1    1    1
[2,]    0    0    1    1    1
[3,]    0    0    0    0    1
[4,]    0    0    0    0    1
[5,]    0    0    0    0    0

# simmap & summarize on tree
mtrees.troph<-make.simmap(tre2,troph,model=loss3,nsim=100)
pd.troph<-summary(mtrees.troph,plot=TRUE)
 
```

















