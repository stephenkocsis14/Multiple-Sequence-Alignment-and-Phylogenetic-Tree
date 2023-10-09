## Setup environment
knitr::opts_chunk$set(echo = TRUE)

## Load packages
library(Biostrings)
library(DECIPHER)
library(msa)
library(seqinr)
library(ggseqlogo)
library(ggmsa)
library(ape)

## Read sequences from FASTA file of the cases
sequence <- readAAStringSet("homologs_seq.FASTA")

## Read the FASTA file
mySequences <- readAAStringSet("homologs_seq.fasta")
mySequences

## Now that we have loaded the sequences, we can run the msa() function which, by default, runs ClustalW with default parameters:
myFirstAlignment <- msa(mySequences)
myFirstAlignment

## Obviously, the default printing function shortens the alignment for the sake of compact output. The print() function provided by the msa package provides some ways for customizing the output, such as, showing the entire alignment split over multiple blocks of sub-sequences:
print(myFirstAlignment, show="complete")

## Write the fasta output of the aligned sequence

# Extract the sequences from the MSA
sequences <- as.character(myFirstAlignment)

# Define the sequence names (headers)
sequence_names <- names(sequences)

# Specify the output file path
output_file <- "output_alignment.fasta"

# Write the sequences to a FASTA file
write.fasta(sequences, names = sequence_names, file.out = output_file)

## Calculate conservation statistics from the alignment
myAlignment <- myFirstAlignment

# Assuming 'myAlignment' is your MsaAAMultipleAlignment object
# Convert it to an AAStringSet
alignment <- AAStringSet(myAlignment)

# Calculate conservation scores
conservation_scores <- consensusMatrix(alignment)

# Calculate summary statistics
mean_score <- mean(conservation_scores)
median_score <- median(conservation_scores)
min_score <- min(conservation_scores)
max_score <- max(conservation_scores)

# Create a custom summary report
cat("Conservation Score Summary Report:\n")
cat("Mean Score:", mean_score, "\n")
cat("Median Score:", median_score, "\n")
cat("Minimum Score:", min_score, "\n")
cat("Maximum Score:", max_score, "\n")

## Visualize the multiple sequence alignment
aln <- myFirstAlignment
class(aln) <- "AAMultipleAlignment"

ggmsa::ggmsa(aln, 
             start = 0, 
             end = 50,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 51, 
             end = 100,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 101, 
             end = 150,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 151, 
             end = 200,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 201, 
             end = 250,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 251, 
             end = 300,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 301, 
             end = 350,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 351, 
             end = 400,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 401, 
             end = 450,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

ggmsa::ggmsa(aln, 
             start = 451, 
             end = 492,
             seq_name = TRUE,
             char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") + 
  geom_msaBar()

## Starting phylogenetic tree analysis
## Compute distance matrix
my_alignment_sequence <- msaConvert(myFirstAlignment, type="seqinr::alignment")
distance_alignment <- dist.alignment(my_alignment_sequence)

## compute phylogenetic tree using neighbor joining
Tree <- ape::bionj(distance_alignment)

## Change the tip.labels to include the protein and organism names
Tree$tip.label <- c("f3a-201 [Gasterosteus aculeatus]", "KAF1372746 [Perca fluviatilis]", "TDG97026 [Perca flavescens]", "TKS88104 [Collichthys lucidus]", "XP_010738138 [Larimichthys crocea]", "XP_030254448 [Sparus aurata]", "XP_029314020 [Cottoperca gobio]", "ACQ58994 [Anoplopoma fimbria]", "f3a-201 [Lates calcarifer]", "XP_018554310 [Lates calcarifer]", "XP_026147516 [Mastacembelus armatus]", "f3a-201 [Amphiprion percula]", "XP_023155773 [Amphiprion ocellaris]", "XP_022069124 [Acanthochromis polyacanthus]", "XP_008303090 [Stegastes partitus]", "XP_028288293 [Parambassis ranga]", "FGENESH CD Start 2148288", "XP_008305582 [Cynoglossus semilaevis]", "CDQ68572 [Oncorhynchus mykiss]")

## Plot the tree
plot(Tree, show.tip.label = TRUE, cex = 0.8, edge.width = 2)
