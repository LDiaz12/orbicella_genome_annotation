## The following code was developed by Laurel C. Diaz on 6/22/2026
## This code is for the genome annotation for Orbicella corals 

## Load libraries 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(here)
library(Biostrings)
library(rtracklayer)
library(biomartr)

#### Read in O franksi files and sort chromosomes by size to create chromosomal reference map ####
### Read in the softmasked fasta file of O franksi annotation (from ASG portal) ###
### Read in gff file for O franksi - this has the structural annotations: "coordinates" of where things are ###
o_fr_fasta <- readDNAStringSet(here("Data", "O_franksi", "O_franksi_softmasked.fa"))
o_fr_gff <- import(here("Data", "O_franksi", "O_franksi_genes.gff3"))

# Sort FASTA file by size: largest to smallest
o_fr_fasta_sorted <- o_fr_fasta[order(width(o_fr_fasta), decreasing = TRUE)] # reorder chromosomes by size largest to smallest
# Clean the names (remove any other symbols)
names(o_fr_fasta_sorted) <- gsub(" .*", "", names(o_fr_fasta_sorted))

# Sort gff3 to match FASTA order 
o_fr_chrom_order <- names(o_fr_fasta_sorted)
o_fr_gff_sorted <- o_fr_gff[order(match(as.character(seqnames(o_fr_gff)), o_fr_chrom_order))]

# Filter to top 15 chromosomes 
top15_ofr <- names(o_fr_fasta_sorted)[1:15]
o_fr_fasta_chrom <- o_fr_fasta_sorted[1:15]
o_fr_gff_chrom <- o_fr_gff_sorted[as.character(seqnames(o_fr_gff_sorted)) %in% top15_ofr]

# Confirm
length(o_fr_fasta_chrom)
sum(width(o_fr_fasta_chrom))
length(unique(seqnames(o_fr_gff_chrom)))

# Write files for sorted FASTA and GFF3
writeXStringSet(o_fr_fasta_chrom, here("Data", "O_franksi", "o_fr_sorted.fasta"))
file.info(here("Data","O_franksi", "o_fr_sorted.fasta"))$size / 1e6

export(o_fr_gff_sorted, here("Data", "O_franksi", "o_fr_sorted.gff3"))

# Create a data frame for downstream visualization 
chrom_sizes <- data.frame(
  chromosome = names(o_fr_fasta_sorted), 
  length_bp = width(o_fr_fasta_sorted)
)

# Add rank column
chrom_sizes$rank <- 1:nrow(chrom_sizes)

o_franksi_chr <- ggplot(chrom_sizes[1:20, ], aes(x = rank, y = length_bp)) + # this plot shows genome sequences by size with chromosomes 1-15
  geom_point(size = 1.5) + 
  geom_vline(xintercept = 15, color = "coral", linetype = "dashed") + 
  scale_x_continuous(breaks = 1:20) + 
  scale_y_continuous(limits = c(min(chrom_sizes$length_bp[1:20]) * 0.9,
                               max(chrom_sizes$length_bp[1:20]) * 1.1)) +
  labs(title = "O. franksi genome sequence sizes",
       x = "Sequence rank", 
       y = "Length (bp)") + 
  theme_classic()
o_franksi_chr # we should see a clear difference in size from the chromosomes and the rest of the sequences
ggsave(plot = o_franksi_chr, here("Outputs", "o_franksi_sequence_sizes.png"))

#### Read in O annularis files, sort, then get pseudochromosomes and align to reference chromosomal map ####
### Sort O annularis chromosomes by size ### 
o_ann_fasta_new <- readDNAStringSet(here("Data", "O_annularis", "o_ann_assembly_new.fna"))

length(o_ann_fasta_new) # 2767 sequences
sum(width(o_ann_fasta_new)) #551315961 genome size

o_ann_gff <- import(here("Data", "O_annularis", "o_ann_annotation.gff3"))
length(o_ann_gff) # 819084 features

names(o_ann_fasta_new) <- gsub(" .*", "", names(o_ann_fasta_new))
o_ann_fasta_new_sorted <- o_ann_fasta_new[order(width(o_ann_fasta_new), decreasing = TRUE)]

# Write data frame
chrom_sizes_o_ann_new <- data.frame(
  chromosome = names(o_ann_fasta_new_sorted),
  length_bp  = width(o_ann_fasta_new_sorted)
)
head(chrom_sizes_o_ann_new, 20)

table(cut(chrom_sizes_o_ann_new$length_bp,
          breaks = c(0, 100000, 1000000, 10000000, Inf),
          labels = c("<100kb", "100kb-1Mb", "1Mb-10Mb", ">10Mb")))
ggplot(chrom_sizes_o_ann_new, aes(x = 1:nrow(chrom_sizes_o_ann_new), y = length_bp)) +
  geom_point(size = 0.5) +
  labs(title = "O. annularis new assembly sequence sizes",
       x = "Sequence rank",
       y = "Length (bp)") +
  theme_classic()

# Write sorted fasta file 
writeXStringSet(o_ann_fasta_new_sorted, here("Data", "O_annularis", "o_ann_new_sorted.fasta"))
o_ann_fasta_new_sorted
length(o_ann_fasta_new_sorted)
head(names(o_ann_fasta_new_sorted), 5)

### After this step, then work in the terminal and use bash with conda and RagTag (which will use minimap2)
### to align O annularis assembly to O franksi 

### Bring the new files created from RagTag in the terminal back into R ###
o_ann_pseudo_new <- readDNAStringSet(here("Data", "O_annularis", "ragtag_output_oann_new2", "ragtag.scaffold.fasta"))

# Sort pseudo chromosomes by size
o_ann_pseudo_new_sorted <- o_ann_pseudo_new[order(width(o_ann_pseudo_new), decreasing = TRUE)]

# Clean the names
names(o_ann_pseudo_new_sorted) <- gsub(" .*", "", names(o_ann_pseudo_new_sorted))

# Write data frame 
chrom_sizes_oann_new <- data.frame(
  chromosome = names(o_ann_pseudo_new_sorted),
  length_bp  = width(o_ann_pseudo_new_sorted)
)
head(chrom_sizes_oann_new, 20)

table(cut(chrom_sizes_oann_new$length_bp,
          breaks = c(0, 100000, 1000000, 10000000, Inf),
          labels = c("<100kb", "100kb-1Mb", "1Mb-10Mb", ">10Mb")))

ggplot(chrom_sizes_oann_new, aes(x = 1:nrow(chrom_sizes_oann_new), y = length_bp)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 15, color = "red", linetype = "dashed") +
  scale_x_continuous(breaks = 1:20, limits = c(1,20)) +
  labs(title = "O. annularis new pseudochromosome sizes",
       x = "Sequence rank",
       y = "Length (bp)") +
  theme_classic()

# Chromosome numbers 1:15 selected
o_ann_chrom_new <- o_ann_pseudo_new_sorted[1:15]
length(o_ann_chrom_new)
sum(width(o_ann_chrom_new))

# Write pseudochromosome sorted FASTA file
writeXStringSet(o_ann_chrom_new, here("Data", "O_annularis", "o_ann_pseudochrom_new_sorted.fasta"))

# Visualize with a plot 
chrom_sizes_oann$rank <- 1:nrow(chrom_sizes_oann)

o_ann_chr_plot <- ggplot(chrom_sizes_oann[1:20, ], aes(x = rank, y = length_bp)) + # this plot shows genome sequences by size with chromosomes 1-15
  geom_point(size = 1.5) + 
  geom_vline(xintercept = 15, color = "coral", linetype = "dashed") + 
  scale_x_continuous(breaks = 1:20) + 
  scale_y_continuous(limits = c(min(chrom_sizes_oann$length_bp[1:20]) * 0.9,
                                max(chrom_sizes_oann$length_bp[1:20]) * 1.1)) +
  labs(title = "O. annularis pseudochromosome sizes",
       x = "Sequence rank", 
       y = "Length (bp)") + 
  theme_classic()
o_ann_chr_plot

ggsave(plot = o_ann_chr_plot, here("Outputs", "o_ann_pseudochr_plot.png"))

## Continue with work in the terminal using Liftoff to align 
## Read in lifted GFF3 file of O annularis 
o_ann_gff_new_pseudo <- import(here("Data", "O_annularis","o_ann_pseudochrom_new.gff3"))
length(o_ann_gff_new_pseudo)
head(unique(as.character(seqnames(o_ann_gff_new_pseudo))), 20)

# Filter top 15 pseudochromosomes
top15_oann_new <- names(o_ann_chrom_new)
o_ann_gff_chrom_new <- o_ann_gff_new_pseudo[as.character(seqnames(o_ann_gff_new_pseudo)) %in% top15_oann_new]

# Confirm 
length(unique(seqnames(o_ann_gff_chrom_new))) # 15 

# Export final GFF3 file for O annularis
export(o_ann_gff_chrom_new, here("Data", "O_annularis", "o_ann_pseudochrom_final_new.gff3"))

#### Read in O faveolata files, sort by size, then save file to work in the terminal ####
# O faveolata FASTA and GFF files are from the publication: Young et al 2024 "Annotated genome and 
# transcriptome of the endangered Caribbean mountainous star coral (Orbicella faveolata) using 
# PacBio long-read sequencing" 
o_fav_pacbio <- readDNAStringSet(here("Data", "O_faveolata", "o_fav_scaffolds.fa"))
length(o_fav_pacbio)
sum(width(o_fav_pacbio))

o_fav_pacbio_gff <- import(here("Data", "O_faveolata", "o_fav.gff3"))
head(names(o_fav_pacbio), 10)
head(unique(as.character(seqnames(o_fav_pacbio_gff))), 10)

intersect(names(o_fav_pacbio), 
          unique(as.character(seqnames(o_fav_pacbio_gff))))

o_fav_pacbio_sorted <- o_fav_pacbio[order(width(o_fav_pacbio), decreasing = TRUE)]
o_fav_chrom_sizes_pacbio <- data.frame(
  chromosome = names(o_fav_pacbio_sorted),
  length_bp  = width(o_fav_pacbio_sorted)
)

head(o_fav_chrom_sizes_pacbio, 20)

table(cut(o_fav_chrom_sizes_pacbio$length_bp,
          breaks = c(0, 100000, 1000000, 10000000, Inf),
          labels = c("<100kb", "100kb-1Mb", "1Mb-10Mb", ">10Mb")))

# Plot to visualize
# Create a rank column 
o_fav_chrom_sizes_pacbio$rank <- 1:nrow(o_fav_chrom_sizes_pacbio)
o_fav_plot <- ggplot(o_fav_chrom_sizes_pacbio, aes(x = rank, y = length_bp)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 15, color = "coral", linetype = "dashed") +
  scale_x_continuous(breaks = 1:51) +
  labs(title = "O. faveolata PacBio assembly sequence sizes",
       x = "Sequence rank",
       y = "Length (bp)") +
  theme_classic()
o_fav_plot
ggsave(plot = o_fav_plot, here("Outputs", "o_fav_chr_plot.png"))

# Get top 15 scaffold names
top15_pacbio <- names(o_fav_pacbio_chrom)
# Filter GFF3 to top 15
o_fav_pacbio_gff_chrom <- o_fav_pacbio_gff[as.character(seqnames(o_fav_pacbio_gff)) %in% top15_pacbio]
# Confirm
length(unique(seqnames(o_fav_pacbio_gff_chrom)))  # should be 15
# Export final files
writeXStringSet(o_fav_pacbio_chrom, here("Data", "O_faveolata", "o_faveolata_pacbio_15chr.fasta"))
export(o_fav_pacbio_gff_chrom, here("Data", "O_faveolata", "o_faveolata_pacbio_15chr.gff3"))

### Repeat process for Cyphastrea salae, which will serve as our outgroup 
## Load in Cyphastrea FASTA file 
cyph_fasta <- readDNAStringSet(here("Data", "C_salae", "cyphastrea_softmasked.fa"))
cyph_gff <- import(here("Data", "C_salae", "cyphastrea.gff3"))

# Check assembly stats 
length(cyph_fasta)
sum(width(cyph_fasta))
length(cyph_gff)
head(unique(as.character(seqnames(cyph_gff))), 10)

# Sort FASTA file by size 
cyph_fasta_sorted <- cyph_fasta[order(width(cyph_fasta), decreasing = TRUE)]

# Clean the names in the file 
names(cyph_fasta_sorted) <- gsub(" .*", "", names(cyph_fasta_sorted))

# Create data frame of chromosome sizes for visualization
cyph_chrom_sizes <- data.frame(
  chromosome = names(cyph_fasta_sorted), 
  length_bp = width(cyph_fasta_sorted)
)

# Sort gff file according to chrom order from sorted FASTA file 
cyph_chrom_order <- names(cyph_fasta_sorted)
cyph_gff_sorted <- cyph_gff[order(match(as.character(seqnames(cyph_gff)), cyph_chrom_order))]

top15_csal <- names(cyph_fasta_sorted)[1:15]
c_sal_fasta_chrom <- cyph_fasta_sorted[1:15]
c_sal_gff_chrom <- cyph_gff_sorted[as.character(seqnames(cyph_gff_sorted)) %in% top15_csal]

# Write sorted FASTA file 
writeXStringSet(cyph_fasta_sorted, here("Data", "C_salae", "cyphastrea_sorted.fasta"))
# Write sorted gff file 650
export(cyph_gff_sorted, here("Data", "C_salae", "cyphastrea_sorted.gff3"))

# Add rank column
cyph_chrom_sizes$rank <- 1:nrow(cyph_chrom_sizes)

cyphastrea_chr <- ggplot(cyph_chrom_sizes[1:20, ], aes(x = rank, y = length_bp)) + # this plot shows genome sequences by size with chromosomes 1-15
  geom_point(size = 1.5) + 
  geom_vline(xintercept = 15, color = "coral", linetype = "dashed") + 
  scale_x_continuous(breaks = 1:20) + 
  scale_y_continuous(limits = c(min(cyph_chrom_sizes$length_bp[1:20]) * 0.9,
                                max(cyph_chrom_sizes$length_bp[1:20]) * 1.1)) +
  labs(title = "C. salae genome sequence sizes",
       x = "Sequence rank", 
       y = "Length (bp)") + 
  theme_classic()
cyphastrea_chr # we should see a clear difference in size from the chromosomes and the rest of the sequences
ggsave(plot = cyphastrea_chr, here("Outputs", "cyphastrea_sequence_sizes.png"))