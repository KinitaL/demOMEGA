#! /usr/bin/Rscript
# First of all, we clean the workspace
rm(list=ls())

# Get the output dir
args <- commandArgs(trailingOnly = FALSE)
dir <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))[2]

# install libraries
install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("ggtree")

library(ggplot2)
library(ggtree)
library(ape)
library(stringr)
library(dplyr)

# Read raw tree text from file
tree_text <- readLines(paste(dir, "tree.treefile", sep="/"))
tree_text <- paste(tree_text, collapse = "")  # in case it's multi-line

# Parse tree with ape
tree <- read.tree(text = tree_text)

# Extract dN/dS values from raw Newick string
matches <- str_match_all(tree_text, "([A-Za-z0-9]+)?#([0-9\\.Ee+-]+)")[[1]]
dn_ds_values <- as.numeric(matches[,3])
labels <- matches[,2]  # NA for internal nodes

# Plot tree with ggtree
p <- ggtree(tree)

# Add empty column for dnds
p$data$dnds <- NA_real_

# Fill tips by matching label
for (i in seq_along(labels)) {
  if (!is.na(labels[i])) {
    idx <- which(str_detect(p$data$label, fixed(labels[i])))
    if (length(idx) > 0) {
      p$data$dnds[idx] <- dn_ds_values[i]
    }
  }
}

# Fill internal nodes in order
internal_nodes <- which(is.na(p$data$label) | p$data$label == "" | str_detect(p$data$label, "^#"))
p$data$dnds[internal_nodes] <- dn_ds_values[is.na(labels)]

# Check mapping
print(p$data[, c("label", "dnds")])

# Plot with branch width reflecting dN/dS
p_final <- p +
  geom_tree(aes(size = dnds, color = dnds)) +
  geom_text2(aes(label = round(dnds, 3))) +
  geom_tiplab(align = TRUE, linetype = "dotted", size = 3) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_color_viridis_c(option = "plasma") +
  theme_tree2() +
  ggtitle("Branch width, color and labels reflect dN/dS values")

print(p_final)

ggsave(paste(dir, "tree.png", sep="/"), plot = p_final, width = 8, height = 6, dpi = 300)
