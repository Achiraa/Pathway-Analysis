setwd("D:/geo_data/Stromal")

# Load required packages
library(gprofiler2)
library(ggplot2)

# Read in your CSV file with the KEGG pathway analysis results
data <- read.csv("Significant stromal genes GFP.csv")

all_genes=read.csv("AllstromalGenesGFP.csv")

# Perform KEGG pathway analysis using gprofiler2
GFP_gost <- gost(
  query = data$X,
  organism = "mmusculus",
  ordered_query = FALSE,
  exclude_iea = FALSE,
  custom_bg = all_genes$X,
  sources = "KEGG",
  user_threshold = 0.05
)

GFP_gost

# Extract all the significantly enriched pathways
top_pathways <- GFP_gost$result[GFP_gost$result$significant,]

# Create the dotplot using ggplot2
ggplot(top_pathways, aes(x = reorder(term_name, -p_value), 
                         y = p_value, 
                         size = intersection_size, 
                         color = source)) +
  geom_point() +
  coord_flip() +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Stromal GFP-KEGG Pathways",
       x = "KEGG Pathway", 
       y = "p-value") +
  ylim(0,0.1)+
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
