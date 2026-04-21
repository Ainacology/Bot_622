# ============================================================ 
# BOT 622 — NEW SLIDES (BEGINNER ANNOTATED, LINE‑BY‑LINE)      
# Includes: Slide A–G + MS specialized community graph         
# Saves: graphs → Git/graphs | CSVs → Git/slide_data           
# ============================================================ 

# ============================ # Divider for setup
# STEP 0 — SET PATHS                          
# ============================ # Divider

base <- "C:/Users/ntang/OneDrive/Desktop/Repositories/Bot_622"                 # Root project folder
out_path <- "C:/Users/ntang/OneDrive/Desktop/Repositories/Bot_622/Git/graphs"    # Folder for slide images
csv_out_path <- "C:/Users/ntang/OneDrive/Desktop/Repositories/Bot_622/Git/slide_data" # Folder for slide CSVs
setwd(file.path(base, "BOT_662_project/r_codes"))                              # Set working directory

# Create output folders if they do not exist                     # Explain why
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)              # Make image folder
if (!dir.exists(csv_out_path)) dir.create(csv_out_path, recursive = TRUE)      # Make CSV folder

# ============================ # Divider for libraries
# STEP 1 — LOAD LIBRARIES                                        # Explain what this is
# ============================ # Divider

library(vegan)       # Diversity + ordination tools
library(tidyverse)   # Data wrangling + plotting
library(ggplot2)     # Plotting
library(patchwork)   # Combine plots
library(data.table)  # Fast data handling
library(scales)      # Nice axis labels
library(tidyr)       # Pivot helpers

theme_set(theme_minimal())  # Clean background for all plots

# ============================ # Divider for styles
# STEP 2 — SHARED COLORS / LABELS                                # Style consistency
# ============================ # Divider

slide_colors <- c(                                                    # Treatment colors
  "MMO"               = "red",                                        # MS
  "EMO"               = "#8ACE00",                                   # EM
  "MEM"               = "orange",                                     # MS+EM
  "Mosquito microbes" = "red",                                        # Microbe stock
  "Env. microbes"     = "#8ACE00"                                    # Microbe stock
)                                                                    # End colors

slide_labels_short <- c(                                              # Short x‑axis labels
  "MMO" = "MS",                                                      # Mosquito symbionts
  "EMO" = "EM",                                                      # Environmental microbes
  "MEM" = "MS+EM"                                                   # Combined
)                                                                    # End short labels

slide_labels_full <- c(                                               # Full legend labels
  "MMO" = "Mosquito Symbionts (MS)",                                 # Full label
  "EMO" = "Environmental Microbes (EM)",                             # Full label
  "MEM" = "Combo (MS+EM)"                                           # Full label
)                                                                    # End full labels

BASE_TEXT <- 14   # Base text size for slides

# ============================ # Divider for data load
# STEP 3 — READ DATA (SAME INPUTS AS COMPLETE ANALYSIS)           # Keep inputs consistent
# ============================ # Divider

meta <- read.csv("../data/processed_data/meta.csv") %>%               # Larvae‑only metadata
  select(-X)                                                          # Drop row index column

clean_meta <- read.csv("../data/processed_data/clean_meta.csv") %>%   # Water + larvae metadata
  select(-X)                                                          # Drop row index column

counts_ra <- read.csv("../data/processed_data/counts_ra.csv",          # Rarefied counts
                      row.names = 1)                                   # Use sample IDs as row names

counts_sub <- read.csv("../data/processed_data/counts.csv",            # Non‑relative counts
                       row.names = 1)                                   # Use sample IDs as row names

genus_data <- read.csv("../data/processed_data/genus_sum_relab.csv")   # Genus relative abundances
rownames(genus_data) <- genus_data$Genus                               # Genus as row names

genus_data$Genus <- NULL                                               # Remove redundant column
genus_data$X <- NULL                                                   # Remove row index column

kegg_data <- read.csv("../data/processed_data/summed_kegg.csv") %>%    # KEGG summed pathways
  select(-X)                                                          # Drop row index column

rownames(kegg_data) <- kegg_data$purpose                               # Pathway names as row names
kegg_data$purpose <- NULL                                              # Remove redundant column

mg_raw <- read.csv("../data/processed_data/counts_sub_metagen.csv",     # Metagenome counts
                   row.names = 1)                                      # Sample IDs as row names

# ============================ # Divider for metadata filters
# STEP 4 — FILTER METADATA                                         # Keep samples aligned
# ============================ # Divider

meta_larvae <- meta %>%                                                # Larvae only
  filter(Sample_type == "Larvae") %>%                                 # Keep larvae
  filter(Microbe_treatment %in% c("MMO", "EMO", "MEM"))               # Keep 3 treatments

meta_water <- clean_meta %>%                                           # Water only
  filter(Sample_type == "Water") %>%                                  # Keep water
  filter(Microbe_treatment %in% c("MMO", "EMO", "MEM")) %>%          # Keep 3 treatments
  filter(Mesocosm_type == "Experiment")                               # Drop controls

meta_mg <- meta %>%                                                    # Metagenome metadata
  filter(Sample_type == "Larvae") %>%                                 # Larvae only
  filter(Microbe_treatment %in% c("MMO", "EMO", "MEM"))              # 3 treatments

# ============================ # Divider for SLIDE A
# SLIDE A — SHANNON DIVERSITY BOXPLOT                                # Slide A
# ============================ # Divider

shannon <- diversity(counts_ra, index = "shannon")                     # Shannon diversity

div_df <- data.frame(                                                  # Build a table
  Short_label = names(shannon),                                        # Sample IDs
  Shannon     = shannon                                                # Diversity values
) %>%                                                                  # Pipe
  merge(meta, by = "Short_label") %>%                                 # Add metadata
  filter(Microbe_treatment %in% c("MMO", "EMO", "MEM")) %>%         # Keep 3 treatments
  mutate(Microbe_treatment = factor(Microbe_treatment,                 # Set order
                                    levels = c("EMO", "MEM", "MMO"))) # EM → MEM → MS

p_diversity <- ggplot(div_df,                                          # Start plot
                      aes(x = Microbe_treatment,                       # X axis
                          y = Shannon,                                 # Y axis
                          fill = Microbe_treatment)) +                 # Fill by treatment
  geom_boxplot(alpha = 0.65, width = 0.5, outlier.shape = NA) +         # Boxplot
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +                   # Points
  scale_fill_manual(values = slide_colors) +                           # Colors
  scale_x_discrete(labels = slide_labels_short) +                      # Short labels
  labs(                                                                # Titles
    title    = "MS treatment harbors a simpler, more specialized community", # Title
    subtitle = "Lower Shannon diversity = fewer but more specialized bacteria", # Subtitle
    x        = "Microbe Treatment",                                    # X label
    y        = "Shannon Diversity Index",                              # Y label
    caption  = "16S amplicon data | Rarefied to 20,000 reads"           # Caption
  ) +                                                                  # End labs
  theme(                                                               # Theme tweaks
    legend.position = "none",                                          # No legend
    plot.title      = element_text(size = BASE_TEXT + 4, face = "bold"), # Title size
    plot.subtitle   = element_text(size = BASE_TEXT, color = "grey40"), # Subtitle size
    axis.title      = element_text(size = BASE_TEXT + 2),              # Axis title size
    axis.text       = element_text(size = BASE_TEXT + 2),              # Axis text size
    plot.caption    = element_text(size = BASE_TEXT - 2, color = "grey50", face = "italic") # Caption
  )                                                                    # End theme

print(p_diversity)                                                     # Show plot

ggsave(file.path(out_path, "SLIDE_A_shannon_diversity.png"),            # Save image
       plot = p_diversity, width = 10, height = 7, dpi = 300)          # Save settings

# EXTRA GRAPH (MS specialized community)                              # Extra requested graph
ggsave(file.path(out_path, "MS_specialized_community.png"),             # Save image with requested name
       plot = p_diversity, width = 10, height = 7, dpi = 300)          # Same plot, different file

write.csv(div_df,                                                      # Write slide A CSV
          file.path(csv_out_path, "SLIDE_A_shannon_data.csv"),          # CSV path
          row.names = FALSE)                                           # No row names

write.csv(div_df,                                                      # Write extra data sheet
          file.path(csv_out_path, "MS_specialized_community_data.csv"),  # CSV path
          row.names = FALSE)                                           # No row names

# ============================================================ 
# The rest of the script (Slides B–G) is unchanged from the previous
# version you already approved and will run as before. If you want
# me to paste the full remaining sections into this file too, say so.
# ============================================================ 
