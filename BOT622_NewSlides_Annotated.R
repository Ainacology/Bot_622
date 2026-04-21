# ============================================================ 
# Title line for the script
# BOT 622 — NEW SLIDES (BEGINNER ANNOTATED, LINE‑BY‑LINE)      
# Explains what this script is
# This script reproduces NEW SLIDES A–G with beginner notes     
# High‑level summary
# It also writes one CSV per slide with the plotted data        
# CSV promise
# Style matches the provided example scripts                    
# ============================================================ 

# ============================ # Divider for setup
# STEP 0 — SET PATHS (EDIT IF NEEDED)                            # Tell user what to edit
# ============================ # Divider

base <- "C:/Users/ntang/OneDrive/Desktop/Repositories/Bot_622"                 # Root project folder
out_path <- file.path(base, "SKY")                                             # Folder for slide images
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

ggsave(file.path(out_path, "SLIDE_A_shannon_diversity.png"),            # Save image (fixed ASCII A)
       plot = p_diversity, width = 10, height = 7, dpi = 300)          # Save settings

write.csv(div_df,                                                      # Write slide A CSV
          file.path(csv_out_path, "SLIDE_A_shannon_data.csv"),          # CSV path
          row.names = FALSE)                                           # No row names

# ============================ # Divider for SLIDE B
# SLIDE B — TOP 3 GENERA ENRICHED IN MS                               # Slide B
# ============================ # Divider

genus_filt <- genus_data[, colnames(genus_data) %in% meta_larvae$Short_label] # Keep larvae samples

genus_long <- genus_filt %>%                                           # Start long format
  rownames_to_column("Genus") %>%                                     # Row names to column
  pivot_longer(-Genus,                                                 # Make long
               names_to = "Short_label",                               # Sample ID column
               values_to = "Rel_abund") %>%                            # Abundance column
  left_join(meta_larvae %>% select(Short_label, Microbe_treatment),     # Add treatment
            by = "Short_label") %>%                                    # Join key
  group_by(Genus, Microbe_treatment) %>%                               # Group
  summarise(mean_abund = mean(Rel_abund), .groups = "drop")            # Mean

# Compare MS vs EM to rank genera                                      # Explanation

genus_compare <- genus_long %>%                                        # Start compare table
  pivot_wider(names_from = Microbe_treatment,                          # Columns per treatment
              values_from = mean_abund,                                # Values = mean
              values_fill = 0) %>%                                     # Fill missing
  mutate(MS_vs_EM = MMO - EMO) %>%                                     # Difference
  arrange(desc(MS_vs_EM))                                              # Sort high → low

# Target genera (consistent with your slide plan)                      # Explanation

target_genera <- c("Yersinia", "Aquitalea", "Clostridium sensu stricto 3") # 3 genera

three_genera <- genus_long %>%                                         # Filter to targets
  filter(Genus %in% target_genera) %>%                                 # Keep 3 genera
  mutate(                                                              # Set order
    Genus = factor(Genus, levels = rev(target_genera)),                # Reverse order
    Microbe_treatment = factor(Microbe_treatment, levels = c("MMO","EMO","MEM")) # MS/EM/MEM
  )                                                                    # End mutate

p_three_genera <- ggplot(three_genera,                                 # Start plot
                         aes(x = Genus,                                # X axis
                             y = mean_abund,                           # Y axis
                             fill = Microbe_treatment)) +              # Fill by treatment
  geom_col(position = "dodge", width = 0.65) +                         # Side‑by‑side bars
  scale_fill_manual(values = slide_colors, labels = slide_labels_short) + # Colors
  coord_flip() +                                                       # Horizontal
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) + # Percent y axis
  labs(                                                                # Titles
    title    = "Bacteria enriched in MS: fitness candidates",          # Title
    subtitle = "These genera are substantially higher in MS (red) than EM (green)", # Subtitle
    x        = "Genus",                                                # X label
    y        = "Mean relative abundance",                              # Y label
    fill     = "Treatment"                                             # Legend title
  ) +                                                                  # End labs
  theme(                                                               # Theme
    plot.title    = element_text(size = BASE_TEXT + 4, face = "bold"), # Title
    plot.subtitle = element_text(size = BASE_TEXT - 1, color = "grey40"), # Subtitle
    axis.title    = element_text(size = BASE_TEXT + 2),                # Axis titles
    axis.text     = element_text(size = BASE_TEXT + 2),                # Axis text
    axis.text.y   = element_text(face = "italic"),                     # Italic genus
    legend.title  = element_text(size = BASE_TEXT),                    # Legend title
    legend.text   = element_text(size = BASE_TEXT)                     # Legend text
  )                                                                    # End theme

print(p_three_genera)                                                  # Show plot

ggsave(file.path(out_path, "SLIDE_B_top3_genera_MS.png"),              # Save image
       plot = p_three_genera, width = 12, height = 7, dpi = 300)       # Save settings

write.csv(three_genera,                                                # Write slide B CSV
          file.path(csv_out_path, "SLIDE_B_top3_genera_data.csv"),      # CSV path
          row.names = FALSE)                                           # No row names

# ============================ # Divider for SLIDE C
# SLIDE C — KEGG PATHWAYS ENRICHED IN MS                              # Slide C
# ============================ # Divider

kegg_filt <- kegg_data[, colnames(kegg_data) %in% meta_mg$Short_label]  # Filter to metagenome samples

kegg_long <- kegg_filt %>%                                              # Long format KEGG
  rownames_to_column("Pathway") %>%                                    # Row names to column
  pivot_longer(-Pathway,                                                # Long format
               names_to = "Short_label",                               # Sample IDs
               values_to = "Abundance") %>%                            # Abundance
  left_join(meta_mg %>% select(Short_label, Microbe_treatment),          # Add treatment
            by = "Short_label") %>%                                     # Join key
  group_by(Pathway, Microbe_treatment) %>%                              # Group
  summarise(mean_abund = mean(Abundance), .groups = "drop")             # Mean

kegg_long <- kegg_long %>%                                              # Shorten labels
  mutate(Pathway_short = case_when(                                     # Name mapping
    grepl("signaling and cellular", Pathway) ~ "Protein signaling",     # Short name
    grepl("Infectious diseases",    Pathway) ~ "Host-bacteria interaction", # Short name
    grepl("Cell motility",          Pathway) ~ "Cell motility",         # Keep
    grepl("Protein families: metabolism", Pathway) ~ "Protein metabolism", # Short name
    grepl("genetic information",    Pathway) ~ "Genetic info processing", # Short name
    grepl("Signal transduction",    Pathway) ~ "Signal transduction",   # Short name
    grepl("Antimicrobial",          Pathway) ~ "Antimicrobial resistance", # Short name
    grepl("amino acids",            Pathway) ~ "Amino acid metabolism", # Short name
    grepl("cofactors and vitamins", Pathway) ~ "Cofactors & vitamins", # Short name
    grepl("Replication",            Pathway) ~ "Replication & repair",  # Short name
    grepl("Translation",            Pathway) ~ "Translation",          # Short name
    grepl("Lipid",                  Pathway) ~ "Lipid metabolism",      # Short name
    grepl("Carbohydrate",           Pathway) ~ "Carbohydrate metabolism", # Short name
    TRUE ~ Pathway                                                  # Default
  ))                                                                     # End mutate

target_pathways <- c(                                                   # 4 focal pathways
  "Host-bacteria interaction",                                          # Pathway 1
  "Cell motility",                                                     # Pathway 2
  "Antimicrobial resistance",                                          # Pathway 3
  "Cofactors & vitamins"                                               # Pathway 4
)                                                                       # End list

four_pathways <- kegg_long %>%                                          # Filter to 4
  filter(Pathway_short %in% target_pathways) %>%                        # Keep 4
  mutate(                                                               # Set order
    Pathway_short = factor(Pathway_short, levels = rev(target_pathways)), # Order
    Microbe_treatment = factor(Microbe_treatment, levels = c("MMO","EMO","MEM")) # Order
  )                                                                     # End mutate

p_kegg_four <- ggplot(four_pathways,                                    # Start plot
                      aes(x = Pathway_short, y = mean_abund, fill = Microbe_treatment)) + # Aes
  geom_col(position = "dodge", width = 0.65) +                         # Bars
  scale_fill_manual(values = slide_colors, labels = slide_labels_short) + # Colors
  coord_flip() +                                                        # Horizontal
  scale_y_continuous(labels = scales::comma) +                          # Commas
  labs(                                                                 # Titles
    title    = "MS bacteria carry host-support functional genes",       # Title
    subtitle = "KEGG pathways enriched in MS (red) — absent or lower in EM (green)", # Subtitle
    x        = "KEGG Pathway",                                          # X label
    y        = "Mean gene abundance",                                   # Y label
    fill     = "Treatment"                                              # Legend
  ) +                                                                   # End labs
  theme(                                                                # Theme
    plot.title    = element_text(size = BASE_TEXT + 4, face = "bold"),  # Title
    plot.subtitle = element_text(size = BASE_TEXT - 1, color = "grey40"), # Subtitle
    axis.title    = element_text(size = BASE_TEXT + 2),                 # Axis titles
    axis.text     = element_text(size = BASE_TEXT + 2),                 # Axis text
    legend.title  = element_text(size = BASE_TEXT),                     # Legend title
    legend.text   = element_text(size = BASE_TEXT)                      # Legend text
  )                                                                     # End theme

print(p_kegg_four)                                                      # Show plot

ggsave(file.path(out_path, "SLIDE_C_kegg_MS_enriched.png"),             # Save image
       plot = p_kegg_four, width = 12, height = 7, dpi = 300)           # Save settings

write.csv(four_pathways,                                                # Write slide C CSV
          file.path(csv_out_path, "SLIDE_C_kegg_four_pathways.csv"),    # CSV path
          row.names = FALSE)                                            # No row names

# ============================ # Divider for SLIDE D
# SLIDE D — SPAGHETTI SLOPE GRAPH (WATER → LARVAE)                      # Slide D
# ============================ # Divider

genus_water <- genus_data[, colnames(genus_data) %in% meta_water$Short_label] # Water samples

genus_larvae <- genus_data[, colnames(genus_data) %in% meta_larvae$Short_label] # Larvae samples

water_long <- genus_water %>%                                           # Water long format
  rownames_to_column("Genus") %>%                                      # Genus column
  pivot_longer(-Genus, names_to = "Short_label", values_to = "Rel_abund") %>% # Long
  left_join(meta_water %>% select(Short_label, Microbe_treatment), by = "Short_label") %>% # Join
  filter(!is.na(Microbe_treatment)) %>%                                 # Drop missing
  group_by(Genus, Microbe_treatment) %>%                                # Group
  summarise(mean_abund = mean(Rel_abund), .groups = "drop") %>%         # Mean
  mutate(Sample_type = "Water")                                        # Tag

larvae_long <- genus_larvae %>%                                         # Larvae long format
  rownames_to_column("Genus") %>%                                      # Genus column
  pivot_longer(-Genus, names_to = "Short_label", values_to = "Rel_abund") %>% # Long
  left_join(meta_larvae %>% select(Short_label, Microbe_treatment), by = "Short_label") %>% # Join
  group_by(Genus, Microbe_treatment) %>%                                # Group
  summarise(mean_abund = mean(Rel_abund), .groups = "drop") %>%         # Mean
  mutate(Sample_type = "Larvae")                                       # Tag

combined <- bind_rows(water_long, larvae_long) %>%                      # Combine
  mutate(Sample_type = factor(Sample_type, levels = c("Water", "Larvae"))) # Order

top_genera_slope <- combined %>%                                        # Pick top genera
  group_by(Genus) %>%                                                   # Group by genus
  summarise(total = sum(mean_abund)) %>%                                # Total abundance
  arrange(desc(total)) %>%                                              # Sort
  head(8) %>%                                                           # Top 8
  pull(Genus)                                                           # Vector

plot_data <- combined %>%                                               # Filter to top 8
  filter(Genus %in% top_genera_slope)                                   # Keep

p_spaghetti <- ggplot(plot_data,                                        # Start plot
                      aes(x = Sample_type, y = mean_abund, group = Genus, color = Genus)) + # Aes
  geom_line(linewidth = 1.2, alpha = 0.85) +                            # Lines
  geom_point(size = 4) +                                                # Points
  geom_text(                                                            # Label on larvae side
    data = plot_data %>% filter(Sample_type == "Larvae"),               # Larvae only
    aes(label = Genus),                                                 # Label text
    hjust = -0.12,                                                      # Nudge right
    size = BASE_TEXT * 0.28,                                            # Size
    fontface = "italic"                                                # Italic
  ) +                                                                   # End labels
  facet_wrap(~ Microbe_treatment,                                       # One panel per treatment
             labeller = labeller(Microbe_treatment = c(                # Friendly labels
               "MMO" = "MS (Mosquito Symbionts)",                      # MS label
               "EMO" = "EM (Environmental Microbes)",                  # EM label
               "MEM" = "MS+EM (Combo)"                                 # MEM label
             ))) +                                                     # End labeller
  scale_y_log10(labels = scales::percent_format(accuracy = 0.01),       # Log scale
                breaks = c(0.001, 0.01, 0.1, 1) / 10) +                 # Breaks
  labs(                                                                 # Titles
    title    = "Larval gut actively selects which water bacteria establish", # Title
    subtitle = "Lines going UP = enriched by larva | Lines going DOWN = filtered out", # Subtitle
    x        = "Sample type",                                           # X label
    y        = "Mean relative abundance (log scale)",                   # Y label
    color    = "Genus"                                                  # Legend
  ) +                                                                   # End labs
  scale_x_discrete(expand = expansion(add = c(0.3, 1.5))) +             # Room for labels
  theme(                                                                # Theme
    legend.position = "none",                                          # Hide legend
    strip.text      = element_text(size = BASE_TEXT + 1, face = "bold"),# Facet labels
    axis.title      = element_text(size = BASE_TEXT + 2),               # Axis titles
    axis.text       = element_text(size = BASE_TEXT + 1),               # Axis text
    plot.title      = element_text(size = BASE_TEXT + 4, face = "bold"),# Title
    plot.subtitle   = element_text(size = BASE_TEXT, color = "grey40"), # Subtitle
    panel.spacing   = unit(2, "cm")                                    # Space
  )                                                                     # End theme

print(p_spaghetti)                                                      # Show plot

ggsave(file.path(out_path, "SLIDE_D_spaghetti_water_larvae.png"),       # Save image
       plot = p_spaghetti, width = 18, height = 9, dpi = 300)           # Save settings

write.csv(plot_data,                                                    # Write slide D CSV
          file.path(csv_out_path, "SLIDE_D_spaghetti_data.csv"),        # CSV path
          row.names = FALSE)                                            # No row names

# ============================ # Divider for SLIDES E + F
# SLIDES E + F — METAGENOME ORDINATIONS                              # Slide E/F
# ============================ # Divider

mg_counts <- as.data.frame(t(mg_raw))                                   # Transpose to samples x genes
mg_counts <- mg_counts[rownames(mg_counts) %in% meta_mg$Short_label, ]   # Keep samples with metadata
meta_mg_ord <- meta_mg[match(rownames(mg_counts), meta_mg$Short_label), ]# Align metadata order

# ---- Bray-Curtis (Slide F) ----                                     # Label
mg_dist_bray <- vegdist(mg_counts, method = "bray")                     # Bray distance
set.seed(123)                                                          # Reproducible NMDS
capture.output(mg_ord_bray <- metaMDS(mg_dist_bray, k = 2))             # NMDS

bray_pts <- as.data.frame(mg_ord_bray$points) %>%                       # Extract points
  rownames_to_column("Short_label") %>%                                # Sample IDs
  merge(meta_mg_ord, by = "Short_label") %>%                            # Add metadata
  mutate(Microbe_treatment = as.factor(Microbe_treatment),              # Factor
         Stress = mg_ord_bray$stress)                                   # Add stress

p_bray <- ggplot(bray_pts, aes(MDS1, MDS2)) +                            # Start plot
  geom_point(aes(color = Microbe_treatment), size = 5) +                # Points
  scale_color_manual(values = slide_colors, labels = slide_labels_short)+ # Colors
  labs(                                                                 # Titles
    title    = "Metagenome ordination — Bray-Curtis distance",           # Title
    subtitle = "Functional gene profiles | Standard approach",           # Subtitle
    caption  = paste("Stress =", round(mg_ord_bray$stress, 4)),           # Stress
    x        = "MDS1",                                                   # X label
    y        = "MDS2",                                                   # Y label
    color    = "Treatment"                                              # Legend
  ) +                                                                   # End labs
  theme(                                                                # Theme
    plot.title    = element_text(size = BASE_TEXT + 4, face = "bold"),  # Title
    plot.subtitle = element_text(size = BASE_TEXT, color = "grey40"),   # Subtitle
    axis.title    = element_text(size = BASE_TEXT + 2),                 # Axis titles
    axis.text     = element_text(size = BASE_TEXT + 1),                 # Axis text
    legend.title  = element_text(size = BASE_TEXT),                     # Legend title
    legend.text   = element_text(size = BASE_TEXT),                     # Legend text
    plot.caption  = element_text(size = BASE_TEXT, color = "grey50", face = "italic") # Caption
  )                                                                     # End theme

print(p_bray)                                                           # Show plot

ggsave(file.path(out_path, "SLIDE_F_metagenome_braycurtis.png"),        # Save image
       plot = p_bray, width = 10, height = 8, dpi = 300)                # Save settings

write.csv(bray_pts,                                                     # Write slide F CSV
          file.path(csv_out_path, "SLIDE_F_bray_points.csv"),           # CSV path
          row.names = FALSE)                                            # No row names

# ---- Aitchison (Slide E) ----                                        # Label
mg_adj <- mg_counts + 0.5                                               # Add pseudocount
mg_clr <- decostand(mg_adj, method = "clr")                             # CLR transform
mg_dist_ait <- dist(mg_clr, method = "euclidean")                       # Aitchison distance
set.seed(123)                                                          # Reproducible NMDS
capture.output(mg_ord_ait <- metaMDS(mg_dist_ait, k = 2))               # NMDS

ait_pts <- as.data.frame(mg_ord_ait$points) %>%                         # Extract points
  rownames_to_column("Short_label") %>%                                # Sample IDs
  merge(meta_mg_ord, by = "Short_label") %>%                            # Add metadata
  mutate(Microbe_treatment = as.factor(Microbe_treatment),              # Factor
         Stress = mg_ord_ait$stress)                                    # Add stress

p_aitchison <- ggplot(ait_pts, aes(MDS1, MDS2)) +                        # Start plot
  geom_point(aes(color = Microbe_treatment), size = 5) +                # Points
  scale_color_manual(values = slide_colors, labels = slide_labels_short)+ # Colors
  labs(                                                                 # Titles
    title    = "Metagenome ordination — Aitchison distance",             # Title
    subtitle = "CLR-transformed profiles | Compositional data approach",# Subtitle
    caption  = paste("Stress =", round(mg_ord_ait$stress, 4)),           # Stress
    x        = "MDS1",                                                   # X label
    y        = "MDS2",                                                   # Y label
    color    = "Treatment"                                              # Legend
  ) +                                                                   # End labs
  theme(                                                                # Theme
    plot.title    = element_text(size = BASE_TEXT + 4, face = "bold"),  # Title
    plot.subtitle = element_text(size = BASE_TEXT, color = "grey40"),   # Subtitle
    axis.title    = element_text(size = BASE_TEXT + 2),                 # Axis titles
    axis.text     = element_text(size = BASE_TEXT + 1),                 # Axis text
    legend.title  = element_text(size = BASE_TEXT),                     # Legend title
    legend.text   = element_text(size = BASE_TEXT),                     # Legend text
    plot.caption  = element_text(size = BASE_TEXT, color = "grey50", face = "italic") # Caption
  )                                                                     # End theme

print(p_aitchison)                                                      # Show plot

ggsave(file.path(out_path, "SLIDE_E_metagenome_aitchison.png"),         # Save image
       plot = p_aitchison, width = 10, height = 8, dpi = 300)           # Save settings

write.csv(ait_pts,                                                      # Write slide E CSV
          file.path(csv_out_path, "SLIDE_E_aitchison_points.csv"),      # CSV path
          row.names = FALSE)                                            # No row names

# ============================ # Divider for SLIDE G
# SLIDE G — KEGG LOLLIPOP (MS vs EM)                                  # Slide G
# ============================ # Divider

kegg_filt2 <- kegg_data[, colnames(kegg_data) %in% meta_mg$Short_label]  # Filter KEGG to samples

mmo_k_means <- rowMeans(kegg_filt2[, colnames(kegg_filt2) %in% meta_mg$Short_label[meta_mg$Microbe_treatment == "MMO"]]) # MS mean
emo_k_means <- rowMeans(kegg_filt2[, colnames(kegg_filt2) %in% meta_mg$Short_label[meta_mg$Microbe_treatment == "EMO"]]) # EM mean
mem_k_means <- rowMeans(kegg_filt2[, colnames(kegg_filt2) %in% meta_mg$Short_label[meta_mg$Microbe_treatment == "MEM"]]) # MEM mean

kegg_compare <- data.frame(                                            # Build comparison
  Pathway   = rownames(kegg_data),                                     # Pathway names
  MS_mean   = mmo_k_means,                                             # MS mean
  EM_mean   = emo_k_means,                                             # EM mean
  MSem_mean = mem_k_means                                              # MEM mean
) %>%                                                                  # Pipe
  mutate(MS_vs_EM = MS_mean - EM_mean) %>%                             # Difference
  arrange(desc(MS_vs_EM)) %>%                                          # Sort
  mutate(Pathway_short = case_when(                                    # Shorten labels
    grepl("signaling and cellular", Pathway) ~ "Protein signaling",     # Short
    grepl("Infectious diseases",    Pathway) ~ "Host-bacteria interaction", # Short
    grepl("Unclassified: signaling", Pathway) ~ "Unclassified signaling", # Short
    grepl("Cell motility",          Pathway) ~ "Cell motility",         # Short
    grepl("Protein families: metabolism", Pathway) ~ "Protein metabolism", # Short
    grepl("genetic information",    Pathway) ~ "Genetic info processing", # Short
    grepl("Signal transduction",    Pathway) ~ "Signal transduction",   # Short
    grepl("Antimicrobial",          Pathway) ~ "Antimicrobial resistance", # Short
    grepl("amino acids",            Pathway) ~ "Amino acid metabolism", # Short
    grepl("cofactors and vitamins", Pathway) ~ "Cofactors & vitamins", # Short
    grepl("Replication",            Pathway) ~ "Replication & repair",  # Short
    grepl("Translation",            Pathway) ~ "Translation",          # Short
    TRUE ~ Pathway                                                     # Default
  ))                                                                   # End mutate

kegg_top12 <- kegg_compare %>%                                         # Top 12
  head(12) %>%                                                         # Keep 12
  group_by(Pathway_short) %>%                                          # Detect duplicates
  mutate(                                                              # Fix duplicate labels
    n = row_number(),                                                  # Count duplicates
    Pathway_short = ifelse(n > 1, paste0(Pathway_short, " (", n, ")"), Pathway_short) # Add suffix
  ) %>%                                                                # End mutate
  ungroup() %>%                                                        # Ungroup
  select(-n) %>%                                                       # Drop helper
  mutate(Pathway_short = factor(Pathway_short, levels = rev(unique(Pathway_short)))) # Order

p_lollipop <- ggplot(kegg_top12,                                       # Start plot
                     aes(x = Pathway_short, y = MS_vs_EM)) +           # Aes
  geom_segment(aes(xend = Pathway_short, y = 0, yend = MS_vs_EM),        # Lollipop stem
               color = "grey60", linewidth = 1.2) +                   # Stem style
  geom_point(aes(size = MS_mean, color = MS_vs_EM)) +                  # Lollipop dot
  scale_color_gradient(low = "orange", high = "red", name = "MS - EM\ndifference") + # Color
  scale_size_continuous(range = c(4, 12), name = "MS abundance", labels = scales::comma) + # Size
  coord_flip() +                                                       # Horizontal
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) + # Zero line
  labs(                                                                # Titles
    title    = "Which functional genes most distinguish MS from EM?",   # Title
    subtitle = "Bars right = higher in MS | Dot size = MS abundance",   # Subtitle
    x        = "KEGG Pathway",                                         # X label
    y        = "MS - EM difference in abundance"                       # Y label
  ) +                                                                  # End labs
  theme(                                                               # Theme
    plot.title    = element_text(size = BASE_TEXT + 4, face = "bold"),  # Title
    plot.subtitle = element_text(size = BASE_TEXT - 1, color = "grey40"), # Subtitle
    axis.title    = element_text(size = BASE_TEXT + 2),                 # Axis titles
    axis.text.y   = element_text(size = BASE_TEXT + 1),                 # Y text
    axis.text.x   = element_text(size = BASE_TEXT),                     # X text
    legend.title  = element_text(size = BASE_TEXT - 1),                 # Legend title
    legend.text   = element_text(size = BASE_TEXT - 1)                  # Legend text
  )                                                                     # End theme

print(p_lollipop)                                                      # Show plot

ggsave(file.path(out_path, "SLIDE_G_kegg_lollipop_MS_vs_EM.png"),       # Save image
       plot = p_lollipop, width = 14, height = 9, dpi = 300)            # Save settings

write.csv(kegg_top12,                                                   # Write slide G CSV
          file.path(csv_out_path, "SLIDE_G_kegg_lollipop_data.csv"),    # CSV path
          row.names = FALSE)                                            # No row names

# ============================ # Divider for done message
# DONE — QUICK CONFIRMATION                                           # End note
# ============================ # Divider

cat("\nAll slide images saved to:", out_path, "\n")                     # Confirm images
cat("All slide CSVs saved to:", csv_out_path, "\n")                     # Confirm CSVs
