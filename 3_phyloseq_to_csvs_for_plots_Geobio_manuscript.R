

### Manipulation of NTU_combined data in order to calculate relative abundance and format data for plotting ###

library(forcats)
library(tidyr)
library(dplyr)
library(tibble)
library(here)


### Load our exported Phyloseq BIOM format file to be handled using the tidyverse

phyloseq_biom <- read.csv("phyloseq_biom.csv", header = TRUE, sep = ",")

names(phyloseq_biom) <- gsub("\\.", "-", names(phyloseq_biom)) # replace periods that are inserted in column names. If you use check_names = FALSE with read.csv it will throw an error in the next step with pivot_longer. This is a silly workaround to preserve the hyphen in the sample name.

phyloseq_biom_long <- phyloseq_biom %>% 
  ## !!! check classes, filter below only necessary 
  # because classes are missing!
  #filter(!is.na(class)) %>% 
  pivot_longer(cols = c("GR19-MIT5", "IS19-13", "IS19-14"),  
               names_to = "Site", 
               values_to = "Abs_vals")


# calculations should be per domain!
# Check ASVs for each domain
table(phyloseq_biom_long$domain, useNA = "always")

phyloseq_biom_long <- phyloseq_biom_long %>% 
  mutate(
    domain = ifelse(domain == "Eukaryota", domain, "Bacteria")
  ) # careful as this will change any Archaea to Bacteria

# Check domains again
table(phyloseq_biom_long$domain, useNA = "always")

# Convert sequencing values to relative abundance and add as a new column (take note of the group and ungroup functions!!)
phyloseq_biom_long <- phyloseq_biom_long %>% 
  group_by(Site, domain) %>% 
  mutate(
    Rel_vals = 100*Abs_vals/sum(Abs_vals, na.rm = TRUE)
  ) %>% 
  ungroup()

# Convert from long data in order to export to CSV for working in Excel.
BONCAT_rel_vals <- phyloseq_biom_long %>%
  select(-Abs_vals) %>% 
  pivot_wider(
    #id_cols = c("GR19-MIT5", "IS19-13", "IS19-14"),
    #names_prefix = c("Abs_", "Rel_"),
    names_from = "Site",
    values_from = c("Rel_vals")
      )

write.csv(BONCAT_rel_vals, file = here("output", "TotalRNA_ALL_rel_abund.csv"))


# Pivot table but with Absolute values as well
BONCAT_rel_abs_vals <- phyloseq_biom_long %>%
  #select(-Abs_vals) %>% 
  pivot_wider(
    #id_cols = c("GR19-MIT5", "IS19-13", "IS19-14"),
    #names_prefix = c("Abs_", "Rel_"),
    names_from = "Site",
    values_from = c("Rel_vals", "Abs_vals")
  )

write.csv(BONCAT_rel_abs_vals, file = here("output", "TotalRNA_ALL_rel_abs_abund.csv"))


table(BONCAT_rel_vals$domain)


BONCAT_rel_vals_16S <- BONCAT_rel_vals %>% 
  filter(domain == "Bacteria")



BONCAT_rel_vals_18S <- BONCAT_rel_vals %>% 
  filter(domain == "Eukaryota")


# Create class groups 
phyloseq_biom_class <- phyloseq_biom_long %>% 
  group_by(domain, Site, class) %>% 
  mutate(class_sum_val = sum(Rel_vals)) %>% 
  ungroup()

# Make a list of the the top classes
top_classes <- phyloseq_biom_class %>% 
  select(domain, Site, class, class_sum_val) %>% 
  unique() %>% 
  group_by(domain, Site) %>% 
  slice_max(class_sum_val, n=4) %>% 
  pull(class) %>% unique()


# now create new class, that aggregrates small classes into an other group:
phyloseq_biom_class <- phyloseq_biom_class %>% 
  mutate(
    class_aggr = fct_other(class, keep = top_classes)
  )

table(phyloseq_biom_class$class_aggr)

# Top genera
top_genus <- phyloseq_biom_class %>% 
  group_by(domain, Site, class_aggr) %>% 
  slice_max(Rel_vals, n = 3) %>% 
  filter(Rel_vals > 0.04) %>% 
  pull(genus) %>% 
  unique()

phyloseq_biom_class <- phyloseq_biom_class %>% 
  ungroup() %>% 
  mutate(
    genus_aggr = fct_other(genus, keep = top_genus),
    class_genus = paste(class_aggr, genus_aggr, sep = " - "),
    class_genus = gsub("^Other - .*$", "Other", class_genus)
  )

table(phyloseq_biom_class$class_genus)


phyloseq_biom_class_agg <- phyloseq_biom_class %>% 
  group_by(domain, Site, class_genus) %>% 
  summarize(Rel_vals = sum(Rel_vals))

# Just to check that values equal 100%
phyloseq_biom_class_agg %>% 
  group_by(domain, Site) %>% 
  summarize(sum(Rel_vals))


phyloseq_biom_class_agg_wide <- phyloseq_biom_class_agg %>%
  pivot_wider(
    names_from = "Site",
    values_from = c("Rel_vals")
  )


phyloseq_biom_class_agg_wide_16S <- phyloseq_biom_class_agg_wide %>% 
  filter(domain == "Bacteria")

write.csv(phyloseq_biom_class_agg_wide_16S, file = here("output", "TotalRNA_class_genus_16S.csv"))



phyloseq_biom_class_agg_wide_18S <- phyloseq_biom_class_agg_wide %>% 
  filter(domain == "Eukaryota")

write.csv(phyloseq_biom_class_agg_wide_18S, file = here("output", "TotalRNA_class_genus_18S.csv"))



