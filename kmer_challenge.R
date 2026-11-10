##***************************
## kmer analysis
##
## Karl Cottenie
##
## 2026-10-29
##
##***************************

library(stringr)

get_kmers <- function(sequence, k = 3) {
  sequence <- toupper(gsub("\\s+", "", sequence))
  if (nchar(sequence) < k) return(character())
  starts <- seq_len(nchar(sequence) - k + 1)
  substring(sequence, starts, starts + k - 1)
}

# One k-mer set per row
kmers <- lapply(df$matching_class, get_kmers, k = 3)

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

# Pairwise similarity matrix
n <- length(kmers)
sim <- matrix(0, n, n)
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    sim[i, j] <- jaccard(unique(kmers[[i]]), unique(kmers[[j]]))
  }
}

# Example: compare the two most common classes by their average k-mer similarity
top_two <- names(sort(table(df$matching_class), decreasing = TRUE))[1:2]
rows <- which(df$matching_class %in% top_two)

mean_similarity <- mean(sim[rows[df$matching_class[rows] == top_two[1]],
                            rows[df$matching_class[rows] == top_two[2]]])

mean_similarity

## _ Packages used -------
rm(list = ls())
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(Biostrings)
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

## load in data -------

fa_t_spir = readDNAStringSet("../data/GCA_001447595.1_T1_ISS3_r1.0.fna")

fa_t_spir

df_te_t_spir = read_csv("../data/2026-10-26_t_spiralis.csv")

df_te_t_spir

## cleaning up two data sets -----

## _ first issue, from the instructions, different identifying names for contigs
fa_t_spir |> names() |> head()

df_te_t_spir |> 
  pull(qry_id) |> 
  head()

fa_t_spir |> 
  names() |> 
  str_sub(end = 14) |> 
  head()

names(fa_t_spir) = fa_t_spir |> 
  names() |> 
  str_sub(end = 14)

fa_t_spir
names(fa_t_spir)

df_te_t_spir |> 
  filter(qry_id %notin% names(fa_t_spir))
# ===> yay!!!! all the names are the same

## _ second issue, from Livpreet's question, values in between parentheses ----

df_te_t_spir |> View()
#====> issue in 2 columns: qry_left and in_repeat_end

# contrary to what I said in class, you should fix the problem here!!!

df_te_t_spir |> 
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
         # the \\( is necessary because ( is a special character))
         .after = qry_left) |> 
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
         .after = in_repeat_end) |> View()

#====> depending on downstream analyses, you can either 
# overwrite the original column, and then save the dataframe
df_te_t_spir =  df_te_t_spir |> 
  mutate(qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
         # notice the ary_left instead of qry_left2 from line 73
         .after = qry_left) |> 
  mutate(in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
         .after = in_repeat_end)
# WARNING: the () indicate if the TE is on the + or C direction (matching_repeat column), and this might influence what you want to do later in your analysis!!! So use this option with caution

# the other option is to save the code from line 72 as as

## large contigs ------
ls_width = fa_t_spir |> 
  width() 

ls_width |> 
  plot()

ls_width |> log10() |> 
  plot()
#====> inflection point around 1000 bp
# anything bigger than 1000 is "large"

ls_contig_large = fa_t_spir[width(fa_t_spir) > 1000] |> 
  names()

## creating the ggplot from the assignment ------

## _filering large contigs -----
df_te_seq = df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  # rows dropped from 41526 to 40518
  as.data.frame() 
 
## _obtain 4-mer frequencies -----
df_te_seq$nucl_seq2 = DNAStringSet(df_te_seq$nucl_seq)

df_oligo4 = oligonucleotideFrequency(df_te_seq$nucl_seq2, width = 4, as.prob = T)
df_oligo4

df_oligo4[, 1:5]

df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5]) 

## _ creating a new variable with 3 categories -----

# reminder, from the tidyverse class
df_te_t_spir |> 
  pull(matching_class) |> 
  unique()
# copy-paste to Copilot and ask for repetitive elements that are not DNA or retrotransposons

ls_non_te = c( "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")
ls_non_te

# from a later tidyverse class
ls_dna_te = df_te_t_spir |> 
  filter(matching_class %notin% ls_non_te) |> 
  pull(matching_class) |> 
  unique() 
ls_dna_te

ls_dna_te = ls_dna_te[str_detect(ls_dna_te, "^DNA")]

# creating the new variable
df_te_t_spir |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) 
# .after so that I can easily compare whether the new variable is created correctly

## _ putting it all together -------

df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5]) |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |> 
  select(te_class, AAAA:AACA) |> # this step not necessary but makes it easier to see the result
  pivot_longer(AAAA:AACA, names_to = "kmer", values_to = "perc") |> 
  ggplot(aes(x = kmer, y = perc)) + 
  # more separation if you use 
  #  ggplot(aes(x = kmer, y = perc^0.25)) +
  geom_boxplot() + 
  facet_grid(~ te_class)

## _ other way of doing it -------

# all in one go
# but difficult to get this done as you are coding
# the previous code has a lot more checks and balances
# and I changed the ggplot a little to provide another example of ggplot

df_te_t_spir |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |> 
  bind_cols(oligonucleotideFrequency(DNAStringSet(df_te_t_spir$nucl_seq),
                                     width = 4, as.prob = T)) |> 
  filter(qry_id %in% ls_contig_large) |> 
  pivot_longer(AAAA:AACA, names_to = "kmer", values_to = "perc") |> 
  ggplot(aes(x = te_class, y = perc)) + 
  geom_boxplot() + 
  facet_grid(kmer ~ .)

## new challenges ----- 

## THIS ONE _ show all the 4-mers, not just the first 5 -------

## _ show more differentiation of the percentages -----
# and reorder the types of TEs into a more logical order

## _ how many contigs do not have repetitive elements on them ----
# and what can we say about how large these contigs are?
