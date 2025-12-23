############################################################
## LGIM analysis : Gravity × Time × Donor (FULL R PIPELINE)
############################################################

## ---- libraries ----
library(tidyverse)
library(lme4)
library(lmerTest)
library(igraph)
library(permute)

## ---- load data ----
df <- read.csv("LGIM.csv", check.names = FALSE)

## ---- metadata / taxa columns ----
meta_cols <- c("SampleID","SampleType","Donor",
               "Gravity","Time","Replicate","TotalReads")
taxa_cols <- setdiff(colnames(df), meta_cols)

## ---- relative abundance ----
df_rel <- df %>%
  mutate(across(all_of(taxa_cols), ~ .x / TotalReads))

## ---- time numeric ----
df_rel$TimeNum <- as.numeric(sub("h","", df_rel$Time))

############################################################
## 1. MIXED EFFECT MODEL (Donor as random effect)
############################################################

mixed_results <- df_rel %>%
  filter(Gravity %in% c("0g","1g")) %>%
  pivot_longer(all_of(taxa_cols),
               names_to="Taxon",
               values_to="Abundance") %>%
  group_by(Taxon) %>%
  summarise(
    model = list(
      lmer(Abundance ~ Gravity * TimeNum + (1|Donor),
           data = cur_data())
    ),
    p_gravity = summary(model[[1]])$coefficients["Gravity0g","Pr(>|t|)"],
    p_interaction = summary(model[[1]])$coefficients["Gravity0g:TimeNum","Pr(>|t|)"],
    .groups="drop"
  )

############################################################
## 2. PERMUTATION TEST (0g vs 1g)
############################################################

perm_test <- function(df, taxon, nperm=999) {
  sub <- df %>%
    filter(Gravity %in% c("0g","1g")) %>%
    select(Gravity, all_of(taxon))

  obs <- mean(sub[[taxon]][sub$Gravity=="0g"]) -
         mean(sub[[taxon]][sub$Gravity=="1g"])

  perm <- replicate(nperm, {
    g <- sample(sub$Gravity)
    mean(sub[[taxon]][g=="0g"]) -
      mean(sub[[taxon]][g=="1g"])
  })

  mean(abs(perm) >= abs(obs))
}

perm_results <- tibble(
  Taxon = taxa_cols,
  p_perm = map_dbl(taxa_cols, ~ perm_test(df_rel, .x))
)

############################################################
## 3. NETWORK ANALYSIS (igraph)
############################################################

## ---- build correlation network (0g only) ----
net_df <- df_rel %>%
  filter(Gravity=="0g") %>%
  select(all_of(taxa_cols))

cor_mat <- cor(net_df, method="spearman")

## threshold
edges <- which(abs(cor_mat) > 0.6 & upper.tri(cor_mat), arr.ind=TRUE)

edge_list <- data.frame(
  from = rownames(cor_mat)[edges[,1]],
  to   = colnames(cor_mat)[edges[,2]],
  weight = cor_mat[edges]
)

g <- graph_from_data_frame(edge_list, directed=FALSE)

network_metrics <- tibble(
  Taxon = V(g)$name,
  Degree = degree(g),
  Betweenness = betweenness(g),
  Eigenvector = eigen_centrality(g)$vector
)

############################################################
## 4. TIME × GRAVITY INTERACTION (fixed-effect summary)
############################################################

interaction_results <- df_rel %>%
  pivot_longer(all_of(taxa_cols),
               names_to="Taxon",
               values_to="Abundance") %>%
  group_by(Taxon) %>%
  summarise(
    interaction_p =
      summary(lm(Abundance ~ Gravity * TimeNum,
                 data = cur_data()))$coefficients["Gravity0g:TimeNum","Pr(>|t|)"],
    .groups="drop"
  )

############################################################
## 5. INTEGRATION OF ALL RESULTS
############################################################

final_results <- mixed_results %>%
  left_join(perm_results, by="Taxon") %>%
  left_join(network_metrics, by="Taxon") %>%
  left_join(interaction_results, by="Taxon")

## ---- top 10 taxa (multi-criteria) ----
top10 <- final_results %>%
  arrange(p_gravity, p_interaction, p_perm, desc(Eigenvector)) %>%
  slice_head(n=10)

############################################################
## OUTPUT
############################################################

print(top10)

write.csv(top10,
          "Top10_Gravity_Time_Donor_Impacted_Taxa.csv",
          row.names=FALSE)

############################################################
## END
############################################################
