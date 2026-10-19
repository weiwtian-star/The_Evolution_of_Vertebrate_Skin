#####speckle#####
library(speckle)
props <- getTransformedProps(sce$celltype, sce$orig.ident, transform="logit")

A<-data.frame('orig.ident'=colnames(props$Counts))
B<-read_xlsx('d:/AAA_Skin_revise/Q2/A.xlsx')

C<-left_join(A,B)

site<-as.factor(C$site)
species<-as.factor(C$species)

design <- model.matrix(~ 0 + species + site)

result_species <- propeller.anova(
  prop.list = props,
  design    = design ,
  coef      = c(1:31),
  robust    = TRUE,
  trend     = FALSE,
  sort      = TRUE
)

library(tidyverse)

csv1<-result_species[,1:31]

df <- csv1 %>%
  tibble::rownames_to_column("CellType")

df_long <- df %>%
  pivot_longer(
    cols = -CellType,
    names_to = "Species",
    values_to = "Proportion"
  )

cell.cols <- c(
  "Differentiated_KC" = "#89C860",
  "Undifferentiated_KC"   = "#3CBBA6",
  "Fibroblast"       = "#F37C27",
  "Endothelial_cells"  = "#F8F197",
  "Melanocytes"   = "#7D1418",
  "Immune_cells"= "#C16CAD",
  "Sebaceous_gland_channel"   = "#D670BA",
  "Sebaceous_gland_secretion"   = "#CE0B95",
  "SMC_Pericytes"   = "#FCB423",
  "Merkel_cell"   = "#d5a7b1",
  "Mucous_cell"   = "#c2bf53",
  "Proliferating_cells"   = "#dcad7f",
  "Ionocytes"   = "#454a80",
  "Hair_follicle_stem_cells"   = "#7FC6E6",
  "Hair_follicle_matrix_cells"   = "#5DA3D4",
  "Goblet_cell"   = "#5f89a7",
  "Schwann"   = "#687286"
)

df_long$CellType <- factor(
  df_long$CellType,
  levels = rev(c(  "Differentiated_KC" ,
                   "Undifferentiated_KC"  ,
                   "Fibroblast" ,
                   "Endothelial_cells" ,
                   "SMC_Pericytes"  ,
                   "Immune_cells",
                   "Schwann",
                   "Proliferating_cells",
                   "Melanocytes" ,
                   "Hair_follicle_stem_cells" ,
                   "Hair_follicle_matrix_cells",
                   "Sebaceous_gland_channel",
                   "Sebaceous_gland_secretion"  ,
                   "Merkel_cell"  ,
                   "Mucous_cell" ,
                   "Goblet_cell",
                   "Ionocytes" 
  )))

df_long$Species<-gsub('PropMean.species','',df_long$Species)
df_long$Species <- factor(
  df_long$Species,
  levels = rev(c(  "Goldfish" ,
                   "Axolotl"  ,
                   "Gecko" ,
                   "Snake" ,
                   "Crocodile"  ,
                   "Pigeon" ,
                   "Duck" ,
                   "Quail",
                   "Turkey",
                   "Chicken",
                   "Hedgehog"  ,
                   'Cat',
                   "Dog",
                   "Horse" ,
                   "Donkey",
                   "Camel"  ,
                   "Yak",
                   "Cattle",
                   "Buffalo"  ,
                   "Goat"  ,
                   "Beozar.ibex.goat" ,  
                   "Sheep" ,
                   "Asian.mouflon",
                   "Reindeer",
                   "Pig"  ,
                   "Rabbit"  ,
                   "Rat" ,
                   "Mouse",
                   "Human",
                   "Monkey",
                   'Tree.shrew'
  )))



ggplot(df_long,
       aes(x = Species,
           y = Proportion,
           fill = CellType)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  scale_fill_manual(values = cell.cols) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Species",
    y = "Cell proportion",
    fill = "Cell type"
  )
