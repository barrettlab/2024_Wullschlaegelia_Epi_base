# 2024_Wullschlaegelia_Epi_base
2024_Wullschlaegelia_Epi_base


## Topology tests

```bash

## Run analyses of unconstrained and constrained data
nohup iqtree2 -s Wull_plastomes_72taxa.fasta -m GTR+FO*H6 --prefix WGonly__H6_plastid -T 32 -o Mapania-palustris &
nohup iqtree2 -s Wull_plastomes_72taxa.fasta -m GTR+FO*H6 -g plastid_new_constraint.tre --prefix WGonly_constrained2_H6_plastid -T 32 -o Mapania-palustris &

nohup iqtree2 -s Mito_aligned_w_Pogoniopsis2.fasta -m GTR+FO*H4 -B 1000 --prefix WGonly_H4_mito -T 32 -o Phoenix_dactylifera &
nohup iqtree2 -s Mito_aligned_w_Pogoniopsis2.fasta -m GTR+FO*H4 -g mito_constraint2.tre --prefix WGonly_constrained2_H4_mito -T 32 -o Phoenix_dactylifera &

## Conduct S-H and AU tests
iqtree2 -s Wull_plastomes_72taxa.fasta -z plastid_trees_test.tre -T 32 -n 0 -zb 10000 -zw -au
iqtree2 -s Mito_aligned_w_Pogoniopsis2.fasta -z mito_trees_topology_test.tre -T 32 -n 0 -zb 10000 -zw -au

```

## Analyses with CoEvol

## Run COEVOL in a UNIX terminal
## need 3 files: tree, alignment, traits

### base_tree.tre
### coevol_wull_concat.phy
### coevol_traits.txt, formatted as below

```bash
#      #TRAITS             
#       6 10 length cds irlength tRNA rttdist gc rep20 rep20dens msat msatdens DCJ
#      Aphyllorchis-montana-KU551262 94559 34 0.001 30 0.1084235 37.1 32 0.033841305 31 0.032783765 0.0001  
#      Dendrobium-acinaciforme-LC636124 150841 67 26316 30 0.1079295 36.4 90 0.059665476 31 0.020551442 0.0001  
#      Didymoplexis-pallens-ON515488 51241 17 21541 14 0.4829715 33.3 20 0.039031245 17 0.033176558 0.0001  
#      Epipogium-aphyllum-KJ946456 30650 17 11267 6 0.6526735 30.2 8 0.026101142 22 0.07177814 1  
#      Epipogium-roseum-KJ946455 18966 18 0.001 7 0.7541375 31 18 0.094906675 15 0.079088896 2  
#      Gastrodia-elata-MN200389 35304 20 0.001 5 0.5295795 26.8 22 0.062315885 38 0.107636528 2  
#      Holcoglossum-flavescens-MK442925 146863 69 25808 30 0.1169135 35.3 63 0.042897122 53 0.036088055 0.0001  
#      Nervilia-fordii-ON515491 162651 80 26449 30 0.1202415 35.4 113 0.069473904 50 0.030740666 0.0001  
#      NY1-Wullschlaegelia-aphylla-PR-final 36750 23 0.001 7 0.2629925 33.9 16 0.043537415 18 0.048979592 1  
#      Pogoniopsis-schenckii-OP425397 14015 12 1508 2 0.9713255 23.7 4 0.028540849 16 0.114163396 2  
#      Sobralia-callosa-KM032623 161430 79 26985 30 0.0992505 35.9 57 0.035309422 47 0.029114787 0.0001  
#      Spiranthes-sinensis-MK936427 152786 78 25701 21 0.050842 34.8 80 0.052360818 99 0.064796513 0.0001  
#      Triphora-trianthophoros 151899 70 34260 30 0.1291115 36.1 64 0.04213326 27 0.017774969 0.0001  
```

## Run 3 instances (use screen or nohup &)

```bash
nohup coevol -dsdn -gc -univ -d coevol_concat_wull.phy -t base_tree.tre -c traits_coevol_wull.txt run1 &
nohup coevol -dsdn -gc -univ -d coevol_concat_wull.phy -t base_tree.tre -c traits_coevol_wull.txt run2 &
nohup coevol -dsdn -gc -univ -d coevol_concat_wull.phy -t base_tree.tre -c traits_coevol_wull.txt run3 &
```

## periodically check the runs for convergence

```bash
tracecomp -x 10000 run1 run2 run3
```

## create the correlation table

```bash
readcoevol -x 100 run1
```

# In R, corrplot with GGally/ggpairs

```{r}

traits <- read.csv("traits_coevol.csv",header=T, row.names = 1)

## exclude grouping variable "trophic"

library(GGally)
library(ggplot2)

pm <- ggpairs(traits, columns = 2:10, ggplot2::aes(colour = trophic))

pm + scale_color_manual(values = c("blue", "red")) + scale_fill_manual(values = c("blue", "red"))

```

# Hyphy analyses

```bash

## In unix, remove all stop codons from fasta files. Will add "_macse_NT" to all

macse_v1.2.jar -prog alignSequences -seq z-accD.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-clpP.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-infA.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-matK.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rpl14.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rpl16.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rpl20.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rpl23.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rpl2.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rpl36.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps11.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps12.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps14.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps15.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps18.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps19.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps2.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps3.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps4.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps7.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-rps8.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-ycf1_aligned.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog alignSequences -seq z-ycf2_aligned.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;

macse_v1.2.jar -prog exportAlignment  -align z-accD_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-clpP_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-infA_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-matK_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rpl14_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rpl16_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rpl20_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rpl23_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rpl2_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rpl36_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps11_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps12_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps14_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps15_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps18_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps19_macse_NT.fasta  -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps2_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps3_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps4_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps7_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-rps8_macse_NT.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-ycf1_aligned.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment  -align z-ycf2_aligned.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;

macse_v1.2.jar -prog exportAlignment  -align z-matK_new.fasta   -codonForFinalStop --- -codonForInternalStop NNN ;


### OR

for f in *.fasta
	do
		macse_v1.2.jar -prog alignSequences -seq {f}.fasta   -codonForFinalStop --- -codonForInternalStop NNN
	done

```

```{r}

library(phylotools)
  
accD   <- read.dna(file = "z-accD_macse_NT_macse_NT.fasta",  format = "fasta", as.character = T)
clpP   <- read.dna(file = "z-clpP_macse_NT_macse_NT.fasta",  format = "fasta", as.character = T)
infA   <- read.dna(file = "z-infA_macse_NT_macse_NT.fasta",  format = "fasta", as.character = T)
matK   <- read.dna(file = "z-matK_macse_NT_macse_NT.fasta",  format = "fasta", as.character = T)
rpl14  <- read.dna(file = "z-rpl14_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rpl16  <- read.dna(file = "z-rpl16_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rpl20  <- read.dna(file = "z-rpl20_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rpl23  <- read.dna(file = "z-rpl23_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rpl2   <- read.dna(file = "z-rpl2_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rpl36  <- read.dna(file = "z-rpl36_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps11  <- read.dna(file = "z-rps11_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps12  <- read.dna(file = "z-rps12_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps14  <- read.dna(file = "z-rps14_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps15  <- read.dna(file = "z-rps15_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps18  <- read.dna(file = "z-rps18_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps19  <- read.dna(file = "z-rps19_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps2   <- read.dna(file = "z-rps2_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps3   <- read.dna(file = "z-rps3_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps4   <- read.dna(file = "z-rps4_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps7   <- read.dna(file = "z-rps7_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
rps8   <- read.dna(file = "z-rps8_macse_NT_macse_NT.fasta", format = "fasta", as.character = T)
ycf1   <- read.dna(file = "z-ycf1_aligned_macse_NT.fasta", format = "fasta", as.character = T)
ycf2   <- read.dna(file = "z-ycf2_aligned_macse_NT.fasta", format = "fasta", as.character = T)
matK   <- read.dna(file = "z-matK_new_macse_NT.fasta",  format = "fasta", as.character = T)


tmp <- name.check(basetree,clpP)

infAtree<-drop.tip(basetree,c("Gastrodia-elata-MN200389","Pogoniopsis-schenckii-OP425397")
write.tree(infAtree,file="infAtree.tre")

matKtree<-drop.tip(basetree,c("Didymoplexis-pallens-ON515488","Epipogium-aphyllum-KJ946456","Epipogium-roseum-KJ946455","Pogoniopsis-schenckii-OP425397"))
write.tree(matKtree,file="matKtree.tre")

rpl14tree<-drop.tip(basetree,c("Pogoniopsis-schenckii-OP425397"))
write.tree(rpl14tree,file="rpl14tree.tre")

rpl16tree<-drop.tip(basetree,c("Pogoniopsis-schenckii-OP425397"))
write.tree(rpl16tree,file="rpl16tree.tre")

rpl20tree<-drop.tip(basetree,c("Didymoplexis-pallens-ON515488","Epipogium-aphyllum-KJ946456","Pogoniopsis-schenckii-OP425397"))
write.tree(rpl20tree,file="rpl20tree.tre")

rpl23tree<-drop.tip(basetree,c("Didymoplexis-pallens-ON515488","Epipogium-aphyllum-KJ946456","Epipogium-roseum-KJ946455","Gastrodia-elata-MN200389","Pogoniopsis-schenckii-OP425397"))
write.tree(rpl23tree,file="rpl23tree.tre")

rpl2tree<-drop.tip(basetree,c("Dendrobium-acinaciforme-LC636124"))
write.tree(rpl2tree,file="rpl2tree.tre")

rpl36tree<-drop.tip(basetree,c("Pogoniopsis-schenckii-OP425397"))
write.tree(rpl36tree,file="rpl36tree.tre")

rps12tree<-drop.tip(basetree,c("Dendrobium-acinaciforme-LC636124"))
write.tree(rps12tree,file="rps12tree.tre")

rps15tree<-drop.tip(basetree,c("Didymoplexis-pallens-ON515488","Epipogium-aphyllum-KJ946456","Epipogium-roseum-KJ946455","Gastrodia-elata-MN200389","Pogoniopsis-schenckii-OP425397"))
write.tree(rps15tree,file="rps15tree.tre")

rps19tree<-drop.tip(basetree,c("Triphora-trianthophoros"))
write.tree(rps19tree,file="rps19tree.tre")

rps2tree<-drop.tip(basetree,c("Didymoplexis-pallens-ON515488","Pogoniopsis-schenckii-OP425397"))
write.tree(rps2tree,file="rps2tree.tre")

ycf1tree<-drop.tip(basetree,c("Dendrobium-acinaciforme-LC636124","Epipogium-aphyllum-KJ946456","Epipogium-roseum-KJ946455","Pogoniopsis-schenckii-OP425397"))
write.tree(ycf1tree,file="ycf1tree.tre")

ycf2tree<-drop.tip(basetree,c("Epipogium-aphyllum-KJ946456","Epipogium-roseum-KJ946455","Pogoniopsis-schenckii-OP425397"))
write.tree(ycf2tree,file="ycf2tree.tre")


tmp <- name.check(clPptree,clpP)
tmp


tre <- read.tree("H4-boot-hyphy.tre")

tmp <- name.check(tre, clpP)

```

## Run Hyphy

```bash

conda activate hyphy

## need to re-insert '{myco}' after ':' following species names for all test branches in each tree!!!
hyphy relax --alignment clpP-new-alignment_macse_NT.fasta --tree clpP.tre

hyphy relax --alignment rpl14-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
[ERROR]

hyphy relax --alignment z-accD_macse_NT_macse_NT.fasta  --tree base_tree.tre
hyphy relax --alignment z-clpP_macse_NT_macse_NT.fasta  --tree base_tree.tre
hyphy relax --alignment z-infA_macse_NT_macse_NT.fasta  --tree infAtree.tre
hyphy relax --alignment z-matK_new_macse_NT.fasta  --tree matKtree.tre
hyphy relax --alignment z-rpl14_macse_NT_macse_NT.fasta --tree rpl14tree.tre
hyphy relax --alignment z-rpl16_macse_NT_macse_NT.fasta --tree rpl16tree.tre
hyphy relax --alignment z-rpl20_macse_NT_macse_NT.fasta --tree rpl20tree.tre
hyphy relax --alignment z-rpl23_macse_NT_macse_NT.fasta --tree rpl23tree.tre
hyphy relax --alignment z-rpl2_macse_NT_macse_NT.fasta  --tree rpl2tree.tre
hyphy relax --alignment z-rpl36_macse_NT_macse_NT.fasta --tree rpl36tree.tre
hyphy relax --alignment z-rps11_macse_NT_macse_NT.fasta --tree base_tree.tre
hyphy relax --alignment z-rps12_macse_NT_macse_NT.fasta --tree rps12tree.tre
hyphy relax --alignment z-rps14_macse_NT_macse_NT.fasta --tree base_tree.tre
hyphy relax --alignment z-rps15_macse_NT_macse_NT.fasta --tree rps15tree.tre
hyphy relax --alignment z-rps18_macse_NT_macse_NT.fasta --tree base_tree.tre
hyphy relax --alignment z-rps19_macse_NT_macse_NT.fasta --tree rps19tree.tre
hyphy relax --alignment z-rps2_macse_NT_macse_NT.fasta  --tree rps2tree.tre
hyphy relax --alignment z-rps3_macse_NT_macse_NT.fasta  --tree base_tree.tre
hyphy relax --alignment z-rps4_macse_NT_macse_NT.fasta  --tree base_tree.tre
hyphy relax --alignment z-rps7_macse_NT_macse_NT.fasta  --tree base_tree.tre
hyphy relax --alignment z-rps8_macse_NT_macse_NT.fasta  --tree base_tree.tre
hyphy relax --alignment z-ycf1_aligned_macse_NT.fasta   --tree ycf1tree.tre
hyphy relax --alignment z-ycf2_aligned_macse_NT.fasta   --tree ycf2tree.tre

```

## Optional plotting in R

```{r}

library(ggplot2)
library(tidyr)

# read in data
k <- read.csv("k_plotting.csv",header=T)

# lock in order
k$species <- factor(k$species, levels = k$species)

# convert wide to long format
pcm <- gather(k, gene, value, accD:ycf2, factor_key=TRUE)
head(pcm)

# plot, flip axes (log-k values)

 xx = ggplot(pcm, aes(x = species, y = gene)) + 
     geom_point(aes(size = log(value),fill= ifelse( log(value) < 0, "relaxed", "intensified")), alpha = 1, shape = 21) + 
     #scale_size_continuous() + 
     labs( x= "", y = "", size = "k-value", fill = "")  + 
     theme(legend.key=element_blank(), 
           axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
           axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
           legend.text = element_text(size = 10, face ="bold", colour ="black"), 
           legend.title = element_text(size = 12, face = "bold"), 
           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
           legend.position = "right") +  
     scale_y_discrete() + 
	 scale_fill_manual(name="K", values = c("dodgerblue","red")) +   
	 coord_flip()

xx +   facet_grid(. ~ trophic, scales="free_x", space="free_x")

```

# matK intron plots in R

```{r}

library(ggplot2)
library(data.table)

# Read in data
matk <- read.csv("matK_intron_retention.csv",header=T)

# convert wide to long format
long <- melt(matk, id.vars = "Species", variable.name = "group_II_intron")

# Lock in y axis order
long$Species <- factor(long$Species)

# plot with geom_tile, with grid, 
ggplot(long, aes(x = group_II_intron, y = fct_inorder(Species), fill = value)) +
    geom_tile(color="black") + 
	scale_fill_viridis(discrete = TRUE) + 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

```





