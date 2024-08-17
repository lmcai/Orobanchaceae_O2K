###############################
#Pairwise sequence divergence
#############################

library(ape)
sequences <- read.dna("your_file.fasta", format = "fasta")
#raw sequence difference p-distance
divergence_matrix <- dist.dna(seqs_dnabin, model = "raw”,pairwise.deletion = TRUE)
divergence_df <- as.data.frame(as.matrix(divergence_matrix))
write.table(divergence_df, file = "distance_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)

#sequence divergence under the F84 model
divergence_matrix <- dist.dna(seqs_dnabin, model = "F84”,pairwise.deletion = TRUE)
divergence_df <- as.data.frame(as.matrix(divergence_matrix))
write.table(divergence_df, file = "distance_matrix_F84.tsv", sep = "\t", quote = FALSE, col.names = NA)


#############################
#Phylogenetic ANOVA
#############################

#Prepare data input

#trait data
library(caper)
data <- read.csv("seq_div_missing_gene.csv")
data$Class_oro=as.factor(data$Class_oro)
data$Class_para=as.factor(data$Class_para)

#tree
phy_tree <- read.tree("Nmt.con.rooted.tre")
phy_tree$node.label<-NULL
test_data <- comparative.data(phy_tree, data, Species)

#Phylogenetic anova test
#make tree and trait data consistant
phylANOVA(phy_tree, data$Class_oro, data$raw, nsim=1000, posthoc=TRUE, p.adj="holm")

##########################
#Raw sequence divergence
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq  F value Pr(>F)
x        0.003043 0.003043 1.086074  0.519
Residual 0.106469 0.002802                

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0         1
0 0.000000 -1.042149
1 1.042149  0.000000

Pairwise corrected P-values:
      0     1
0 1.000 0.519
1 0.519 1.000
---------


###########################
#corrected divergence under the F84 model
phylANOVA(phy_tree, data$Class_oro, data$F84, nsim=1000, posthoc=TRUE, p.adj="holm")

Warning: no labels for x. Assuming order of tree$tip.label.

Warning: no labels for y. Assuming order of tree$tip.label.

ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq  F value Pr(>F)
x        0.006351 0.006351 1.209861  0.548
Residual 0.199463 0.005249                

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0         1
0 0.000000 -1.099937
1 1.099937  0.000000

Pairwise corrected P-values:
      0     1
0 1.000 0.548
1 0.548 1.000
---------

#########################
#plant busco in parasite
phylANOVA(phy_tree, data$Class_para, data$missing_plant_busco, nsim=1000, posthoc=TRUE, p.adj="holm")
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq    Mean Sq   F value Pr(>F)
x        371584.4 371584.411 64.552049  0.001
Residual 224497.8   5756.353                 

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0         1
0 0.000000 -8.034429
1 8.034429  0.000000

Pairwise corrected P-values:
      0     1
0 1.000 0.001
1 0.001 1.000
---------

#########################
#eukaryote busco in parasite

ANOVA table: Phylogenetic ANOVA

Response: y
            Sum Sq   Mean Sq F value Pr(>F)
x         142.9734 142.97336 2.86793  0.266
Residual 1944.2462  49.85247               

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0         1
0 0.000000 -1.693496
1 1.693496  0.000000

Pairwise corrected P-values:
      0     1
0 1.000 0.266
1 0.266 1.000
---------

