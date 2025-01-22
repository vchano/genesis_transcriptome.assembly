#setwd("/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL")

library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggthemes)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(ggVennDiagram)
library(pheatmap)
#library(ComplexHeatmap)
library(dendextend)
library(dplyr)
library(ggplotify)
library(patchwork)
library(UpSetR)

# (a)
length <- read.table("1a.seqs.by.length.txt", header = TRUE)
df.length <- as.data.frame(length)
f1a <- ggplot(data = length, aes(x = Length, y = Seqs)) +
  theme_bw() +
  labs(title = "(a)", x = "Length (Avg.: 1154; Total Symbols: 59935272)", y = "Number of sequences") +
  geom_bar(stat = "identity", fill = "tomato2", width = 5) +
  theme(
    plot.title = element_text(size = 1.4 * 20),  # increased font size by 40%
    axis.title.x = element_text(size = 1.4 * 15),
    axis.title.y = element_text(size = 1.4 * 15),
    axis.text = element_text(size = 1.4 * 10),
    panel.border = element_blank(),  # remove the frame
    axis.line = element_line(color = "black")
  )
f1a
# (b)
annot.dist <- read.table("1b.data.distrib.txt", header = TRUE)

# Create the pie chart
f1b <- ggplot(annot.dist, aes(x = "", y = Seqs, fill = Annotation_progress)) +
  geom_bar(stat = "identity", color = NA) +  # Remove the border by setting color to NA
  coord_polar("y", start = 0) +
  labs(title = "(b)") +
  theme_void() +  # Use theme_void() to remove all background elements
  theme(
    plot.title = element_text(size = 1.4 * 20),  # increased font size by 40%
    legend.title = element_text(size = 1.4 * 10),  # increased font size by 40%
    legend.text = element_text(size = 1.4 * 10),  # increased font size by 40%
    legend.key = element_blank()  # Remove the legend key frame
  ) +
  geom_text(
    aes(x = 1.6, label = c("No hits\n25056 (48.3%)", "Blast hits\n1426 (2.8%)",
                           "GO mapping\n8819 (17%)", "GO-Slim\n16603 (32%)")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    angle = c(-90, 0, 40, -45),
    size = 1.4 * 3.5  # increased font size by 40%
  ) +
  scale_fill_discrete(name = "Annotation Progress",  # Replace underscore with whitespace in legend title
                      labels = function(x) gsub("_", " ", x))  # Replace underscores with spaces in legend labels
f1b
# (c)
go.dist <- read.table("1c.GO.level.distrib.txt", header = TRUE)

# Replace Component, Function, and Process with CC, MF, and BP
go.dist$Category <- gsub("Component", "CC", go.dist$Category)
go.dist$Category <- gsub("Function", "MF", go.dist$Category)
go.dist$Category <- gsub("Process", "BP", go.dist$Category)

f1c <- ggplot(go.dist, aes(GO_Level, Annotations, fill = Category)) +
  geom_bar(width = 0.5, stat = 'identity', position = position_dodge(0.7)) +
  theme_bw() +
  labs(title = "(c)", x = "GO Level", y = "Number of annotations") +
  theme(
    plot.title = element_text(size = 1.4 * 20),  # increased font size by 40%
    axis.title.x = element_text(size = 1.4 * 15),
    axis.title.y = element_text(size = 1.4 * 15),
    axis.text.x = element_text(angle = 0, size = 1.4 * 10),
    axis.text.y = element_text(size = 1.4 * 10),
    legend.title = element_text(size = 1.4 * 10),
    legend.text = element_text(size = 1.4 * 10),
    panel.border = element_blank(),  # remove the frame
    axis.line = element_line(color = "black")
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limit = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"),
                   labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")) +
  scale_fill_brewer('Categories', palette = 'Set2')
f1c
## Arrange the plots
combined_plot <- grid.arrange(f1a, f1b, f1c, ncol = 2, layout_matrix = rbind(c(1, 2), c(3, 3)))

# Save the combined figure
ggsave("ulmi.sci.data.figure.1.transcriptome.tiff", plot = combined_plot, width = 18, height = 12, dpi = 300, units = "in")

## FIGURE 1. A) Seqs vs length. B) Annottation data distribution. C) GO-level distribution
# (a)
length<-read.table("1a.seqs.by.length.txt",header=T)
df.length<-as.data.frame(length)
tiff(file="ulmi.figure.1a.my_seqs.vs.length.tiff",width=8,height=8,units="in",res=300)
ggplot(data=length, aes(x=Length, y=Seqs)) +
  theme_bw() + labs(#title="(a)",
    x="Lenght (Average: 1,154 bp)",y="Number of sequences") +
  geom_bar(stat="identity",fill="tomato2",width=5)+
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size = 15))
dev.off()

# (b)
annot.dist<-read.table("1b.data.distrib.txt",header=T) #,check.names = F)
#f2b<-xx
tiff(file="ulmi.figure.1.b.my_piechart_annot.tiff",width=8,height=8,units="in",res=300)
ggplot(annot.dist, aes(x = "", y = Seqs, fill = Annotation_progress)) +
  geom_bar(stat = "identity", color = "black") +
  coord_polar("y", start=0) + theme_bw() + #labs(title="(b)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(),
        axis.title.y=element_blank(), 
        axis.title = element_text(size = 15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        axis.text.y=element_blank(),strip.text = element_text(size = 15),
        axis.line.y = element_blank()) +
  geom_text(aes(x = 1.6,label = c("No hits\n25056 (48.3%)","Blast hits\n1426 (2.8%)",
                                  "GO mapping\n8819 (17%)","GO-Slim\n16603 (32%)")),
            position = position_stack(vjust = 0.5,reverse = F),
            angle = c(-90, 0, 40, -45),size=6) 
dev.off()

# (c)
go.dist<-read.table("1c.GO.level.distrib.txt",header=T)
tiff(file="ulmi.figure.1c.my_go_level_barplot.tiff",width=15,height=8,units="in",res=300)
ggplot(go.dist, aes(GO_Level, Annotations, fill=Category)) + 
  geom_bar(width=.5, stat='identity', position=position_dodge(.7)) + 
  theme_bw() + labs(#title="(c)",
    x="GO Level",y="Number of annotations") +
  #theme(axis.text.x = element_text(angle = 0, size = 10))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limit = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),
                   labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))+
  scale_fill_brewer('Categories', palette='Set2')+
  theme(axis.text=element_text(size=15),
        axis.text.x = element_text(angle = 0, size = 12),
        axis.title = element_text(size = 15),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
dev.off()

### Grid plots
## Move to a new page
#grid.newpage()
## Create layout : nrow = 2, ncol = 2
#pushViewport(viewport(layout = grid.layout(2, 2)))
## A helper function to define a region on the layout
#define_region <- function(row, col){viewport(layout.pos.row = row, layout.pos.col = col)} 
## Arrange the plots
#print(f2a, vp=define_region(1,1))
#print(f2b, vp = define_region(1, 2))
#print(f2c, vp = define_region(2, 1:2))

##### ##### LOAD LOCAL RESPONSE DATA FROM THE NEW TABLE         ####

local<-read.table("/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL/ulmi.counts.local.dm.04.txt",header=TRUE, sep="\t",dec=".")

# BY USING DESEQ2, WE CREATE THE ANALYSIS DESIGN, TRANSFORM THE DATA, AND PLOT PCA AND HEATMAP-DISTANCES          ####

coldata.local<-c(#samples 
  "MDV1.1","MDV1.3","MDV1.36","MDV1.47","MDV1.4","MDV1.20","MDV1.22","MDV1.43",
  "MDV1.2","MDV1.11","MDV1.34","MDV1.42","MDV1.8","MDV1.21","MDV1.37","MDV1.39",
  "MDV1.13","MDV1.15","MDV1.23","MDV1.50","MDV1.16","MDV1.33","MDV1.41","MDV1.44",
  "MDV1.18","MDV1.31","MDV1.46","MDV1.49","MDV1.9","MDV1.12","MDV1.24","MDV1.45",
  "MDV2.3.8","MDV2.3.24","MDV2.3.29","MDV2.3.38","MDV2.3.5","MDV2.3.6","MDV2.3.7","MDV2.3.43",
  "MDV2.3.3","MDV2.3.11","MDV2.3.41","MDV2.3.47","MDV2.3.10","MDV2.3.27","MDV2.3.40","MDV2.3.42",
  "MDV2.3.18","MDV2.3.21","MDV2.3.33","MDV2.3.49","MDV2.3.13","MDV2.3.17","MDV2.3.30","MDV2.3.45",
  "MDV2.3.15","MDV2.3.28","MDV2.3.32","MDV2.3.44","MDV2.3.19","MDV2.3.26","MDV2.3.46","MDV2.3.50",
  "VAD2.4","VAD2.6","VAD2.7","VAD2.37","VAD2.1","VAD2.11","VAD2.19","VAD2.42",
  "VAD2.2","VAD2.9","VAD2.27","VAD2.41","VAD2.8","VAD2.36","VAD2.40","VAD2.43",
  "VAD2.14","VAD2.15","VAD2.29","VAD2.44","VAD2.25","VAD2.26","VAD2.28","VAD2.45",
  "VAD2.17","VAD2.22","VAD2.33","VAD2.48","VAD2.18","VAD2.23","VAD2.35","VAD2.47",
  #genotypes
  "MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1",
  "MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1",
  "MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3",
  "MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3",
  "VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2",
  "VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2",
  # condition per genotype and sample
  "control","control","control","control","control","control","control","control",
  "control","control","control","control","control","control","control","control",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "control","control","control","control","control","control","control","control",
  "control","control","control","control","control","control","control","control",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "control","control","control","control","control","control","control","control",
  "control","control","control","control","control","control","control","control",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  # time per condition, genotype and sample
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  # GENOTYPE AND TREATMENT
  "MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR",
  "MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR",
  "MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO",
  "MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO",
  "MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR",
  "MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR",
  "MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO",
  "MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO",
  "VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR",
  "VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR",
  "VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO",
  "VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO",
  # GENOTYPE AND TIME
  "MDV1_6","MDV1_6","MDV1_6","MDV1_6","MDV1_24","MDV1_24","MDV1_24","MDV1_24",
  "MDV1_72","MDV1_72","MDV1_72","MDV1_72","MDV1_144","MDV1_144","MDV1_144","MDV1_144",
  "MDV1_6","MDV1_6","MDV1_6","MDV1_6","MDV1_24","MDV1_24","MDV1_24","MDV1_24",
  "MDV1_72","MDV1_72","MDV1_72","MDV1_72","MDV1_144","MDV1_144","MDV1_144","MDV1_144",
  "MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_24","MDV2.3_24","MDV2.3_24","MDV2.3_24",
  "MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_144","MDV2.3_144","MDV2.3_144","MDV2.3_144",
  "MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_24","MDV2.3_24","MDV2.3_24","MDV2.3_24",
  "MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_144","MDV2.3_144","MDV2.3_144","MDV2.3_144",
  "VAD2_6","VAD2_6","VAD2_6","VAD2_6","VAD2_24","VAD2_24","VAD2_24","VAD2_24",
  "VAD2_72","VAD2_72","VAD2_72","VAD2_72","VAD2_144","VAD2_144","VAD2_144","VAD2_144",
  "VAD2_6","VAD2_6","VAD2_6","VAD2_6","VAD2_24","VAD2_24","VAD2_24","VAD2_24",
  "VAD2_72","VAD2_72","VAD2_72","VAD2_72","VAD2_144","VAD2_144","VAD2_144","VAD2_144",
  # TREATMENT AND TIME
  "CTR_6","CTR_6","CTR_6","CTR_6","CTR_24","CTR_24","CTR_24","CTR_24",
  "CTR_72","CTR_72","CTR_72","CTR_72","CTR_144","CTR_144","CTR_144","CTR_144",
  "INO_6","INO_6","INO_6","INO_6","INO_24","INO_24","INO_24","INO_24",
  "INO_72","INO_72","INO_72","INO_72","INO_144","INO_144","INO_144","INO_144",
  "CTR_6","CTR_6","CTR_6","CTR_6","CTR_24","CTR_24","CTR_24","CTR_24",
  "CTR_72","CTR_72","CTR_72","CTR_72","CTR_144","CTR_144","CTR_144","CTR_144",
  "INO_6","INO_6","INO_6","INO_6","INO_24","INO_24","INO_24","INO_24",
  "INO_72","INO_72","INO_72","INO_72","INO_144","INO_144","INO_144","INO_144",
  "CTR_6","CTR_6","CTR_6","CTR_6","CTR_24","CTR_24","CTR_24","CTR_24",
  "CTR_72","CTR_72","CTR_72","CTR_72","CTR_144","CTR_144","CTR_144","CTR_144",
  "INO_6","INO_6","INO_6","INO_6","INO_24","INO_24","INO_24","INO_24",
  "INO_72","INO_72","INO_72","INO_72","INO_144","INO_144","INO_144","INO_144")

coldata.local<-matrix(coldata.local,nrow=96,ncol=7,byrow=FALSE)

rownames(coldata.local)<-c("MDV1.1_C_L_6","MDV1.3_C_L_6","MDV1.36_C_L_6","MDV1.47_C_L_6",
                           "MDV1.4_C_L_24","MDV1.20_C_L_24","MDV1.22_C_L_24","MDV1.43_C_L_24",
                           "MDV1.2_C_L_72","MDV1.11_C_L_72","MDV1.34_C_L_72","MDV1.42_C_L_72",
                           "MDV1.8_C_L_144","MDV1.21_C_L_144","MDV1.37_C_L_144","MDV1.39_C_L_144",
                           "MDV1.13_O_L_6","MDV1.15_O_L_6","MDV1.23_O_L_6","MDV1.50_O_L_6",
                           "MDV1.16_O_L_24","MDV1.33_O_L_24","MDV1.41_O_L_24","MDV1.44_O_L_24",
                           "MDV1.18_O_L_72","MDV1.31_O_L_72","MDV1.46_O_L_72","MDV1.49_O_L_72",
                           "MDV1.9_O_L_144","MDV1.12_O_L_144","MDV1.24_O_L_144","MDV1.45_O_L_144",
                           "MDV2.3.8_C_L_6","MDV2.3.24_C_L_6","MDV2.3.29_C_L_6","MDV2.3.38_C_L_6",
                           "MDV2.3.5_C_L_24","MDV2.3.6_C_L_24","MDV2.3.7_C_L_24","MDV2.3.43_C_L_24",
                           "MDV2.3.3_C_L_72","MDV2.3.11_C_L_72","MDV2.3.41_C_L_72","MDV2.3.47_C_L_72",
                           "MDV2.3.10_C_L_144","MDV2.3.27_C_L_144","MDV2.3.40_C_L_144","MDV2.3.42_C_L_144",
                           "MDV2.3.18_O_L_6","MDV2.3.21_O_L_6","MDV2.3.33_O_L_6","MDV2.3.49_O_L_6",
                           "MDV2.3.13_O_L_24","MDV2.3.17_O_L_24","MDV2.3.30_O_L_24","MDV2.3.45_O_L_24",
                           "MDV2.3.15_O_L_72","MDV2.3.28_O_L_72","MDV2.3.32_O_L_72","MDV2.3.44_O_L_72",
                           "MDV2.3.19_O_L_144","MDV2.3.26_O_L_144","MDV2.3.46_O_L_144","MDV2.3.50_O_L_144",
                           "VAD2.4_C_L_6","VAD2.6_C_L_6","VAD2.7_C_L_6","VAD2.37_C_L_6",
                           "VAD2.1_C_L_24","VAD2.11_C_L_24","VAD2.19_C_L_24","VAD2.42_C_L_24",
                           "VAD2.2_C_L_72","VAD2.9_C_L_72","VAD2.27_C_L_72","VAD2.41_C_L_72",
                           "VAD2.8_C_L_144","VAD2.36_C_L_144","VAD2.40_C_L_144","VAD2.43_C_L_144",
                           "VAD2.14_O_L_6","VAD2.15_O_L_6","VAD2.29_O_L_6","VAD2.44_O_L_6",
                           "VAD2.25_O_L_24","VAD2.26_O_L_24","VAD2.28_O_L_24","VAD2.45_O_L_24",
                           "VAD2.17_O_L_72","VAD2.22_O_L_72","VAD2.33_O_L_72","VAD2.48_O_L_72",
                           "VAD2.18_O_L_144","VAD2.23_O_L_144","VAD2.35_O_L_144","VAD2.47_O_L_144")

colnames(coldata.local)<-c("SAMPLE","GENOTYPE","TREATMENT","TIME","GENOT.TREATMENT","GENOT.TIME","TREAT.TIME")

dds.local <- DESeqDataSetFromMatrix(countData = local,
                                    colData = coldata.local,
                                    design = ~ GENOTYPE + TIME + TREATMENT + TIME:TREATMENT)

## DATA VISUALIZATION: TRANSFORMATION, PCAs, HEATMAPS...
t_data.local<-vst(dds.local)
#t_data.local <- assay(vst(dds.local, blind=FALSE)) # create a table with vst values and save it for filtering later DEGs and make heatmap
#write.table(as.data.frame   (t_data.local), file="ulmi.t_data.vst.local.txt",    sep="\t",dec =".",row.names = TRUE)

head(assay(t_data.local))
distances.local<-dist(t(assay(t_data.local)))
distances_matrix.local<-as.matrix(distances.local)
rownames(distances_matrix.local)<-paste(t_data.local$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.local<-hclust(distances.local)

heatmap.2(distances_matrix.local,Rowv = as.dendrogram(hc.local),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE)

#plotPCA(t_data.local,intgroup="SAMPLE")
plotPCA(t_data.local,intgroup="GENOTYPE")
plotPCA(t_data.local,intgroup="TREATMENT")
plotPCA(t_data.local,intgroup="TIME")
plotPCA(t_data.local,intgroup="GENOT.TREATMENT")
plotPCA(t_data.local,intgroup="GENOT.TIME")
plotPCA(t_data.local,intgroup="TREAT.TIME")

pcaData.l <- plotPCA(t_data.local, intgroup=c("TREATMENT", "GENOTYPE"), returnData=TRUE)
percentVar.l <- round(100 * attr(pcaData.l, "percentVar"))

tiff(file="ulmi.sci.data.figure.2.pca.tiff",width=8,height=8,units="in",res=300)
ggplot(pcaData.l, aes(PC1, PC2, shape=TREATMENT, color=GENOTYPE)) +
  xlab(paste0("PC1: ",percentVar.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.l[2],"% variance")) + 
  geom_point(size=3.5) +# guides(color = FALSE, shape = FALSE) +
  theme_bw() + #labs(title="(a)") + 
  xlim(-15, 15) + ylim(-15, 15) +
  theme(axis.title = element_text(size = 15)) + 
  theme(legend.title=element_text(size=15), 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=15)) +
  coord_fixed()
dev.off()

#####

##### ##### SEPARATE GENOTYPES AND TIMES AND PERFORM PAIRWISE ANALYSIS BY COMPARING CONTROL VS INFECTED PLANTS         ####
#####
##### ##### MDV1  6 hpi       ####

mdv1.l.6<-cbind(local$MDV1.1_C_L_6, local$MDV1.3_C_L_6,  local$MDV1.36_C_L_6, local$MDV1.47_C_L_6,
                local$MDV1.13_O_L_6,local$MDV1.15_O_L_6, local$MDV1.23_O_L_6, local$MDV1.50_O_L_6)

colnames(mdv1.l.6)<-c("MDV1.1_C_L_6","MDV1.3_C_L_6","MDV1.36_C_L_6","MDV1.47_C_L_6",
                      "MDV1.13_O_L_6","MDV1.15_O_L_6","MDV1.23_O_L_6","MDV1.50_O_L_6")
rownames(mdv1.l.6)<-rownames(local)

coldata.mdv1.l.6<-c("MDV1.1",  "MDV1.3",  "MDV1.36",  "MDV1.47","MDV1.13","MDV1.15","MDV1.23","MDV1.50",
                    "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.mdv1.l.6<-matrix(coldata.mdv1.l.6,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv1.l.6)<-colnames(mdv1.l.6)
colnames(coldata.mdv1.l.6)<-c("SAMPLE","TREATMENT")
print(coldata.mdv1.l.6)
dds.mdv1.l.6 <- DESeqDataSetFromMatrix(countData = mdv1.l.6,
                                       colData = coldata.mdv1.l.6,
                                       design = ~ TREATMENT)
dds.mdv1.l.6<-DESeq(dds.mdv1.l.6)
res.mdv1.l.6<-results(dds.mdv1.l.6)
res.mdv1.l.6
resOrdered.mdv1.l.6 <- res.mdv1.l.6[order(res.mdv1.l.6$pvalue),]
summary(res.mdv1.l.6)
sum(res.mdv1.l.6$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.6 <- results(dds.mdv1.l.6, alpha=0.05)
summary(res05.mdv1.l.6)
sum(res05.mdv1.l.6$padj < 0.05, na.rm=TRUE)

#res01.mdv1.l.6 <- results(dds.mdv1.l.6, alpha=0.01)
#summary(res01.mdv1.l.6)
#sum(res01.mdv1.l.6$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.6 = subset(res05.mdv1.l.6, padj<0.05)
print(resSig05.mdv1.l.6)
write.table(res.mdv1.l.6,file="res.mdv1.l.6.txt")
write.table(resSig05.mdv1.l.6,file="resSig05.mdv1.l.6.txt")

# Volcano plot

vp01<-EnhancedVolcano(res.mdv1.l.6,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      subtitle = NULL,
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,
                      col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 6 hpi") + theme(plot.title=element_text(size=20)) 
vp01


##### ##### MDV1  24 hpi       ####

mdv1.l.24<-cbind(local$MDV1.4_C_L_24,  local$MDV1.20_C_L_24, local$MDV1.22_C_L_24, local$MDV1.43_C_L_24,
                 local$MDV1.16_O_L_24, local$MDV1.33_O_L_24, local$MDV1.41_O_L_24, local$MDV1.44_O_L_24)
colnames(mdv1.l.24)<-c("MDV1.4_C_L_24",  "MDV1.20_C_L_24",  "MDV1.22_C_L_24",  "MDV1.43_C_L_24",
                       "MDV1.16_O_L_24", "MDV1.33_O_L_24", "MDV1.41_O_L_24", "MDV1.44_O_L_24")
rownames(mdv1.l.24)<-rownames(local)
coldata.mdv1.l.24<-c("MDV1.4",  "MDV1.20",  "MDV1.22",  "MDV1.43","MDV1.16","MDV1.33","MDV1.41","MDV1.44",
                     "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.mdv1.l.24<-matrix(coldata.mdv1.l.24,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv1.l.24)<-colnames(mdv1.l.24)
colnames(coldata.mdv1.l.24)<-c("SAMPLE","TREATMENT")
print(coldata.mdv1.l.24)

dds.mdv1.l.24 <- DESeqDataSetFromMatrix(countData = mdv1.l.24,
                                        colData = coldata.mdv1.l.24,
                                        design = ~ TREATMENT)
dds.mdv1.l.24<-DESeq(dds.mdv1.l.24)
res.mdv1.l.24<-results(dds.mdv1.l.24)
res.mdv1.l.24
resOrdered.mdv1.l.24 <- res.mdv1.l.24[order(res.mdv1.l.24$pvalue),]
summary(res.mdv1.l.24)
sum(res.mdv1.l.24$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.24 <- results(dds.mdv1.l.24, alpha=0.05)
summary(res05.mdv1.l.24)
sum(res05.mdv1.l.24$padj < 0.05, na.rm=TRUE)

#res01.mdv1.l.24 <- results(dds.mdv1.l.24, alpha=0.01)
#summary(res01.mdv1.l.24)
#sum(res01.mdv1.l.24$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.24 = subset(res05.mdv1.l.24, padj<0.05)
print(resSig05.mdv1.l.24)
write.table(res.mdv1.l.24,file="res.mdv1.l.24.txt")
write.table(resSig05.mdv1.l.24,file="resSig05.mdv1.l.24.txt")

# Volcano plot
vp02<-EnhancedVolcano(res.mdv1.l.24,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 24 hpi") + theme(plot.title=element_text(size=20)) 
vp02

##### ##### MDV1  72 hpi       ####

mdv1.l.72<-cbind(local$MDV1.2_C_L_72,  local$MDV1.11_C_L_72,local$MDV1.34_C_L_72,local$MDV1.42_C_L_72,
                 local$MDV1.18_O_L_72, local$MDV1.31_O_L_72,local$MDV1.46_O_L_72,local$MDV1.49_O_L_72)

colnames(mdv1.l.72)<-c("MDV1.2_C_L_72",  "MDV1.11_C_L_72",  "MDV1.34_C_L_72",  "MDV1.42_C_L_72",
                       "MDV1.18_O_L_72", "MDV1.31_O_L_72", "MDV1.46_O_L_72", "MDV1.49_O_L_72")

rownames(mdv1.l.72)<-rownames(local)

coldata.mdv1.l.72<-c("MDV1.2",  "MDV1.11",  "MDV1.34",  "MDV1.42","MDV1.18","MDV1.31","MDV1.46","MDV1.49",
                     "control","control","control","control","inoculated","inoculated","inoculated","inoculated")

coldata.mdv1.l.72<-matrix(coldata.mdv1.l.72,nrow=8,ncol=2,byrow = FALSE)

rownames(coldata.mdv1.l.72)<-colnames(mdv1.l.72)

colnames(coldata.mdv1.l.72)<-c("SAMPLE","TREATMENT")

print(coldata.mdv1.l.72)

dds.mdv1.l.72 <- DESeqDataSetFromMatrix(countData = mdv1.l.72,
                                        colData = coldata.mdv1.l.72,
                                        design = ~ TREATMENT)
dds.mdv1.l.72<-DESeq(dds.mdv1.l.72)
res.mdv1.l.72<-results(dds.mdv1.l.72)
res.mdv1.l.72
resOrdered.mdv1.l.72 <- res.mdv1.l.72[order(res.mdv1.l.72$pvalue),]
summary(res.mdv1.l.72)
sum(res.mdv1.l.72$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.72 <- results(dds.mdv1.l.72, alpha=0.05)
summary(res05.mdv1.l.72)
sum(res05.mdv1.l.72$padj < 0.05, na.rm=TRUE)

#res01.mdv1.l.72 <- results(dds.mdv1.l.72, alpha=0.01)
#summary(res01.mdv1.l.72)
#sum(res01.mdv1.l.72$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.72 = subset(res05.mdv1.l.72, padj<0.05)
print(resSig05.mdv1.l.72)
write.table(res.mdv1.l.72,file="res.mdv1.l.72.txt")
write.table(resSig05.mdv1.l.72,file="resSig05.mdv1.l.72.txt")

# Volcano plot
vp03<-EnhancedVolcano(res.mdv1.l.72,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 72 hpi") + theme(plot.title=element_text(size=20)) 
vp03

##### ##### MDV1  144 hpi       ####

mdv1.l.144<-cbind(local$MDV1.8_C_L_144,local$MDV1.21_C_L_144,local$MDV1.37_C_L_144,local$MDV1.39_C_L_144,
                  local$MDV1.9_O_L_144,local$MDV1.12_O_L_144,local$MDV1.24_O_L_144,local$MDV1.45_O_L_144)

colnames(mdv1.l.144)<-c("MDV1.8_C_L_144",  "MDV1.21_C_L_144",  "MDV1.37_C_L_144",  "MDV1.39_C_L_144",
                        "MDV1.9_O_L_144", "MDV1.12_O_L_144","MDV1.24_O_L_144","MDV1.45_O_L_144")

rownames(mdv1.l.144)<-rownames(local)

coldata.mdv1.l.144<-c("MDV1.8",  "MDV1.21",  "MDV1.37",  "MDV1.39","MDV1.9","MDV1.12","MDV1.24","MDV1.45",
                      "control","control","control","control","inoculated","inoculated","inoculated","inoculated")

coldata.mdv1.l.144<-matrix(coldata.mdv1.l.144,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv1.l.144)<-colnames(mdv1.l.144)
colnames(coldata.mdv1.l.144)<-c("SAMPLE","TREATMENT")
print(coldata.mdv1.l.144)

dds.mdv1.l.144 <- DESeqDataSetFromMatrix(countData = mdv1.l.144,
                                         colData = coldata.mdv1.l.144,
                                         design = ~ TREATMENT)

dds.mdv1.l.144<-DESeq(dds.mdv1.l.144)
res.mdv1.l.144<-results(dds.mdv1.l.144)
res.mdv1.l.144
resOrdered.mdv1.l.144 <- res.mdv1.l.144[order(res.mdv1.l.144$pvalue),]
summary(res.mdv1.l.144)
sum(res.mdv1.l.144$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.144 <- results(dds.mdv1.l.144, alpha=0.05)
summary(res05.mdv1.l.144)
sum(res05.mdv1.l.144$padj < 0.05, na.rm=TRUE)

#res01.mdv1.l.144 <- results(dds.mdv1.l.144, alpha=0.01)
#summary(res01.mdv1.l.144)
#sum(res01.mdv1.l.144$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.144 = subset(res05.mdv1.l.144, padj<0.05)
print(resSig05.mdv1.l.144)
write.table(res.mdv1.l.144,file="res.mdv1.l.144.txt")
write.table(resSig05.mdv1.l.144,file="resSig05.mdv1.l.144.txt")

# Volcano plot
vp04<-EnhancedVolcano(res.mdv1.l.144,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 144 hpi") + theme(plot.title=element_text(size=20)) 
vp04


##### ##### MDV2.3  6 hpi       ####

mdv2.3.l.6<-cbind(local$MDV2.3.8_C_L_6, local$MDV2.3.24_C_L_6,local$MDV2.3.29_C_L_6,local$MDV2.3.38_C_L_6,
                  local$MDV2.3.18_O_L_6,local$MDV2.3.21_O_L_6,local$MDV2.3.33_O_L_6,local$MDV2.3.49_O_L_6)
colnames(mdv2.3.l.6)<-c("MDV2.3.8_C_L_6","MDV2.3.24_C_L_6","MDV2.3.29_C_L_6","MDV2.3.38_C_L_6",
                        "MDV2.3.18_O_L_6",  "MDV2.3.21_O_L_6",  "MDV2.3.33_O_L_6",  "MDV2.3.49_O_L_6")
rownames(mdv2.3.l.6)<-rownames(local)
coldata.mdv2.3.l.6<-c("MDV2.3.8","MDV2.3.24","MDV2.3.29","MDV2.3.38","MDV2.3.18","MDV2.3.21","MDV2.3.33","MDV2.3.49",
                      "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.mdv2.3.l.6<-matrix(coldata.mdv2.3.l.6,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.6)<-colnames(mdv2.3.l.6)
colnames(coldata.mdv2.3.l.6)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.6)

dds.mdv2.3.l.6 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.6,
                                         colData = coldata.mdv2.3.l.6,
                                         design = ~ TREATMENT)
dds.mdv2.3.l.6<-DESeq(dds.mdv2.3.l.6)
res.mdv2.3.l.6<-results(dds.mdv2.3.l.6)
res.mdv2.3.l.6
resOrdered.mdv2.3.l.6 <- res.mdv2.3.l.6[order(res.mdv2.3.l.6$pvalue),]
summary(res.mdv2.3.l.6)
sum(res.mdv2.3.l.6$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.6 <- results(dds.mdv2.3.l.6, alpha=0.05)
summary(res05.mdv2.3.l.6)
sum(res05.mdv2.3.l.6$padj < 0.05, na.rm=TRUE)

#res01.mdv2.3.l.6 <- results(dds.mdv2.3.l.6, alpha=0.01)
#summary(res01.mdv2.3.l.6)
#sum(res01.mdv2.3.l.6$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.6 = subset(res05.mdv2.3.l.6, padj<0.05)
print(resSig05.mdv2.3.l.6)
write.table(resSig05.mdv2.3.l.6,file="resSig05.mdv2.3.l.6.txt")
write.table(res.mdv2.3.l.6,file="res.mdv2.3.l.6.txt")

# Volcano plot
vp05<-EnhancedVolcano(res.mdv2.3.l.6,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 6 hpi") + theme(plot.title=element_text(size=20)) 
vp05

##### ##### MDV2.3  24 hpi       ####

mdv2.3.l.24<-cbind(local$MDV2.3.5_C_L_24,local$MDV2.3.6_C_L_24,local$MDV2.3.7_C_L_24,local$MDV2.3.43_C_L_24,
                   local$MDV2.3.13_O_L_24, local$MDV2.3.17_O_L_24, local$MDV2.3.30_O_L_24, local$MDV2.3.45_O_L_24)
colnames(mdv2.3.l.24)<-c("MDV2.3.5_C_L_24","MDV2.3.6_C_L_24","MDV2.3.7_C_L_24","MDV2.3.43_C_L_24",
                         "MDV2.3.15_O_L_24", "MDV2.3.17_O_L_24", "MDV2.3.30_O_L_24", "MDV2.3.45_O_L_24")
rownames(mdv2.3.l.24)<-rownames(local)
coldata.mdv2.3.l.24<-c("MDV2.3.5","MDV2.3.6","MDV2.3.7","MDV2.3.43","MDV2.3.13","MDV2.3.17","MDV2.3.30","MDV2.3.45",
                       "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.mdv2.3.l.24<-matrix(coldata.mdv2.3.l.24,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.24)<-colnames(mdv2.3.l.24)
colnames(coldata.mdv2.3.l.24)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.24)

dds.mdv2.3.l.24 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.24,
                                          colData = coldata.mdv2.3.l.24,
                                          design = ~ TREATMENT)
dds.mdv2.3.l.24<-DESeq(dds.mdv2.3.l.24)
res.mdv2.3.l.24<-results(dds.mdv2.3.l.24)
res.mdv2.3.l.24
resOrdered.mdv2.3.l.24 <- res.mdv2.3.l.24[order(res.mdv2.3.l.24$pvalue),]
summary(res.mdv2.3.l.24)
sum(res.mdv2.3.l.24$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.24 <- results(dds.mdv2.3.l.24, alpha=0.05)
summary(res05.mdv2.3.l.24)
sum(res05.mdv2.3.l.24$padj < 0.05, na.rm=TRUE)

#res01.mdv2.3.l.24 <- results(dds.mdv2.3.l.24, alpha=0.01)
#summary(res01.mdv2.3.l.24)
#sum(res01.mdv2.3.l.24$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.24 = subset(res05.mdv2.3.l.24, padj<0.05)
print(resSig05.mdv2.3.l.24)
write.table(resSig05.mdv2.3.l.24,file="resSig05.mdv2.3.l.24.txt")
write.table(res.mdv2.3.l.24,file="res.mdv2.3.l.24.txt")

# Volcano plot
vp06<-EnhancedVolcano(res.mdv2.3.l.24,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 24 hpi") + theme(plot.title=element_text(size=20)) 
vp06
##### ##### MDV2.3  72 hpi       ####

mdv2.3.l.72<-cbind(local$MDV2.3.3_C_L_72,local$MDV2.3.11_C_L_72,local$MDV2.3.41_C_L_72,local$MDV2.3.47_C_L_72,
                   local$MDV2.3.15_O_L_72, local$MDV2.3.28_O_L_72, local$MDV2.3.32_O_L_72, local$MDV2.3.44_O_L_72)
colnames(mdv2.3.l.72)<-c("MDV2.3.3_C_L_72","MDV2.3.11_C_L_72","MDV2.3.41_C_L_72","MDV2.3.47_C_L_72",
                         "MDV2.3.15_O_L_72", "MDV2.3.28_O_L_72", "MDV2.3.32_O_L_72", "MDV2.3.44_O_L_72")
rownames(mdv2.3.l.72)<-rownames(local)
coldata.mdv2.3.l.72<-c("MDV2.3.3","MDV2.3.11","MDV2.3.41","MDV2.3.47","MDV2.3.15","MDV2.3.28","MDV2.3.32","MDV2.3.44",
                       "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.mdv2.3.l.72<-matrix(coldata.mdv2.3.l.72,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.72)<-colnames(mdv2.3.l.72)
colnames(coldata.mdv2.3.l.72)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.72)

dds.mdv2.3.l.72 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.72,
                                          colData = coldata.mdv2.3.l.72,
                                          design = ~ TREATMENT)
dds.mdv2.3.l.72<-DESeq(dds.mdv2.3.l.72)
res.mdv2.3.l.72<-results(dds.mdv2.3.l.72)
res.mdv2.3.l.72
resOrdered.mdv2.3.l.72 <- res.mdv2.3.l.72[order(res.mdv2.3.l.72$pvalue),]
summary(res.mdv2.3.l.72)
sum(res.mdv2.3.l.72$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.72 <- results(dds.mdv2.3.l.72, alpha=0.05)
summary(res05.mdv2.3.l.72)
sum(res05.mdv2.3.l.72$padj < 0.05, na.rm=TRUE)

#res01.mdv2.3.l.72 <- results(dds.mdv2.3.l.72, alpha=0.01)
#summary(res01.mdv2.3.l.72)
#sum(res01.mdv2.3.l.72$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.72 = subset(res05.mdv2.3.l.72, padj<0.05)
print(resSig05.mdv2.3.l.72)
write.table(resSig05.mdv2.3.l.72,file="resSig05.mdv2.3.l.72.txt")
write.table(res.mdv2.3.l.72,file="res.mdv2.3.l.72.txt")

# Volcano plot
vp07<-EnhancedVolcano(res.mdv2.3.l.72,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 72 hpi") + theme(plot.title=element_text(size=20)) 
vp07
##### ##### MDV2.3  144 hpi       ####

mdv2.3.l.144<-cbind(local$MDV2.3.10_C_L_144,local$MDV2.3.27_C_L_144,local$MDV2.3.40_C_L_144,local$MDV2.3.42_C_L_144,
                    local$MDV2.3.19_O_L_144,local$MDV2.3.26_O_L_144,local$MDV2.3.46_O_L_144,local$MDV2.3.50_O_L_144)
colnames(mdv2.3.l.144)<-c("MDV2.3.10_C_L_144","MDV2.3.27_C_L_144","MDV2.3.40_C_L_144","MDV2.3.42_C_L_144",
                          "MDV2.3.19_O_L_144","MDV2.3.26_O_L_144","MDV2.3.46_O_L_144","MDV2.3.50_O_L_144")
rownames(mdv2.3.l.144)<-rownames(local)
coldata.mdv2.3.l.144<-c("MDV2.3.10","MDV2.3.27","MDV2.3.40","MDV2.3.42","MDV2.3.19","MDV2.3.26","MDV2.3.46","MDV2.3.50",
                        "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.mdv2.3.l.144<-matrix(coldata.mdv2.3.l.144,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.144)<-colnames(mdv2.3.l.144)
colnames(coldata.mdv2.3.l.144)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.144)

dds.mdv2.3.l.144 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.144,
                                           colData = coldata.mdv2.3.l.144,
                                           design = ~ TREATMENT)
dds.mdv2.3.l.144<-DESeq(dds.mdv2.3.l.144)
res.mdv2.3.l.144<-results(dds.mdv2.3.l.144)
res.mdv2.3.l.144
resOrdered.mdv2.3.l.144 <- res.mdv2.3.l.144[order(res.mdv2.3.l.144$pvalue),]
summary(res.mdv2.3.l.144)
sum(res.mdv2.3.l.144$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.144 <- results(dds.mdv2.3.l.144, alpha=0.05)
summary(res05.mdv2.3.l.144)
sum(res05.mdv2.3.l.144$padj < 0.05, na.rm=TRUE)

#res01.mdv2.3.l.144 <- results(dds.mdv2.3.l.144, alpha=0.01)
#summary(res01.mdv2.3.l.144)
#sum(res01.mdv2.3.l.144$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.144 = subset(res05.mdv2.3.l.144, padj<0.05)
print(resSig05.mdv2.3.l.144)
write.table(resSig05.mdv2.3.l.144,file="resSig05.mdv2.3.l.144.txt")
write.table(res.mdv2.3.l.144,file="res.mdv2.3.l.144.txt")

# Volcano plot
vp08<-EnhancedVolcano(res.mdv2.3.l.144,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 144 hpi") + theme(plot.title=element_text(size=20)) 
vp08
##### ##### VAD2  6 hpi       ####

vad2.l.6<-cbind(local$VAD2.4_C_L_6,local$VAD2.6_C_L_6,local$VAD2.7_C_L_6,local$VAD2.37_C_L_6,
                local$VAD2.14_O_L_6,local$VAD2.15_O_L_6,local$VAD2.29_O_L_6,local$VAD2.44_O_L_6)
colnames(vad2.l.6)<-c("VAD2.4_C_L_6","VAD2.6_C_L_6","VAD2.7_C_L_6","VAD2.37_C_L_6",
                      "VAD2.14_O_L_6",  "VAD2.15_O_L_6",  "VAD2.29_O_L_6",  "VAD2.44_O_L_6")
rownames(vad2.l.6)<-rownames(local)
coldata.vad2.l.6<-c("VAD2.4","VAD2.6","VAD2.7","VAD2.37","VAD2.14","VAD2.15","VAD2.29","VAD2.44",
                    "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.vad2.l.6<-matrix(coldata.vad2.l.6,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.6)<-colnames(vad2.l.6)
colnames(coldata.vad2.l.6)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.6)

dds.vad2.l.6 <- DESeqDataSetFromMatrix(countData = vad2.l.6,
                                       colData = coldata.vad2.l.6,
                                       design = ~ TREATMENT)
dds.vad2.l.6<-DESeq(dds.vad2.l.6)
res.vad2.l.6<-results(dds.vad2.l.6)
res.vad2.l.6
resOrdered.vad2.l.6 <- res.vad2.l.6[order(res.vad2.l.6$pvalue),]
summary(res.vad2.l.6)
sum(res.vad2.l.6$padj < 0.1, na.rm=TRUE)

res05.vad2.l.6 <- results(dds.vad2.l.6, alpha=0.05)
summary(res05.vad2.l.6)
sum(res05.vad2.l.6$padj < 0.05, na.rm=TRUE)

#res01.vad2.l.6 <- results(dds.vad2.l.6, alpha=0.01)
#summary(res01.vad2.l.6)
#sum(res01.vad2.l.6$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.6 = subset(res05.vad2.l.6, padj<0.05)
print(resSig05.vad2.l.6)
write.table(resSig05.vad2.l.6,file="resSig05.vad2.l.6.txt")
write.table(res.vad2.l.6,file="res.vad2.l.6.txt")

# Volcano plot
vp09<-EnhancedVolcano(res.vad2.l.6,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 6 hpi") + theme(plot.title=element_text(size=20)) 
vp09

##### ##### VAD2  24 hpi       ####

vad2.l.24<-cbind(local$VAD2.1_C_L_24,local$VAD2.11_C_L_24,local$VAD2.19_C_L_24,local$VAD2.42_C_L_24,
                 local$VAD2.25_O_L_24,local$VAD2.26_O_L_24,local$VAD2.28_O_L_24,local$VAD2.45_O_L_24)
colnames(vad2.l.24)<-c("VAD2.1_C_L_24","VAD2.11_C_L_24","VAD2.19_C_L_24","VAD2.42_C_L_24",
                       "VAD2.25_O_L_24","VAD2.26_O_L_24","VAD2.28_O_L_24","VAD2.45_O_L_24")
rownames(vad2.l.24)<-rownames(local)
coldata.vad2.l.24<-c("VAD2.1","VAD2.11","VAD2.19","VAD2.42","VAD2.25","VAD2.26","VAD2.28","VAD2.45",
                     "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.vad2.l.24<-matrix(coldata.vad2.l.24,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.24)<-colnames(vad2.l.24)
colnames(coldata.vad2.l.24)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.24)

dds.vad2.l.24 <- DESeqDataSetFromMatrix(countData = vad2.l.24,
                                        colData = coldata.vad2.l.24,
                                        design = ~ TREATMENT)
dds.vad2.l.24<-DESeq(dds.vad2.l.24)
res.vad2.l.24<-results(dds.vad2.l.24)
res.vad2.l.24
resOrdered.vad2.l.24 <- res.vad2.l.24[order(res.vad2.l.24$pvalue),]
summary(res.vad2.l.24)
sum(res.vad2.l.24$padj < 0.1, na.rm=TRUE)

res05.vad2.l.24 <- results(dds.vad2.l.24, alpha=0.05)
summary(res05.vad2.l.24)
sum(res05.vad2.l.24$padj < 0.05, na.rm=TRUE)

#res01.vad2.l.24 <- results(dds.vad2.l.24, alpha=0.01)
#summary(res01.vad2.l.24)
#sum(res01.vad2.l.24$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.24 = subset(res05.vad2.l.24, padj<0.05)
print(resSig05.vad2.l.24)
write.table(resSig05.vad2.l.24,file="resSig05.vad2.l.24.txt")
write.table(res.vad2.l.24,file="res.vad2.l.24.txt")

# Volcano plot
vp10<-EnhancedVolcano(res.vad2.l.24,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 24 hpi") + theme(plot.title=element_text(size=20)) 
vp10
##### ##### VAD2  72 hpi       ####

vad2.l.72<-cbind(local$VAD2.2_C_L_72,local$VAD2.9_C_L_72,local$VAD2.27_C_L_72,local$VAD2.41_C_L_72,
                 local$VAD2.17_O_L_72,local$VAD2.22_O_L_72,local$VAD2.33_O_L_72, local$VAD2.48_O_L_72)
colnames(vad2.l.72)<-c("VAD2.2_C_L_72","VAD2.9_C_L_72","VAD2.27_C_L_72","VAD2.41_C_L_72",
                       "VAD2.17_O_L_72","VAD2.22_O_L_72","VAD2.33_O_L_72","VAD2.48_O_L_72")
rownames(vad2.l.72)<-rownames(local)
coldata.vad2.l.72<-c("VAD2.2","VAD2.9","VAD2.27","VAD2.41","VAD2.17","VAD2.22","VAD2.33","VAD2.48",
                     "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.vad2.l.72<-matrix(coldata.vad2.l.72,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.72)<-colnames(vad2.l.72)
colnames(coldata.vad2.l.72)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.72)

dds.vad2.l.72 <- DESeqDataSetFromMatrix(countData = vad2.l.72,
                                        colData = coldata.vad2.l.72,
                                        design = ~ TREATMENT)
dds.vad2.l.72<-DESeq(dds.vad2.l.72)
res.vad2.l.72<-results(dds.vad2.l.72)
res.vad2.l.72
resOrdered.vad2.l.72 <- res.vad2.l.72[order(res.vad2.l.72$pvalue),]
summary(res.vad2.l.72)
sum(res.vad2.l.72$padj < 0.1, na.rm=TRUE)

res05.vad2.l.72 <- results(dds.vad2.l.72, alpha=0.05)
summary(res05.vad2.l.72)
sum(res05.vad2.l.72$padj < 0.05, na.rm=TRUE)

#res01.vad2.l.72 <- results(dds.vad2.l.72, alpha=0.01)
#summary(res01.vad2.l.72)
#sum(res01.vad2.l.72$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.72 = subset(res05.vad2.l.72, padj<0.05)
print(resSig05.vad2.l.72)
write.table(resSig05.vad2.l.72,file="resSig05.vad2.l.72.txt")
write.table(res.vad2.l.72,file="res.vad2.l.72.txt")

# Volcano plot
vp11<-EnhancedVolcano(res.vad2.l.72,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 72 hpi") + theme(plot.title=element_text(size=20)) 
vp11

##### ##### VAD2  144 hpi       ####

vad2.l.144<-cbind(local$VAD2.8_C_L_144,local$VAD2.36_C_L_144,local$VAD2.40_C_L_144,local$VAD2.43_C_L_144,
                  local$VAD2.18_O_L_144,  local$VAD2.23_O_L_144,  local$VAD2.35_O_L_144,local$VAD2.47_O_L_144)
colnames(vad2.l.144)<-c("VAD2.8_C_L_144","VAD2.36_C_L_144","VAD2.40_C_L_144","VAD2.43_C_L_144",
                        "VAD2.18_O_L_144","VAD2.23_O_L_144","VAD2.35_O_L_144","VAD2.47_O_L_144")
rownames(vad2.l.144)<-rownames(local)
coldata.vad2.l.144<-c("VAD2.8","VAD2.36","VAD2.40","VAD2.43","VAD2.18","VAD2.23", "VAD2.35", "VAD2.47",
                      "control","control","control","control","inoculated","inoculated","inoculated","inoculated")
coldata.vad2.l.144<-matrix(coldata.vad2.l.144,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.144)<-colnames(vad2.l.144)
colnames(coldata.vad2.l.144)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.144)
dds.vad2.l.144 <- DESeqDataSetFromMatrix(countData = vad2.l.144,
                                         colData = coldata.vad2.l.144,
                                         design = ~ TREATMENT)

dds.vad2.l.144<-DESeq(dds.vad2.l.144)
res.vad2.l.144<-results(dds.vad2.l.144)
res.vad2.l.144
resOrdered.vad2.l.144 <- res.vad2.l.144[order(res.vad2.l.144$pvalue),]
summary(res.vad2.l.144)
sum(res.vad2.l.144$padj < 0.1, na.rm=TRUE)

res05.vad2.l.144 <- results(dds.vad2.l.144, alpha=0.05)
summary(res05.vad2.l.144)
sum(res05.vad2.l.144$padj < 0.05, na.rm=TRUE)

#res01.vad2.l.144 <- results(dds.vad2.l.144, alpha=0.01)
#summary(res01.vad2.l.144)
#sum(res01.vad2.l.144$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.144 = subset(res05.vad2.l.144, padj<0.05)
print(resSig05.vad2.l.144)
write.table(resSig05.vad2.l.144,file="resSig05.vad2.l.144.txt")
write.table(res.vad2.l.144,file="res.vad2.l.144.txt")

# Volcano plot
vp12<-EnhancedVolcano(res.vad2.l.144,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 144 hpi") + theme(plot.title=element_text(size=20)) 
vp12

tiff(file="ulmi.sci.data.figure.3.volcano.plot.tiff",width=16,height=24,units="in",res=300)
grid.arrange(vp01, vp05, vp09, vp02, vp06, vp10, vp03, vp07, vp11, vp04, vp08, vp12,
             nrow = 4)
dev.off()

# BARPLOTS

degs.bp <- data.frame(X = rep(c("MDV1 006hpi","MDV1 024hpi","MDV1 072hpi","MDV1 144hpi",
                                "MDV2.3 006hpi","MDV2.3 024hpi","MDV2.3 072hpi","MDV2.3 144hpi",
                                "VAD2 006hpi","VAD2 024hpi","VAD2 072hpi","VAD2 144hpi"),each=2),
                      Genotype = rep(c("MDV1","MDV2.3","VAD2"),each=8),
                      Time= rep(c("6 hpi","24 hpi","72 hpi", "144 hpi"),each=2,times=3),
                      Expression = rep(c("Induced","Repressed"),times=12),
                      DEGs = c(524,-270,935,-364,789,-210,656,-337,
                               370,-247,577,-716,656,-550,559,-205,
                               636,-759,573,-383,770,-355,714,-477))

tiff(file="ulmi.sci.data.figure.4.barplot.tiff",width=8,height=8,units="in",res=300)

ggplot(degs.bp, aes(x=fct_inorder(Time), y=DEGs,fill=Expression)) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  labs(x = "Time (hours post inoculation, hpi)", y = "Number of DEGs", fill = "") +
  #  theme(axis.text.x = element_text(angle = 45))
  facet_wrap(~Genotype, nrow = 1) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=DEGs), vjust=-0.3, size=4)+
  scale_fill_manual(values=c("forestgreen","firebrick1"))+
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size = 15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
dev.off()

### Venn diagrams and upset plot

list.mdv1.6<-    as.matrix(read.table("LISTS/list.mdv1.6.txt"))
list.mdv1.24<-   as.matrix(read.table("LISTS/list.mdv1.24.txt"))
list.mdv1.72<-   as.matrix(read.table("LISTS/list.mdv1.72.txt"))
list.mdv1.144<-  as.matrix(read.table("LISTS/list.mdv1.144.txt"))
list.mdv2.3.6<-  as.matrix(read.table("LISTS/list.mdv2.3.6.txt"))
list.mdv2.3.24<- as.matrix(read.table("LISTS/list.mdv2.3.24.txt"))
list.mdv2.3.72<- as.matrix(read.table("LISTS/list.mdv2.3.72.txt"))
list.mdv2.3.144<-as.matrix(read.table("LISTS/list.mdv2.3.144.txt"))
list.vad2.6<-    as.matrix(read.table("LISTS/list.vad2.6.txt"))
list.vad2.24<-   as.matrix(read.table("LISTS/list.vad2.24.txt"))
list.vad2.72<-   as.matrix(read.table("LISTS/list.vad2.72.txt"))
list.vad2.144<-  as.matrix(read.table("LISTS/list.vad2.144.txt"))

listInput <- list(MDV1_6hpi = list.mdv1.6, MDV1_24hpi = list.mdv1.24, MDV1_72hpi = list.mdv1.72, MDV1_144hpi = list.mdv1.144, 
                  MDV2.3_6hpi = list.mdv2.3.6,MDV2.3_24hpi = list.mdv2.3.24, MDV2.3_72hpi = list.mdv2.3.72,MDV2.3_144hpi = list.mdv2.3.144,
                  VAD2_6hpi = list.vad2.6,VAD2_24hpi = list.vad2.24,VAD2_72hpi = list.vad2.72,VAD2_144hpi = list.vad2.144)

tiff(file="ulmi.sci.data.figure.5.upset.plot.tiff",width=16,height=8,units="in",res=300)
upset(fromList(listInput), 
      nsets = 12,
      nintersects = 40, 
      sets=NULL,
      order.by = "freq", 
      #order.by = "degree",
      #empty.intersections = "on",
      mb.ratio= c(0.6, 0.4), 
      #group.by = "sets", 
      cutoff = 5,
      keep.order = T,
      number.angles = -30, 
      point.size = 3, line.size = 1.5, 
      #mainbar.y.label = "Intersection size", sets.x.label = "Set Size", 
      text.scale = c(2, # vertical barplot label
                     1.5, # vertical barplot axis
                     2, # horizontal barplot label
                     1.5, # horizontal barplot axis 
                     1.6, # set names 
                     1)) # bar numbers
dev.off()

# Save intersection lists
list_filter=listInput
str(list_filter)

#From this list, I generate two tables: One with the unique gene names

df2 <- data.frame(gene=unique(unlist(list_filter)))

head(df2)

dim(df2)


#One is simply a "dataframe" version of the list. 
#With every gene in the signature and the name of every signature (set).
library(dplyr)
df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

head(df1)
dim(df1)

#now I iterate the search of each unique gene name and save the identity of the signatures in a column.

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(V1==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

head(df_int,n=20)

dim(df_int)
write.table(df_int,file="upset_intersections.txt",sep="\t",row.names = FALSE, col.names = TRUE,quote=FALSE)

# lists
list1<-list(
  'MDV1 local - 6 hpi' = as.character(list.mdv1.6),
  'MDV1 local - 24 hpi' = as.character(list.mdv1.24),
  'MDV1 local - 72 hpi' = as.character(list.mdv1.72),
  'MDV1 local - 144 hpi' = as.character(list.mdv1.144))
list2<-list(
  'MDV2.3 local - 6 hpi' =   as.character(list.mdv2.3.6),
  'MDV2.3 local - 24 hpi' =  as.character(list.mdv2.3.24),
  'MDV2.3 local - 72 hpi' =  as.character(list.mdv2.3.72),
  'MDV2.3 local - 144 hpi' = as.character(list.mdv2.3.144))
list3<-list(
  'VAD2 local - 6 hpi' =   as.character(list.vad2.6),
  'VAD2 local - 24 hpi' =  as.character(list.vad2.24),
  'VAD2 local - 72 hpi' =  as.character(list.vad2.72),
  'VAD2 local - 144 hpi' = as.character(list.vad2.144))

# venns

venn1<-ggVennDiagram(list1[1:4], label_alpha = 0.9,
                     category.names = c("6 hpi (794 DEGs)","24 hpi (1299 DEGs)",
                                        "72 hpi (999 DEGs)","144 hpi (993 DEGs)"),
                     show_intersect = FALSE,set_size=3,label_size = 2.5,label_percent_digit = 1) + 
  scale_x_continuous(expand = expansion(mult=.13)) + scale_fill_distiller(palette = 9) + 
  theme_void() + #labs(title="(b)",element_text(size = 35)) +
  theme(text = element_text(size=rel(12)),legend.title=element_text(size=12), 
        legend.text=element_text(size=10),strip.text = element_text(size = 12))

tiff(file="ulmi.figure.3b.my_venn1.tiff",width=7,height=7,units="in",res=300)
venn1
dev.off()
venn2<-ggVennDiagram(list2[1:4], label_alpha = 0.9,
                     category.names = c("6 hpi (617 DEGs)","24 hpi (1293 DEGs)",
                                        "72 hpi (1206 DEGs)","144 hpi (764 DEGs)"),
                     show_intersect = FALSE,set_size=3,label_size = 2.5,label_percent_digit = 1) +
  scale_x_continuous(expand = expansion(mult=.13)) + scale_fill_distiller(palette = 15) + 
  theme_void() + #labs(title="(c)",element_text(size = 35)) +
  theme(text = element_text(size=rel(12)),legend.title=element_text(size=12), 
        legend.text=element_text(size=10),strip.text = element_text(size = 12))
tiff(file="ulmi.figure.3c.my_venn2.tiff",width=7,height=7,units="in",res=300)
venn2
dev.off()
venn3<-ggVennDiagram(list3[1:4], label_alpha = 0.9,
                     category.names = c("6 hpi (1395 DEGs)","24 hpi (956 DEGs)",
                                        "72 hpi (1125 DEGs)","144 hpi (1191 DEGs)"),
                     show_intersect = FALSE,set_size=3,label_size = 2.5,label_percent_digit = 1) +
  scale_x_continuous(expand = expansion(mult=.13)) + scale_fill_distiller(palette = 17) + 
  theme_void() + #labs(title="(d)",element_text(size = 35)) +
  theme(text = element_text(size=rel(12)),legend.title=element_text(size=12), 
        legend.text=element_text(size=10),strip.text = element_text(size = 12))
tiff(file="ulmi.figure.3d.my_venn3.tiff",width=7,height=7,units="in",res=300)
venn3
dev.off()
tiff(file="ulmi.figure.3a-d.tiff",width=15,height=8,units="in",res=300)
grid.arrange(barplot,venn1, venn2, venn3,nrow = 2)
dev.off()



################################################################################################
##########  5. A. Time-course of MDV1 LOCAL                              ##########

attach(local)
mdv1.l<-as.data.frame(cbind(MDV1.1_C_L_6, MDV1.3_C_L_6, MDV1.36_C_L_6, MDV1.47_C_L_6,
                            MDV1.4_C_L_24,MDV1.20_C_L_24,MDV1.22_C_L_24,MDV1.43_C_L_24,
                            MDV1.2_C_L_72,MDV1.11_C_L_72,MDV1.34_C_L_72,MDV1.42_C_L_72,
                            MDV1.8_C_L_144,MDV1.21_C_L_144,MDV1.37_C_L_144,MDV1.39_C_L_144,
                            MDV1.13_O_L_6, MDV1.15_O_L_6, MDV1.23_O_L_6, MDV1.50_O_L_6,
                            MDV1.16_O_L_24, MDV1.33_O_L_24,MDV1.41_O_L_24,MDV1.44_O_L_24,
                            MDV1.18_O_L_72,MDV1.31_O_L_72,MDV1.46_O_L_72,MDV1.49_O_L_72,
                            MDV1.9_O_L_144,MDV1.12_O_L_144,MDV1.24_O_L_144,MDV1.45_O_L_144))
detach()
rownames(mdv1.l)<-rownames(local)

coldata.mdv1.l<-c("MDV1.1","MDV1.3","MDV1.36","MDV1.47","MDV1.4","MDV1.20","MDV1.22","MDV1.43",
                  "MDV1.2","MDV1.11","MDV1.34","MDV1.42","MDV1.8","MDV1.21","MDV1.37","MDV1.39",
                  "MDV1.13","MDV1.15","MDV1.23","MDV1.50","MDV1.16","MDV1.33","MDV1.41","MDV1.44",
                  "MDV1.18","MDV1.31","MDV1.46","MDV1.49","MDV1.9","MDV1.12","MDV1.24","MDV1.45",
                  "control","control","control","control","control","control","control","control",
                  "control","control","control","control","control","control","control","control",
                  "infected","infected","infected","infected","infected","infected","infected","infected",
                  "infected","infected","infected","infected","infected","infected","infected","infected",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi",
                  "72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi",
                  "72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi",
                  "control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi",
                  "infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi",
                  "infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi")
coldata.mdv1.l<-matrix(coldata.mdv1.l,nrow=32,ncol=4,byrow=FALSE)
rownames(coldata.mdv1.l)<-colnames(mdv1.l)
colnames(coldata.mdv1.l)<-c("SAMPLE","TREATMENT","TIME","TREAT.TIME")

dds.mdv1.l <- DESeqDataSetFromMatrix(countData = mdv1.l,
                                     colData = coldata.mdv1.l,
                                     design = ~ TIME + TREATMENT + TIME:TREATMENT)

t_data.mdv1.l<-vst(dds.mdv1.l)
head(assay(t_data.mdv1.l))
distances.mdv1.l<-dist(t(assay(t_data.mdv1.l)))
distances_matrix.mdv1.l<-as.matrix(distances.mdv1.l)
rownames(distances_matrix.mdv1.l)<-paste(t_data.mdv1.l$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.mdv1.l<-hclust(distances.mdv1.l)
heatmap.2(distances_matrix.mdv1.l,Rowv = as.dendrogram(hc.mdv1.l),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE,
          key.title = "Color Key and Histogram",keysize = 1.5)

pcaData.mdv1.l <- plotPCA(t_data.mdv1.l, intgroup=c("TIME", "TREATMENT"), returnData=TRUE)
percentVar.mdv1.l <- round(100 * attr(pcaData.mdv1.l, "percentVar"))
pca.mdv1.l<-ggplot(pcaData.mdv1.l, aes(PC1, PC2, color=TREATMENT, shape=TIME)) +
  xlab(paste0("PC1: ",percentVar.mdv1.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.mdv1.l[2],"% variance")) + 
  geom_point(size=3) + # guides(color = FALSE, shape = FALSE) +
  theme_bw() + labs(title="(a)") + 
  xlim(-18, 18) + ylim(-18, 18) +
  theme(axis.title = element_text(size = 15)) + 
  theme(legend.title=element_text(size=13), 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=15),title = element_text(size = 15)) +
  coord_fixed()

dds.mdv1.l<-DESeq(dds.mdv1.l,test = "LRT", 
                  reduced = ~ TIME + TREATMENT)
res.mdv1.l<-results(dds.mdv1.l)
res.mdv1.l
resOrdered.mdv1.l <- res.mdv1.l[order(res.mdv1.l$pvalue),]
summary(res.mdv1.l)
sum(res.mdv1.l$padj < 0.1, na.rm=TRUE)

res05.mdv1.l <- results(dds.mdv1.l, alpha=0.05)
summary(res05.mdv1.l)
sum(res05.mdv1.l$padj < 0.05, na.rm=TRUE)

res01.mdv1.l <- results(dds.mdv1.l, alpha=0.01)
summary(res01.mdv1.l)
sum(res01.mdv1.l$padj < 0.01, na.rm=TRUE)

resSig05.mdv1.l = subset(res05.mdv1.l, padj<0.05)
print(resSig05.mdv1.l)
write.table(resSig05.mdv1.l,file="resSig05.mdv1.local.txt")



################################################################################################
##########  5. B. Time-course of MDV2.3 LOCAL                              ##########

attach(local)
mdv2.3.l<-as.data.frame(cbind(MDV2.3.8_C_L_6, MDV2.3.24_C_L_6, MDV2.3.29_C_L_6, MDV2.3.38_C_L_6,
                              MDV2.3.5_C_L_24,MDV2.3.6_C_L_24,MDV2.3.7_C_L_24,MDV2.3.43_C_L_24,
                              MDV2.3.3_C_L_72,MDV2.3.11_C_L_72,MDV2.3.41_C_L_72,MDV2.3.47_C_L_72,
                              MDV2.3.10_C_L_144,MDV2.3.27_C_L_144,MDV2.3.40_C_L_144,MDV2.3.42_C_L_144,
                              MDV2.3.18_O_L_6, MDV2.3.21_O_L_6, MDV2.3.33_O_L_6, MDV2.3.49_O_L_6,
                              MDV2.3.13_O_L_24, MDV2.3.17_O_L_24,MDV2.3.30_O_L_24,MDV2.3.45_O_L_24,
                              MDV2.3.15_O_L_72,MDV2.3.28_O_L_72,MDV2.3.32_O_L_72,MDV2.3.44_O_L_72,
                              MDV2.3.19_O_L_144,MDV2.3.26_O_L_144,MDV2.3.46_O_L_144,MDV2.3.50_O_L_144))
detach()
rownames(mdv2.3.l)<-rownames(local)
coldata.mdv2.3.l<-c("MDV2.3.8","MDV2.3.24","MDV2.3.29","MDV2.3.38","MDV2.3.5","MDV2.3.6","MDV2.3.7","MDV2.3.43",
                    "MDV2.3.3","MDV2.3.11","MDV2.3.41","MDV2.3.47","MDV2.3.10","MDV2.3.27","MDV2.3.40","MDV2.3.42",
                    "MDV2.3.18","MDV2.3.21","MDV2.3.33","MDV2.3.49","MDV2.3.13","MDV2.3.17","MDV2.3.30","MDV2.3.45",
                    "MDV2.3.15","MDV2.3.28","MDV2.3.32","MDV2.3.44","MDV2.3.19","MDV2.3.26","MDV2.3.46","MDV2.3.50",
                    "control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",
                    "infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected",
                    "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                    "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                    "control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi",
                    "control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi",
                    "infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi",
                    "infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi")

coldata.mdv2.3.l<-matrix(coldata.mdv2.3.l,nrow=32,ncol=4,byrow=FALSE)
rownames(coldata.mdv2.3.l)<-colnames(mdv2.3.l)
colnames(coldata.mdv2.3.l)<-c("SAMPLE","TREATMENT","TIME","TREAT.TIME")

dds.mdv2.3.l <- DESeqDataSetFromMatrix(countData = mdv2.3.l,
                                       colData = coldata.mdv2.3.l,
                                       design = ~ TIME + TREATMENT + TIME:TREATMENT)
t_data.mdv2.3.l<-vst(dds.mdv2.3.l)
head(assay(t_data.mdv2.3.l))
distances.mdv2.3.l<-dist(t(assay(t_data.mdv2.3.l)))
distances_matrix.mdv2.3.l<-as.matrix(distances.mdv2.3.l)
rownames(distances_matrix.mdv2.3.l)<-paste(t_data.mdv2.3.l$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.mdv2.3.l<-hclust(distances.mdv2.3.l)
heatmap.2(distances_matrix.mdv2.3.l,Rowv = as.dendrogram(hc.mdv2.3.l),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE)
pcaData.mdv2.3.l <- plotPCA(t_data.mdv2.3.l, intgroup=c("TIME", "TREATMENT"), returnData=TRUE)
percentVar.mdv2.3.l <- round(100 * attr(pcaData.mdv2.3.l, "percentVar"))
pca.madv2.3.l<-ggplot(pcaData.mdv2.3.l, aes(PC1, PC2, color=TREATMENT, shape=TIME)) +
  xlab(paste0("PC1: ",percentVar.mdv2.3.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.mdv2.3.l[2],"% variance")) + 
  geom_point(size=3) + # guides(color = FALSE, shape = FALSE) +
  theme_bw() + labs(title="(b)") + 
  xlim(-18, 18) + ylim(-18, 18) +
  theme(axis.title = element_text(size = 15)) + 
  theme(legend.title=element_text(size=13), 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=15),title = element_text(size = 15)) +
  coord_fixed()

dds.mdv2.3.l<-DESeq(dds.mdv2.3.l,test = "LRT", 
                    reduced = ~ TIME + TREATMENT)
res.mdv2.3.l<-results(dds.mdv2.3.l)
res.mdv2.3.l
resOrdered.mdv2.3.l <- res.mdv2.3.l[order(res.mdv2.3.l$pvalue),]
summary(res.mdv2.3.l)
sum(res.mdv2.3.l$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l <- results(dds.mdv2.3.l, alpha=0.05)
summary(res05.mdv2.3.l)
sum(res05.mdv2.3.l$padj < 0.05, na.rm=TRUE)

res01.mdv2.3.l <- results(dds.mdv2.3.l, alpha=0.01)
summary(res01.mdv2.3.l)
sum(res01.mdv2.3.l$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l = subset(res05.mdv2.3.l, padj<0.05)
print(resSig05.mdv2.3.l)
write.table(resSig05.mdv2.3.l,file="resSig05.mdv2.3.local.txt")





################################################################################################
##########  5. C. Time-course of VAD2 LOCAL                             ##########

attach(local)
vad2.l<-as.data.frame(cbind(VAD2.4_C_L_6, VAD2.6_C_L_6, VAD2.7_C_L_6, VAD2.37_C_L_6,
                            VAD2.1_C_L_24,VAD2.11_C_L_24,VAD2.19_C_L_24,VAD2.42_C_L_24,
                            VAD2.2_C_L_72,VAD2.9_C_L_72,VAD2.27_C_L_72,VAD2.41_C_L_72,
                            VAD2.8_C_L_144,VAD2.36_C_L_144,VAD2.40_C_L_144,VAD2.43_C_L_144,
                            VAD2.14_O_L_6, VAD2.15_O_L_6, VAD2.29_O_L_6, VAD2.44_O_L_6,
                            VAD2.25_O_L_24, VAD2.26_O_L_24,VAD2.28_O_L_24,VAD2.45_O_L_24,
                            VAD2.17_O_L_72,VAD2.22_O_L_72,VAD2.33_O_L_72,VAD2.48_O_L_72,
                            VAD2.18_O_L_144,VAD2.23_O_L_144,VAD2.35_O_L_144,VAD2.47_O_L_144))
detach()
rownames(vad2.l)<-rownames(local)
coldata.vad2.l<-c("VAD2.4","VAD2.6","VAD2.7","VAD2.37","VAD2.1","VAD2.11","VAD2.19","VAD2.42",
                  "VAD2.2","VAD2.9","VAD2.27","VAD2.41","VAD2.8","VAD2.36","VAD2.40","VAD2.43",
                  "VAD2.14","VAD2.15","VAD2.29","VAD2.44","VAD2.25","VAD2.26","VAD2.28","VAD2.45",
                  "VAD2.17","VAD2.22","VAD2.33","VAD2.48","VAD2.18","VAD2.23","VAD2.35","VAD2.47",
                  "control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",
                  "infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi",
                  "control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi",
                  "infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi",
                  "infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi")

coldata.vad2.l<-matrix(coldata.vad2.l,nrow=32,ncol=4,byrow=FALSE)
rownames(coldata.vad2.l)<-colnames(vad2.l)
colnames(coldata.vad2.l)<-c("SAMPLE","TREATMENT","TIME","TREAT.TIME")

dds.vad2.l <- DESeqDataSetFromMatrix(countData = vad2.l,
                                     colData = coldata.vad2.l,
                                     design = ~ TIME + TREATMENT + TIME:TREATMENT)
t_data.vad2.l<-vst(dds.vad2.l)
head(assay(t_data.vad2.l))
distances.vad2.l<-dist(t(assay(t_data.vad2.l)))
distances_matrix.vad2.l<-as.matrix(distances.vad2.l)
rownames(distances_matrix.vad2.l)<-paste(t_data.vad2.l$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.vad2.l<-hclust(distances.vad2.l)
heatmap.2(distances_matrix.vad2.l,Rowv = as.dendrogram(hc.vad2.l),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE)
pcaData.vad2.l <- plotPCA(t_data.vad2.l, intgroup=c("TREATMENT", "TIME"), returnData=TRUE)
percentVar.vad2.l <- round(100 * attr(pcaData.vad2.l, "percentVar"))

pca.vad2.l<-ggplot(pcaData.vad2.l, aes(PC1, PC2, color=TREATMENT, shape=TIME)) +
  xlab(paste0("PC1: ",percentVar.vad2.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vad2.l[2],"% variance")) + 
  geom_point(size=3) + # guides(color = FALSE, shape = FALSE) +
  theme_bw() + labs(title="(c)") + 
  xlim(-18, 18) + ylim(-18, 18) +
  theme(axis.title = element_text(size = 15)) + 
  theme(legend.title=element_text(size=13), 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=15),title = element_text(size = 15)) +
  coord_fixed() #+ scale_fill_discrete(breaks=c('6 hpi', '24 hpi', '72 hpi', "144 hpi"))

dds.vad2.l<-DESeq(dds.vad2.l,test = "LRT",
                  reduced = ~ TIME + TREATMENT)
res.vad2.l<-results(dds.vad2.l)
res.vad2.l
resOrdered.vad2.l <- res.vad2.l[order(res.vad2.l$pvalue),]
summary(res.vad2.l)
sum(res.vad2.l$padj < 0.1, na.rm=TRUE)

res05.vad2.l <- results(dds.vad2.l, alpha=0.05)
summary(res05.vad2.l)
sum(res05.vad2.l$padj < 0.05, na.rm=TRUE)

res01.vad2.l <- results(dds.vad2.l, alpha=0.01)
summary(res01.vad2.l)
sum(res01.vad2.l$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l = subset(res05.vad2.l, padj<0.05)
print(resSig05.vad2.l)
write.table(resSig05.vad2.l,file="resSig05.vad2.local.txt",sep="\t",dec=".")


tiff(file="ulmi.sci.data.supp.figure.s1.pcas.tiff",width=16,height=6,units="in",res=300)
grid.arrange(pca.mdv1.l, pca.madv2.3.l, pca.vad2.l,nrow = 1)
dev.off()



##########  7.C. Clustering MDV1 DEGs and plotting heatmap ####
## DEGs are considered when pvalue < 0.05 for MDV1 genotype and LFC>|2| in at least one time

mdv1.degs<-read.delim("degs_mdv1_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
mdv1.degs<-as.matrix(mdv1.degs)
mdv1.degs.my_hclust_gene <- hclust(dist(mdv1.degs), method = "ward.D")
as.dendrogram(mdv1.degs.my_hclust_gene) %>% plot(horiz = TRUE)
mdv1.degs.my_gene_col<- cutree(mdv1.degs.my_hclust_gene,k=4)
mdv1.degs.my_gene_col<- data.frame(Cluster=mdv1.degs.my_gene_col)
table(mdv1.degs.my_gene_col)

mdv1.degs.my_sample_col<-data.frame(Time=c("6 hpi","24 hpi","72 hpi"))#,"144 hpi"))
#mdv1.degs.my_sample_col<-data.frame(Time=c(c("6 hpi","24 hpi","72 hpi","144 hpi")))#,c(1,1,1,1)))
mdv1.degs.my_gene_col$Cluster <- factor(mdv1.degs.my_gene_col$Cluster,levels=c("1","2","3","4")) #,"5","6")) 

my_color=list(Cluster=c("1"="violet","2"="cyan","3"="dodgerblue2","4"="coral"))#,"5"="pink","6"="yellow"))#, 
paletteLength <- 50
myColor <- colorRampPalette(c("firebrick1", "lemonchiffon", "forestgreen"))(paletteLength)

myBreaks.mdv1 <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
                   seq(10/paletteLength, 10, length.out=floor(paletteLength/2)))

p1<-pheatmap(mdv1.degs,show_rownames = FALSE,cluster_cols = F,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D",
             main="DEGs MDV1",
             annotation_row=mdv1.degs.my_gene_col,
             #         annotation_col=mdv1.degs.my_sample_col,
             cutree_rows = 4,
             annotation_colors = my_color,
             cellwidth = 30, 
             cex = 1,
             color = myColor, 
             breaks = myBreaks.mdv1,
             legend_breaks = c(-12,-6,0,6,12),
             annotation_legend = TRUE,
             angle_col = 45,
             annotation_names_row = F,annotation_names_col = F)
p1
mdv1.degs.my_gene_col
write.table(mdv1.degs.my_gene_col, file="mdv1.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")

##########  7.D. Clustering MDV2.3 DEGs and plotting heatmap ####
## DEGs are considered when pvalue < 0.05 for mdv2.3 genotype and LFC>|2| in at least one time

mdv2.3.degs<-read.delim("degs_mdv2.3_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
mdv2.3.degs<-as.matrix(mdv2.3.degs)
mdv2.3.degs.my_hclust_gene <- hclust(dist(mdv2.3.degs), method = "ward.D")
as.dendrogram(mdv2.3.degs.my_hclust_gene) %>% plot(horiz = TRUE)
mdv2.3.degs.my_gene_col<- cutree(mdv2.3.degs.my_hclust_gene,k=4)
mdv2.3.degs.my_gene_col<- data.frame(Cluster=mdv2.3.degs.my_gene_col)
table(mdv2.3.degs.my_gene_col)

mdv2.3.degs.my_sample_col<-data.frame(Time=c("6 hpi","24 hpi","72 hpi"))#,"144 hpi"))
mdv2.3.degs.my_gene_col$Cluster <- factor(mdv2.3.degs.my_gene_col$Cluster,levels=c("1","2","3","4")) #,"5","6")) 

my_color=list(Cluster=c("1"="violet","2"="cyan","3"="dodgerblue2","4"="coral"))#,"5"="pink","6"="yellow"))#, 
paletteLength <- 50
myColor <- colorRampPalette(c("firebrick1", "lemonchiffon", "forestgreen"))(paletteLength)

myBreaks.mdv2.3 <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
                     seq(10/paletteLength, 10, length.out=floor(paletteLength/2))) ## setting max and min in 10 and -10

p2<-pheatmap(mdv2.3.degs,show_rownames = FALSE,cluster_cols = F,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D",
             main="DEGs MDV2.3",
             annotation_row=mdv2.3.degs.my_gene_col,
             #         annotation_col=mdv1.degs.my_sample_col,
             cutree_rows = 4,
             annotation_colors = my_color,
             cellwidth = 30, 
             cex = 1,
             color = myColor, 
             breaks = myBreaks.mdv2.3,
             legend_breaks = c(-12,-6,0,6,12),
             annotation_legend = TRUE,
             angle_col = 45,
             annotation_names_row = F,annotation_names_col = F)
p2
mdv2.3.degs.my_gene_col
write.table(mdv2.3.degs.my_gene_col, file="mdv2.3.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")


##########  7.E. Clustering VAD2 DEGs and plotting heatmap ####
## DEGs are considered when pvalue < 0.05 for VAD2 genotype and LFC>|2| in at least one time

vad2.degs<-read.delim("degs_vad2_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
vad2.degs<-as.matrix(vad2.degs)
vad2.degs.my_hclust_gene <- hclust(dist(vad2.degs), method = "ward.D")
as.dendrogram(vad2.degs.my_hclust_gene) %>% plot(horiz = TRUE)
vad2.degs.my_gene_col<- cutree(vad2.degs.my_hclust_gene,k=4)
vad2.degs.my_gene_col<- data.frame(Cluster=vad2.degs.my_gene_col)
table(vad2.degs.my_gene_col)

vad2.degs.my_sample_col<-data.frame(Time=c("6 hpi","24 hpi","72 hpi"))#,"144 hpi"))
vad2.degs.my_gene_col$Cluster <- factor(vad2.degs.my_gene_col$Cluster,levels=c("1","2","3","4")) #,"5","6")) 

my_color=list(Cluster=c("1"="violet","2"="cyan","3"="dodgerblue2","4"="coral"))#,"5"="pink","6"="yellow"))#, 
paletteLength <- 50
myColor <- colorRampPalette(c("firebrick1", "lemonchiffon", "forestgreen"))(paletteLength)

myBreaks.vad2 <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
                   seq(10/paletteLength, 10, length.out=floor(paletteLength/2))) ## setting max and min in 10 and -10
p3<-pheatmap(vad2.degs,show_rownames = FALSE,cluster_cols = F,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D",
             main="DEGs VAD2",
             annotation_row=vad2.degs.my_gene_col,
             #         annotation_col=vad2.degs.my_sample_col,
             cutree_rows = 4,
             annotation_colors = my_color,
             cellwidth = 30, 
             cex = 1,
             color = myColor, 
             breaks = myBreaks.vad2,
             legend_breaks = c(-12,-6,0,6,12),
             annotation_legend = TRUE,
             angle_col = 45,
             annotation_names_row = F,annotation_names_col = F)
p3
vad2.degs.my_gene_col
write.table(vad2.degs.my_gene_col, file="vad2.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")



#tiff(file="ulmi.my_heatmap_uw2.tiff",width=15,height=7,units="in",res=300)
#p4
#dev.off()


##########  7.F. Plotting Heatmaps of 3 genotypes together in grid ####
## DEGs are considered when pvalue < 0.05 for each genotype and LFC>|2| in at least one time
## There are three heatmaps together, one per genotype, not one heatmap for all three genotypes together

plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2']]=p2[[4]]
plot_list[['p3']]=p3[[4]]
plot_list

tiff(file="ulmi.figure.4.my_heatmaps.tiff",width=15,height=7,units="in",res=300)
h3<-grid.arrange(grobs=plot_list, ncol=3)
dev.off()


# Figure 5 Venn's diagram TIME-COURSE

mdv1.l.tc.degs<-  read.table("6a.mdv1.tc.list",stringsAsFactors = FALSE)
mdv2.3.l.tc.degs<-read.table("6b.mdv2.3.tc.list",stringsAsFactors = FALSE)
vad2.l.tc.degs<-  read.table("6c.vad2.tc.list",stringsAsFactors = FALSE)

mdv1.l.tc.degs<-  as.matrix(mdv1.l.tc.degs)
mdv2.3.l.tc.degs<-as.matrix(mdv2.3.l.tc.degs)
vad2.l.tc.degs<-  as.matrix(vad2.l.tc.degs)

list4<-list(
  'MDV2.3 local - time course'=as.character(mdv2.3.l.tc.degs),
  'MDV1 local - time course' = as.character(mdv1.l.tc.degs),
  'VAD2 local - time course' = as.character(vad2.l.tc.degs))

#ggVennDiagram(list3[1:4])
#ggVennDiagram(list3)

#f6<-ff

venn4<-ggVennDiagram(list4[1:3], label_alpha = 0.9,
                     category.names = c("MDV2.3\n(879 DEGs)","MDV1\n(1539 DEGs)",
                                        "VAD2\n(1500 DEGs)"),
                     show_intersect = FALSE,set_size=5,
                     label_size = 5,label_percent_digit = 1) + 
  scale_x_continuous(expand = expansion(mult=.13)) +
  scale_fill_distiller(palette = 3) +# labs(title="Time-course",element_text(size = 35)) +
  theme_void() +
  theme(text = element_text(size=rel(15)),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))

tiff(file="ulmi.figure.5.my_venn_time-course.tiff",width=8,height=8,units="in",res=300)
venn4
dev.off()

##########  8. Plotting cluster graphs for the 3 genotypes ####
## DEGs for each cluster. Clusters are numbered from 1 to 6, and genotypes as A, B and C for MDV1, MDV2.3 and VAD2 respectively
## There are three heatmaps together, one per genotype, not one heatmap for all three genotypes together

#### Cluster A1 ####
a1<-read.delim("mdv1.cluster1.txt",header=F,dec=".",sep="\t")
rnames.a1<-a1[,1]
cl.a1<-data.frame(a1[,2:ncol(a1)])
rownames(cl.a1)<-rnames.a1
colnames(cl.a1)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A2 ####
a2<-read.delim("mdv1.cluster2.txt",header=F,dec=".",sep="\t")
rnames.a2<-a2[,1]
cl.a2<-data.frame(a2[,2:ncol(a2)])
rownames(cl.a2)<-rnames.a2
colnames(cl.a2)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A3 ####
a3<-read.delim("mdv1.cluster3.txt",header=F,dec=".",sep="\t")
rnames.a3<-a3[,1]
cl.a3<-data.frame(a3[,2:ncol(a3)])
rownames(cl.a3)<-rnames.a3
colnames(cl.a3)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A4 ####
a4<-read.delim("mdv1.cluster4.txt",header=F,dec=".",sep="\t")
rnames.a4<-a4[,1]
cl.a4<-data.frame(a4[,2:ncol(a4)])
rownames(cl.a4)<-rnames.a4
colnames(cl.a4)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A5 ####
a5<-read.delim("mdv1.cluster5.txt",header=F,dec=".",sep="\t")
rnames.a5<-a5[,1]
cl.a5<-data.frame(a5[,2:ncol(a5)])
rownames(cl.a5)<-rnames.a5
colnames(cl.a5)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A6 ####
a6<-read.delim("mdv1.cluster6.txt",header=F,dec=".",sep="\t")
rnames.a6<-a6[,1]
cl.a6<-data.frame(a6[,2:ncol(a6)])
rownames(cl.a6)<-rnames.a6
colnames(cl.a6)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B1 ####
b1<-read.delim("mdv2.3.cluster1.txt",header=F,dec=".",sep="\t")
rnames.b1<-b1[,1]
cl.b1<-data.frame(b1[,2:ncol(b1)])
rownames(cl.b1)<-rnames.b1
colnames(cl.b1)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B2 ####
b2<-read.delim("mdv2.3.cluster2.txt",header=F,dec=".",sep="\t")
rnames.b2<-b2[,1]
cl.b2<-data.frame(b2[,2:ncol(b2)])
rownames(cl.b2)<-rnames.b2
colnames(cl.b2)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B3 ####
b3<-read.delim("mdv2.3.cluster3.txt",header=F,dec=".",sep="\t")
rnames.b3<-b3[,1]
cl.b3<-data.frame(b3[,2:ncol(b3)])
rownames(cl.b3)<-rnames.b3
colnames(cl.b3)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B4 ####
b4<-read.delim("mdv2.3.cluster4.txt",header=F,dec=".",sep="\t")
rnames.b4<-b4[,1]
cl.b4<-data.frame(b4[,2:ncol(b4)])
rownames(cl.b4)<-rnames.b4
colnames(cl.b4)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B5 ####
b5<-read.delim("mdv2.3.cluster5.txt",header=F,dec=".",sep="\t")
rnames.b5<-b5[,1]
cl.b5<-data.frame(b5[,2:ncol(b5)])
rownames(cl.b5)<-rnames.b5
colnames(cl.b5)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B6 ####
b6<-read.delim("mdv2.3.cluster6.txt",header=F,dec=".",sep="\t")
rnames.b6<-b6[,1]
cl.b6<-data.frame(b6[,2:ncol(b6)])
rownames(cl.b6)<-rnames.b6
colnames(cl.b6)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C1 ####
c1<-read.delim("vad2.cluster1.txt",header=F,dec=".",sep="\t")
rnames.c1<-c1[,1]
cl.c1<-data.frame(c1[,2:ncol(c1)])
rownames(cl.c1)<-rnames.c1
colnames(cl.c1)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C2 ####
c2<-read.delim("vad2.cluster2.txt",header=F,dec=".",sep="\t")
rnames.c2<-c2[,1]
cl.c2<-data.frame(c2[,2:ncol(c2)])
rownames(cl.c2)<-rnames.c2
colnames(cl.c2)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C3 ####
c3<-read.delim("vad2.cluster3.txt",header=F,dec=".",sep="\t")
rnames.c3<-c3[,1]
cl.c3<-data.frame(c3[,2:ncol(c3)])
rownames(cl.c3)<-rnames.c3
colnames(cl.c3)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C4 ####
c4<-read.delim("vad2.cluster4.txt",header=F,dec=".",sep="\t")
rnames.c4<-c4[,1]
cl.c4<-data.frame(c4[,2:ncol(c4)])
rownames(cl.c4)<-rnames.c4
colnames(cl.c4)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C5 ####
c5<-read.delim("vad2.cluster5.txt",header=F,dec=".",sep="\t")
rnames.c5<-c5[,1]
cl.c5<-data.frame(c5[,2:ncol(c5)])
rownames(cl.c5)<-rnames.c5
colnames(cl.c5)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C6 ####
c6<-read.delim("vad2.cluster6.txt",header=F,dec=".",sep="\t")
rnames.c6<-c6[,1]
cl.c6<-data.frame(c6[,2:ncol(c6)])
rownames(cl.c6)<-rnames.c6
colnames(cl.c6)<-c("6 hpi","24 hpi","72 hpi","144 hpi")

#### Plot 18 clusters ####

par(mfrow = c(2, 3))
ga1<-matplot(t(cl.a1),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a1)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 1 (100 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a1)))
ga1<-matplot(t(cl.a1),type='l',ylab=bquote(paste('Log'['2']*' Fold Change')),col="grey",
             lwd=1,xaxt='n',cex.lab=2,cex.axis=2)+
  lines(rowMeans(t(cl.a1)),col="red",lwd=2.5)+
  title("Cluster A1 (62 DEGs)",adj  = 0,cex.main=2) +
  axis(1,1:4,labels=rownames(t(cl.a1)),cex.axis=2)

ga2<-matplot(t(cl.a2),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a2)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 2 (79 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a2)))

ga3<-matplot(t(cl.a3),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a3)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 3 (78 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a3)))

ga4<-matplot(t(cl.a4),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a4)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 4 (149 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a4)))

ga5<-matplot(t(cl.a5),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a5)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 5 (110 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a5)))

ga6<-matplot(t(cl.a6),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a6)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 6 (101 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a6)))

gb1<-matplot(t(cl.b1),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b1)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 1 (148 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b1)))

gb2<-matplot(t(cl.b2),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b2)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 2 (267 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b2)))

gb3<-matplot(t(cl.b3),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b3)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 3 (178 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b3)))

gb4<-matplot(t(cl.b4),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b4)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 4 (55 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b4)))

gb5<-matplot(t(cl.b5),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b5)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 5 (79 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b5)))

gb6<-matplot(t(cl.b6),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b6)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 6 (90 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b6)))

gc1<-matplot(t(cl.c1),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c1)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 1 (171 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c1)))

gc2<-matplot(t(cl.c2),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c2)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 2 (118 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c2)))

gc3<-matplot(t(cl.c3),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c3)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 3 (185 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c3)))

gc4<-matplot(t(cl.c4),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c4)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 4 (120 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c4)))

gc5<-matplot(t(cl.c5),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c5)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 5 (127 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c5)))

gc6<-matplot(t(cl.c6),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c6)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 6 (56 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c6)))



## Figure 8 SEA LOCAL

mdv1.sea<-read.table("mdv1.sea",header=T)
g.sea1<-ggplot(mdv1.sea, aes(fill=Series, y=Percent, x=GO)) + 
  geom_bar(position="dodge", stat="identity",width=.5) +
  theme_bw()+ labs(x="GO Term",y="% of sequences") +
  theme(axis.text.x = element_text(angle = 90, size = 10))+
  scale_fill_brewer('Series', palette='Set2')+
  facet_grid(~Category,scale="free",space="free_x") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size = 15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
g.sea1

mdv2.3.sea<-read.table("mdv2.3.sea",header=T)
g.sea2<-ggplot(mdv2.3.sea, aes(fill=Series, y=Percent, x=GO)) + 
  geom_bar(position="dodge", stat="identity",width=.5) +
  theme_bw()+ labs(x="GO Term",y="% of sequences") +
  theme(axis.text.x = element_text(angle = 90, size = 10))+
  scale_fill_brewer('Series', palette='Set2')+
  facet_grid(~Category,scale="free",space="free_x") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size = 15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
g.sea2

vad2.sea<-read.table("vad2.sea",header=T)
g.sea3<-ggplot(vad2.sea, aes(fill=Series, y=Percent, x=GO)) + 
  geom_bar(position="dodge", stat="identity",width=.5) +
  theme_bw()+ labs(x="GO Term",y="% of sequences") +
  theme(axis.text.x = element_text(angle = 90, size = 10))+
  scale_fill_brewer('Series', palette='Set2')+
  facet_grid(~Category,scale="free",space="free_x") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size = 15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
g.sea3

tiff(file="ulmi.figure.6a-c.my_barplos_sea.tiff",width=15,height=15,units="in",res=300)
grid.arrange(g.sea1,g.sea2,g.sea3,nrow = 3)
dev.off()


g.sea1
g.sea2
g.sea3






#setwd("/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL")

library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggthemes)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(ggVennDiagram)
library(pheatmap)
#library(ComplexHeatmap)
library(dendextend)
library(dplyr)
library(ggplotify)
library(patchwork)
library(cowplot)
library(gridExtra)
library(UpSetR)

######## LOAD DATA, GENERATE OBJECTS, SORT BY NAMES, COMBINE INTO NEW TABLES AND SAVE THEM

##### LOAD CONTROL        ####

ctr.ulmi<-read.table("/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL/ulmi.counts.ctr.01.txt",header=TRUE, sep="\t",dec=".")
attach(ctr.ulmi)
ctr.ulmi<-as.data.frame(cbind(est_counts,    est_counts.1,  est_counts.2,  est_counts.3,
                              est_counts.4,  est_counts.5,  est_counts.6,  est_counts.7,  est_counts.8,
                              est_counts.9,  est_counts.10, est_counts.11, est_counts.12, est_counts.13,
                              est_counts.14, est_counts.15, est_counts.16, est_counts.17, est_counts.18,
                              est_counts.19, est_counts.20, est_counts.21, est_counts.22, est_counts.23,
                              est_counts.24, est_counts.25, est_counts.26, est_counts.27, est_counts.28,
                              est_counts.29, est_counts.30, est_counts.31, est_counts.32, est_counts.33,
                              est_counts.34, est_counts.35, est_counts.36, est_counts.37, est_counts.38,
                              est_counts.39, est_counts.40, est_counts.41, est_counts.42, est_counts.43,
                              est_counts.44, est_counts.45, est_counts.46, est_counts.47, est_counts.48,
                              est_counts.49, est_counts.50, est_counts.51, est_counts.52, est_counts.53,
                              est_counts.54, est_counts.55, est_counts.56, est_counts.57, est_counts.58,
                              est_counts.59, est_counts.60, est_counts.61, est_counts.62, est_counts.63,
                              est_counts.64, est_counts.65, est_counts.66, est_counts.67, est_counts.68,
                              est_counts.69, est_counts.70, est_counts.71, est_counts.72, est_counts.73,
                              est_counts.74, est_counts.75, est_counts.76, est_counts.77, est_counts.78,
                              est_counts.79, est_counts.80, est_counts.81, est_counts.82, est_counts.83,
                              est_counts.84, est_counts.85, est_counts.86, est_counts.87, est_counts.88,
                              est_counts.89, est_counts.90, est_counts.91, est_counts.92, est_counts.93,
                              est_counts.94, est_counts.95),row.names = target_id)
detach(ctr.ulmi)
colnames(ctr.ulmi)<-c("VAD2.1_C_L_24" ,"MDV2.3.8_C_L_6" ,"MDV1.43_C_D_24" ,"MDV1.20_C_L_24",
                      "VAD2.43_C_L_144" ,"VAD2.37_C_D_6" ,"VAD2.4_C_D_6" ,"VAD2.42_C_D_24",
                      "MDV1.34_C_D_72" ,"VAD2.6_C_D_6" ,"VAD2.11_C_L_24" ,"MDV1.37_C_D_144",
                      "VAD2.9_C_L_72" ,"MDV2.3.43_C_D_24" ,"VAD2.2_C_D_72" ,"VAD2.36_C_D_144",
                      "VAD2.41_C_L_72" ,"MDV1.36_C_D_6" ,"VAD2.8_C_L_144" ,"MDV1.42_C_L_72",
                      "MDV2.3.42_C_L_144" ,"MDV1.4_C_L_24" ,"MDV2.3.5_C_L_24" ,"MDV1.43_C_L_24",
                      "MDV1.3_C_L_6" ,"MDV1.3_C_D_6" ,"MDV2.3.3_C_L_72" ,"MDV2.3.10_C_D_144",
                      "MDV1.1_C_D_6" ,"MDV2.3.40_C_D_144" ,"MDV2.3.47_C_D_72" ,"MDV2.3.7_C_D_24", 
                      "VAD2.7_C_D_6" ,"MDV2.3.24_C_L_6" ,"MDV1.8_C_L_144" ,"VAD2.41_C_D_72",
                      "MDV1.47_C_D_6" ,"MDV1.21_C_L_144" ,"MDV2.3.3_C_D_72" ,"VAD2.8_C_D_144", 
                      "MDV2.3.47_C_L_72" ,"MDV1.39_C_D_144" ,"MDV1.20_C_D_24" ,"MDV1.11_C_L_72", 
                      "MDV2.3.8_C_D_6" ,"MDV1.47_C_L_6" ,"MDV2.3.27_C_L_144" ,"VAD2.40_C_L_144",
                      "VAD2.1_C_D_24" ,"MDV2.3.41_C_L_72" ,"MDV2.3.38_C_D_6" ,"MDV2.3.29_C_D_6", 
                      "MDV2.3.10_C_L_144" ,"VAD2.40_C_D_144" ,"MDV1.1_C_L_6" ,"MDV1.22_C_L_24",
                      "VAD2.11_C_D_24" ,"MDV2.3.6_C_L_24" ,"MDV2.3.11_C_L_72" ,"MDV1.8_C_D_144", 
                      "VAD2.7_C_L_6" ,"MDV2.3.38_C_L_6" ,"VAD2.27_C_L_72" ,"MDV1.2_C_D_72",
                      "MDV2.3.40_C_L_144" ,"VAD2.36_C_L_144" ,"VAD2.19_C_D_24" ,"VAD2.19_C_L_24",
                      "VAD2.4_C_L_6" ,"MDV1.2_C_L_72" ,"MDV2.3.5_C_D_24" ,"MDV1.42_C_D_72",
                      "MDV2.3.41_C_D_72" ,"MDV1.21_C_D_144" ,"MDV1.36_C_L_6" ,"VAD2.43_C_D_144",
                      "MDV2.3.29_C_L_6" ,"MDV1.22_C_D_24" ,"MDV1.39_C_L_144" ,"MDV1.37_C_L_144",
                      "MDV1.34_C_L_72" ,"MDV2.3.7_C_L_24" ,"VAD2.9_C_D_72" ,"VAD2.27_C_D_72",
                      "VAD2.2_C_L_72" ,"MDV2.3.11_C_D_72" ,"MDV1.11_C_D_72" ,"MDV2.3.6_C_D_24",
                      "MDV1.4_C_D_24" ,"MDV2.3.43_C_L_24" ,"VAD2.42_C_L_24" ,"MDV2.3.24_C_D_6",
                      "MDV2.3.27_C_D_144" ,"MDV2.3.42_C_D_144" ,"VAD2.6_C_L_6" ,"VAD2.37_C_L_6")
ctr.ulmi<-round(ctr.ulmi)


##### LOAD INFECTED        ####

exp.ulmi<-read.table("/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL/ulmi.counts.exp.01.txt",header=TRUE, sep="\t",dec=".")
attach(exp.ulmi)
exp.ulmi<-as.data.frame(cbind(est_counts,est_counts.1,est_counts.2,est_counts.3,
                              est_counts.4,est_counts.5,est_counts.6,est_counts.7,est_counts.8,
                              est_counts.9,est_counts.10,est_counts.11,est_counts.12,est_counts.13,
                              est_counts.14,est_counts.15,est_counts.16,est_counts.17,est_counts.18,
                              est_counts.19,est_counts.20,est_counts.21,est_counts.22,est_counts.23,
                              est_counts.24,est_counts.25,est_counts.26,est_counts.27,est_counts.28,
                              est_counts.29,est_counts.30,est_counts.31,est_counts.32,est_counts.33,
                              est_counts.34,est_counts.35,est_counts.36,est_counts.37,est_counts.38,
                              est_counts.39,est_counts.40,est_counts.41,est_counts.42,est_counts.43,
                              est_counts.44,est_counts.45,est_counts.46,est_counts.47,est_counts.48,
                              est_counts.49,est_counts.50,est_counts.51,est_counts.52,est_counts.53,
                              est_counts.54,est_counts.55,est_counts.56,est_counts.57,est_counts.58,
                              est_counts.59,est_counts.60,est_counts.61,est_counts.62,est_counts.63,
                              est_counts.64,est_counts.65,est_counts.66,est_counts.67,est_counts.68,
                              est_counts.69,est_counts.70,est_counts.71,est_counts.72,est_counts.73,
                              est_counts.74,est_counts.75,est_counts.76,est_counts.77,est_counts.78,
                              est_counts.79,est_counts.80,est_counts.81,est_counts.82,est_counts.83,
                              est_counts.84,est_counts.85,est_counts.86,est_counts.87,est_counts.88,
                              est_counts.89,est_counts.90,est_counts.91,est_counts.92,est_counts.93,
                              est_counts.94,est_counts.95),row.names = target_id)
detach(exp.ulmi)
colnames(exp.ulmi)<-c("VAD2.14_O_D_6","VAD2.29_O_L_6","MDV1.12_O_D_144","MDV1.50_O_D_6",
                      "VAD2.45_O_L_24","MDV1.23_O_D_6","MDV2.3.33_O_D_6","MDV1.13_O_D_6",
                      "VAD2.17_O_D_72","MDV1.16_O_L_24","MDV1.9_O_D_144","MDV2.3.26_O_L_144",
                      "VAD2.18_O_D_144","MDV2.3.32_O_L_72","MDV1.9_O_L_144","MDV2.3.18_O_L_6",
                      "VAD2.23_O_D_144","MDV2.3.32_O_D_72","VAD2.28_O_D_24","MDV1.18_O_D_72",
                      "MDV2.3.33_O_L_6","MDV1.41_O_L_24","MDV2.3.21_O_D_6","MDV1.12_O_L_144",
                      "MDV2.3.46_O_D_144","MDV2.3.17_O_D_24","MDV1.46_O_L_72","VAD2.33_O_D_72",
                      "VAD2.17_O_L_72","VAD2.44_O_D_6","VAD2.47_O_L_144","MDV1.15_O_D_6",
                      "MDV1.46_O_D_72","MDV1.31_O_D_72","MDV1.50_O_L_6","MDV2.3.30_O_D_24",
                      "MDV2.3.49_O_L_6","MDV2.3.17_O_L_24","MDV2.3.19_O_L_144","MDV2.3.50_O_D_144",
                      "MDV1.13_O_L_6","MDV2.3.45_O_L_24","MDV2.3.13_O_L_24","MDV2.3.28_O_D_72",
                      "MDV2.3.28_O_L_72","VAD2.45_O_D_24","MDV2.3.18_O_D_6","VAD2.25_O_D_24",
                      "VAD2.26_O_L_24","MDV1.33_O_L_24","MDV1.24_O_L_144","MDV1.45_O_D_144",
                      "MDV1.18_O_L_72","VAD2.18_O_L_144","VAD2.48_O_L_72","MDV1.49_O_D_72",
                      "MDV2.3.50_O_L_144","MDV1.15_O_L_6","VAD2.14_O_L_6","VAD2.35_O_L_144",
                      "MDV2.3.44_O_L_72","VAD2.35_O_D_144","VAD2.28_O_L_24","VAD2.23_O_L_144",
                      "VAD2.47_O_D_144","VAD2.15_O_L_6","MDV1.44_O_D_24","MDV1.16_O_D_24",
                      "VAD2.25_O_L_24","VAD2.15_O_D_6","MDV1.44_O_L_24","VAD2.29_O_D_6",
                      "VAD2.44_O_L_6","MDV2.3.15_O_D_72","VAD2.22_O_L_72","MDV2.3.30_O_L_24",
                      "VAD2.33_O_L_72","MDV1.31_O_L_72","MDV2.3.15_O_L_72","MDV2.3.44_O_D_72",
                      "MDV1.24_O_D_144","MDV2.3.49_O_D_6","MDV2.3.19_O_D_144","MDV1.23_O_L_6",
                      "MDV2.3.46_O_L_144","MDV1.33_O_D_24","VAD2.22_O_D_72","MDV2.3.26_O_D_144",
                      "MDV1.45_O_L_144","MDV1.41_O_D_24","VAD2.48_O_D_72","MDV2.3.21_O_L_6",
                      "MDV2.3.13_O_D_24","MDV2.3.45_O_D_24","MDV1.49_O_L_72","VAD2.26_O_D_24")
exp.ulmi<-round(exp.ulmi)


##### SORT BY NAME AND COMBINE CONTROL AND INFECTED IN ONE DATAMATRIX        ####

ctr.ordered <- ctr.ulmi[order(row.names(ctr.ulmi)),]
exp.ordered <- exp.ulmi[order(row.names(exp.ulmi)),]
full.dm<-cbind(ctr.ordered,exp.ordered)
rownames(full.dm)<-row.names(ctr.ordered)
#write.table(ctr.ordered,file="ulmi.counts.ctr.o.02.txt",sep="\t",dec =".",row.names = TRUE)
#write.table(exp.ordered,file="ulmi.counts.exp.o.02.txt",sep="\t",dec =".",row.names = TRUE)
#write.table(full.dm,file="ulmi.counts.full.dm.03.txt",sep="\t",dec =".",row.names = TRUE)
rm(ctr.ulmi)
rm(exp.ulmi)
rm(ctr.ordered)
rm(exp.ordered)


##### CREATE LOCAL DATAMATRIX        ####

detach()
attach(full.dm)
full.local <-as.data.frame(cbind(MDV1.1_C_L_6, MDV1.3_C_L_6, MDV1.36_C_L_6, MDV1.47_C_L_6,
                                 MDV1.4_C_L_24,MDV1.20_C_L_24,MDV1.22_C_L_24,MDV1.43_C_L_24,
                                 MDV1.2_C_L_72,MDV1.11_C_L_72,MDV1.34_C_L_72,MDV1.42_C_L_72,
                                 MDV1.8_C_L_144,MDV1.21_C_L_144,MDV1.37_C_L_144,MDV1.39_C_L_144,
                                 MDV1.13_O_L_6, MDV1.15_O_L_6, MDV1.23_O_L_6, MDV1.50_O_L_6,
                                 MDV1.16_O_L_24, MDV1.33_O_L_24,MDV1.41_O_L_24,MDV1.44_O_L_24,
                                 MDV1.18_O_L_72,MDV1.31_O_L_72,MDV1.46_O_L_72,MDV1.49_O_L_72,
                                 MDV1.9_O_L_144,MDV1.12_O_L_144,MDV1.24_O_L_144,MDV1.45_O_L_144,
                                 MDV2.3.8_C_L_6, MDV2.3.24_C_L_6, MDV2.3.29_C_L_6, MDV2.3.38_C_L_6,
                                 MDV2.3.5_C_L_24,MDV2.3.6_C_L_24,MDV2.3.7_C_L_24,MDV2.3.43_C_L_24,
                                 MDV2.3.3_C_L_72,MDV2.3.11_C_L_72,MDV2.3.41_C_L_72,MDV2.3.47_C_L_72,
                                 MDV2.3.10_C_L_144,MDV2.3.27_C_L_144,MDV2.3.40_C_L_144,MDV2.3.42_C_L_144,
                                 MDV2.3.18_O_L_6, MDV2.3.21_O_L_6, MDV2.3.33_O_L_6, MDV2.3.49_O_L_6,
                                 MDV2.3.13_O_L_24, MDV2.3.17_O_L_24,MDV2.3.30_O_L_24,MDV2.3.45_O_L_24,
                                 MDV2.3.15_O_L_72,MDV2.3.28_O_L_72,MDV2.3.32_O_L_72,MDV2.3.44_O_L_72,
                                 MDV2.3.19_O_L_144,MDV2.3.26_O_L_144,MDV2.3.46_O_L_144,MDV2.3.50_O_L_144,
                                 VAD2.4_C_L_6, VAD2.6_C_L_6, VAD2.7_C_L_6, VAD2.37_C_L_6,
                                 VAD2.1_C_L_24,VAD2.11_C_L_24,VAD2.19_C_L_24,VAD2.42_C_L_24,
                                 VAD2.2_C_L_72,VAD2.9_C_L_72,VAD2.27_C_L_72,VAD2.41_C_L_72,
                                 VAD2.8_C_L_144,VAD2.36_C_L_144,VAD2.40_C_L_144,VAD2.43_C_L_144,
                                 VAD2.14_O_L_6, VAD2.15_O_L_6, VAD2.29_O_L_6, VAD2.44_O_L_6,
                                 VAD2.25_O_L_24, VAD2.26_O_L_24,VAD2.28_O_L_24,VAD2.45_O_L_24,
                                 VAD2.17_O_L_72,VAD2.22_O_L_72,VAD2.33_O_L_72,VAD2.48_O_L_72,
                                 VAD2.18_O_L_144,VAD2.23_O_L_144,VAD2.35_O_L_144,VAD2.47_O_L_144))
rownames(full.local)<-rownames(full.dm)


###### CREATE DISTAL DATAMATRIX        ####
#
#full.distal <-as.data.frame(cbind(MDV1.1_C_D_6, MDV1.3_C_D_6, MDV1.36_C_D_6, MDV1.47_C_D_6,
#                                 MDV1.4_C_D_24,MDV1.20_C_D_24,MDV1.22_C_D_24,MDV1.43_C_D_24,
#                                 MDV1.2_C_D_72,MDV1.11_C_D_72,MDV1.34_C_D_72,MDV1.42_C_D_72,
#                                 MDV1.8_C_D_144,MDV1.21_C_D_144,MDV1.37_C_D_144,MDV1.39_C_D_144,
#                                 MDV1.13_O_D_6, MDV1.15_O_D_6, MDV1.23_O_D_6, MDV1.50_O_D_6,
#                                 MDV1.16_O_D_24, MDV1.33_O_D_24,MDV1.41_O_D_24,MDV1.44_O_D_24,
#                                 MDV1.18_O_D_72,MDV1.31_O_D_72,MDV1.46_O_D_72,MDV1.49_O_D_72,
#                                 MDV1.9_O_D_144,MDV1.12_O_D_144,MDV1.24_O_D_144,MDV1.45_O_D_144,
#                                 MDV2.3.8_C_D_6, MDV2.3.24_C_D_6, MDV2.3.29_C_D_6, MDV2.3.38_C_D_6,
#                                 MDV2.3.5_C_D_24,MDV2.3.6_C_D_24,MDV2.3.7_C_D_24,MDV2.3.43_C_D_24,
#                                 MDV2.3.3_C_D_72,MDV2.3.11_C_D_72,MDV2.3.41_C_D_72,MDV2.3.47_C_D_72,
#                                 MDV2.3.10_C_D_144,MDV2.3.27_C_D_144,MDV2.3.40_C_D_144,MDV2.3.42_C_D_144,
#                                 MDV2.3.18_O_D_6, MDV2.3.21_O_D_6, MDV2.3.33_O_D_6, MDV2.3.49_O_D_6,
#                                 MDV2.3.13_O_D_24, MDV2.3.17_O_D_24,MDV2.3.30_O_D_24,MDV2.3.45_O_D_24,
#                                 MDV2.3.15_O_D_72,MDV2.3.28_O_D_72,MDV2.3.32_O_D_72,MDV2.3.44_O_D_72,
#                                 MDV2.3.19_O_D_144,MDV2.3.26_O_D_144,MDV2.3.46_O_D_144,MDV2.3.50_O_D_144,
#                                 VAD2.4_C_D_6, VAD2.6_C_D_6, VAD2.7_C_D_6, VAD2.37_C_D_6,
#                                 VAD2.1_C_D_24,VAD2.11_C_D_24,VAD2.19_C_D_24,VAD2.42_C_D_24,
#                                 VAD2.2_C_D_72,VAD2.9_C_D_72,VAD2.27_C_D_72,VAD2.41_C_D_72,
#                                 VAD2.8_C_D_144,VAD2.36_C_D_144,VAD2.40_C_D_144,VAD2.43_C_D_144,
#                                 VAD2.14_O_D_6, VAD2.15_O_D_6, VAD2.29_O_D_6, VAD2.44_O_D_6,
#                                 VAD2.25_O_D_24, VAD2.26_O_D_24,VAD2.28_O_D_24,VAD2.45_O_D_24,
#                                 VAD2.17_O_D_72,VAD2.22_O_D_72,VAD2.33_O_D_72,VAD2.48_O_D_72,
#                                 VAD2.18_O_D_144,VAD2.23_O_D_144,VAD2.35_O_D_144,VAD2.47_O_D_144))
#rownames(full.distal)<-rownames(full.dm)
#
#rm(full.dm)
##write.table(full.local,file="ulmi.counts.local.dm.04.txt",sep="\t",dec =".",row.names = TRUE)
##write.table(full.distal,file="ulmi.counts.distal.dm.04.txt",sep="\t",dec =".",row.names = TRUE)
#rm(full.local)
#rm(full.distal)




#####

##### ##### LOAD LOCAL RESPONSE DATA FROM THE NEW TABLE         ####

local<-read.table("/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL/ulmi.counts.local.dm.04.txt",header=TRUE, sep="\t",dec=".")

# BY USING DESEQ2, WE CREATE THE ANALYSIS DESIGN, TRANSFORM THE DATA, AND PLOT PCA AND HEATMAP-DISTANCES          ####

coldata.local<-c(#samples 
  "MDV1.1","MDV1.3","MDV1.36","MDV1.47","MDV1.4","MDV1.20","MDV1.22","MDV1.43",
  "MDV1.2","MDV1.11","MDV1.34","MDV1.42","MDV1.8","MDV1.21","MDV1.37","MDV1.39",
  "MDV1.13","MDV1.15","MDV1.23","MDV1.50","MDV1.16","MDV1.33","MDV1.41","MDV1.44",
  "MDV1.18","MDV1.31","MDV1.46","MDV1.49","MDV1.9","MDV1.12","MDV1.24","MDV1.45",
  "MDV2.3.8","MDV2.3.24","MDV2.3.29","MDV2.3.38","MDV2.3.5","MDV2.3.6","MDV2.3.7","MDV2.3.43",
  "MDV2.3.3","MDV2.3.11","MDV2.3.41","MDV2.3.47","MDV2.3.10","MDV2.3.27","MDV2.3.40","MDV2.3.42",
  "MDV2.3.18","MDV2.3.21","MDV2.3.33","MDV2.3.49","MDV2.3.13","MDV2.3.17","MDV2.3.30","MDV2.3.45",
  "MDV2.3.15","MDV2.3.28","MDV2.3.32","MDV2.3.44","MDV2.3.19","MDV2.3.26","MDV2.3.46","MDV2.3.50",
  "VAD2.4","VAD2.6","VAD2.7","VAD2.37","VAD2.1","VAD2.11","VAD2.19","VAD2.42",
  "VAD2.2","VAD2.9","VAD2.27","VAD2.41","VAD2.8","VAD2.36","VAD2.40","VAD2.43",
  "VAD2.14","VAD2.15","VAD2.29","VAD2.44","VAD2.25","VAD2.26","VAD2.28","VAD2.45",
  "VAD2.17","VAD2.22","VAD2.33","VAD2.48","VAD2.18","VAD2.23","VAD2.35","VAD2.47",
  #genotypes
  "MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1",
  "MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1","MDV1",
  "MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3",
  "MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3","MDV2.3",
  "VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2",
  "VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2","VAD2",
  # condition per genotype and sample
  "control","control","control","control","control","control","control","control",
  "control","control","control","control","control","control","control","control",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "control","control","control","control","control","control","control","control",
  "control","control","control","control","control","control","control","control",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "control","control","control","control","control","control","control","control",
  "control","control","control","control","control","control","control","control",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  "inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated","inoculated",
  # time per condition, genotype and sample
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  "006","006","006","006","024","024","024","024","072","072","072","072","144","144","144","144",
  # GENOTYPE AND TREATMENT
  "MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR",
  "MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR","MDV1_CTR",
  "MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO",
  "MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO","MDV1_INO",
  "MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR",
  "MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR","MDV2.3_CTR",
  "MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO",
  "MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO","MDV2.3_INO",
  "VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR",
  "VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR","VAD2_CTR",
  "VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO",
  "VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO","VAD2_INO",
  # GENOTYPE AND TIME
  "MDV1_6","MDV1_6","MDV1_6","MDV1_6","MDV1_24","MDV1_24","MDV1_24","MDV1_24",
  "MDV1_72","MDV1_72","MDV1_72","MDV1_72","MDV1_144","MDV1_144","MDV1_144","MDV1_144",
  "MDV1_6","MDV1_6","MDV1_6","MDV1_6","MDV1_24","MDV1_24","MDV1_24","MDV1_24",
  "MDV1_72","MDV1_72","MDV1_72","MDV1_72","MDV1_144","MDV1_144","MDV1_144","MDV1_144",
  "MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_24","MDV2.3_24","MDV2.3_24","MDV2.3_24",
  "MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_144","MDV2.3_144","MDV2.3_144","MDV2.3_144",
  "MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_6","MDV2.3_24","MDV2.3_24","MDV2.3_24","MDV2.3_24",
  "MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_72","MDV2.3_144","MDV2.3_144","MDV2.3_144","MDV2.3_144",
  "VAD2_6","VAD2_6","VAD2_6","VAD2_6","VAD2_24","VAD2_24","VAD2_24","VAD2_24",
  "VAD2_72","VAD2_72","VAD2_72","VAD2_72","VAD2_144","VAD2_144","VAD2_144","VAD2_144",
  "VAD2_6","VAD2_6","VAD2_6","VAD2_6","VAD2_24","VAD2_24","VAD2_24","VAD2_24",
  "VAD2_72","VAD2_72","VAD2_72","VAD2_72","VAD2_144","VAD2_144","VAD2_144","VAD2_144",
  # TREATMENT AND TIME
  "CTR_6","CTR_6","CTR_6","CTR_6","CTR_24","CTR_24","CTR_24","CTR_24",
  "CTR_72","CTR_72","CTR_72","CTR_72","CTR_144","CTR_144","CTR_144","CTR_144",
  "INO_6","INO_6","INO_6","INO_6","INO_24","INO_24","INO_24","INO_24",
  "INO_72","INO_72","INO_72","INO_72","INO_144","INO_144","INO_144","INO_144",
  "CTR_6","CTR_6","CTR_6","CTR_6","CTR_24","CTR_24","CTR_24","CTR_24",
  "CTR_72","CTR_72","CTR_72","CTR_72","CTR_144","CTR_144","CTR_144","CTR_144",
  "INO_6","INO_6","INO_6","INO_6","INO_24","INO_24","INO_24","INO_24",
  "INO_72","INO_72","INO_72","INO_72","INO_144","INO_144","INO_144","INO_144",
  "CTR_6","CTR_6","CTR_6","CTR_6","CTR_24","CTR_24","CTR_24","CTR_24",
  "CTR_72","CTR_72","CTR_72","CTR_72","CTR_144","CTR_144","CTR_144","CTR_144",
  "INO_6","INO_6","INO_6","INO_6","INO_24","INO_24","INO_24","INO_24",
  "INO_72","INO_72","INO_72","INO_72","INO_144","INO_144","INO_144","INO_144")

coldata.local<-matrix(coldata.local,nrow=96,ncol=7,byrow=FALSE)

rownames(coldata.local)<-c("MDV1.1_C_L_6","MDV1.3_C_L_6","MDV1.36_C_L_6","MDV1.47_C_L_6",
                           "MDV1.4_C_L_24","MDV1.20_C_L_24","MDV1.22_C_L_24","MDV1.43_C_L_24",
                           "MDV1.2_C_L_72","MDV1.11_C_L_72","MDV1.34_C_L_72","MDV1.42_C_L_72",
                           "MDV1.8_C_L_144","MDV1.21_C_L_144","MDV1.37_C_L_144","MDV1.39_C_L_144",
                           "MDV1.13_O_L_6","MDV1.15_O_L_6","MDV1.23_O_L_6","MDV1.50_O_L_6",
                           "MDV1.16_O_L_24","MDV1.33_O_L_24","MDV1.41_O_L_24","MDV1.44_O_L_24",
                           "MDV1.18_O_L_72","MDV1.31_O_L_72","MDV1.46_O_L_72","MDV1.49_O_L_72",
                           "MDV1.9_O_L_144","MDV1.12_O_L_144","MDV1.24_O_L_144","MDV1.45_O_L_144",
                           "MDV2.3.8_C_L_6","MDV2.3.24_C_L_6","MDV2.3.29_C_L_6","MDV2.3.38_C_L_6",
                           "MDV2.3.5_C_L_24","MDV2.3.6_C_L_24","MDV2.3.7_C_L_24","MDV2.3.43_C_L_24",
                           "MDV2.3.3_C_L_72","MDV2.3.11_C_L_72","MDV2.3.41_C_L_72","MDV2.3.47_C_L_72",
                           "MDV2.3.10_C_L_144","MDV2.3.27_C_L_144","MDV2.3.40_C_L_144","MDV2.3.42_C_L_144",
                           "MDV2.3.18_O_L_6","MDV2.3.21_O_L_6","MDV2.3.33_O_L_6","MDV2.3.49_O_L_6",
                           "MDV2.3.13_O_L_24","MDV2.3.17_O_L_24","MDV2.3.30_O_L_24","MDV2.3.45_O_L_24",
                           "MDV2.3.15_O_L_72","MDV2.3.28_O_L_72","MDV2.3.32_O_L_72","MDV2.3.44_O_L_72",
                           "MDV2.3.19_O_L_144","MDV2.3.26_O_L_144","MDV2.3.46_O_L_144","MDV2.3.50_O_L_144",
                           "VAD2.4_C_L_6","VAD2.6_C_L_6","VAD2.7_C_L_6","VAD2.37_C_L_6",
                           "VAD2.1_C_L_24","VAD2.11_C_L_24","VAD2.19_C_L_24","VAD2.42_C_L_24",
                           "VAD2.2_C_L_72","VAD2.9_C_L_72","VAD2.27_C_L_72","VAD2.41_C_L_72",
                           "VAD2.8_C_L_144","VAD2.36_C_L_144","VAD2.40_C_L_144","VAD2.43_C_L_144",
                           "VAD2.14_O_L_6","VAD2.15_O_L_6","VAD2.29_O_L_6","VAD2.44_O_L_6",
                           "VAD2.25_O_L_24","VAD2.26_O_L_24","VAD2.28_O_L_24","VAD2.45_O_L_24",
                           "VAD2.17_O_L_72","VAD2.22_O_L_72","VAD2.33_O_L_72","VAD2.48_O_L_72",
                           "VAD2.18_O_L_144","VAD2.23_O_L_144","VAD2.35_O_L_144","VAD2.47_O_L_144")

colnames(coldata.local)<-c("SAMPLE","GENOTYPE","TREATMENT","TIME","GENOT.TREATMENT","GENOT.TIME","TREAT.TIME")

dds.local <- DESeqDataSetFromMatrix(countData = local,
                                     colData = coldata.local,
                                     design = ~ GENOTYPE + TIME + TREATMENT + TIME:TREATMENT)

## DATA VISUALIZATION: TRANSFORMATION, PCAs, HEATMAPS...
t_data.local<-vst(dds.local)
#t_data.local <- assay(vst(dds.local, blind=FALSE)) # create a table with vst values and save it for filtering later DEGs and make heatmap
#write.table(as.data.frame   (t_data.local), file="ulmi.t_data.vst.local.txt",    sep="\t",dec =".",row.names = TRUE)

head(assay(t_data.local))
distances.local<-dist(t(assay(t_data.local)))
distances_matrix.local<-as.matrix(distances.local)
rownames(distances_matrix.local)<-paste(t_data.local$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.local<-hclust(distances.local)

heatmap.2(distances_matrix.local,Rowv = as.dendrogram(hc.local),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE)

#plotPCA(t_data.local,intgroup="SAMPLE")
plotPCA(t_data.local,intgroup="GENOTYPE")
plotPCA(t_data.local,intgroup="TREATMENT")
plotPCA(t_data.local,intgroup="TIME")
plotPCA(t_data.local,intgroup="GENOT.TREATMENT")
plotPCA(t_data.local,intgroup="GENOT.TIME")
plotPCA(t_data.local,intgroup="TREAT.TIME")

pcaData.l <- plotPCA(t_data.local, intgroup=c("TREATMENT", "GENOTYPE"), returnData=TRUE)
percentVar.l <- round(100 * attr(pcaData.l, "percentVar"))

tiff(file="ulmi.my_pca_plot.tiff",width=8,height=8,units="in",res=300)
ggplot(pcaData.l, aes(PC1, PC2, shape=TREATMENT, color=GENOTYPE)) +
  xlab(paste0("PC1: ",percentVar.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.l[2],"% variance")) + 
  geom_point(size=3.5) +# guides(color = FALSE, shape = FALSE) +
  theme_bw() + #labs(title="(a)") + 
  xlim(-15, 15) + ylim(-15, 15) +
  theme(axis.title = element_text(size = 15)) + 
  theme(legend.title=element_text(size=15), 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=15)) +
  coord_fixed()
dev.off()

#####

##### ##### SEPARATE GENOTYPES AND TIMES AND PERFORM PAIRWISE ANALYSIS BY COMPARING CONTROL VS INFECTED PLANTS         ####
#####
##### ##### MDV1  6 hpi       ####

mdv1.l.6<-cbind(local$MDV1.1_C_L_6, local$MDV1.3_C_L_6,  local$MDV1.36_C_L_6, local$MDV1.47_C_L_6,
                local$MDV1.13_O_L_6,local$MDV1.15_O_L_6, local$MDV1.23_O_L_6, local$MDV1.50_O_L_6)

colnames(mdv1.l.6)<-c("MDV1.1_C_L_6","MDV1.3_C_L_6","MDV1.36_C_L_6","MDV1.47_C_L_6",
                      "MDV1.13_O_L_6","MDV1.15_O_L_6","MDV1.23_O_L_6","MDV1.50_O_L_6")
rownames(mdv1.l.6)<-rownames(local)

coldata.mdv1.l.6<-c("MDV1.1",  "MDV1.3",  "MDV1.36",  "MDV1.47","MDV1.13","MDV1.15","MDV1.23","MDV1.50",
                    "control","control","control","control","infected","infected","infected","infected")
coldata.mdv1.l.6<-matrix(coldata.mdv1.l.6,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv1.l.6)<-colnames(mdv1.l.6)
colnames(coldata.mdv1.l.6)<-c("SAMPLE","TREATMENT")
print(coldata.mdv1.l.6)
dds.mdv1.l.6 <- DESeqDataSetFromMatrix(countData = mdv1.l.6,
                                       colData = coldata.mdv1.l.6,
                                       design = ~ TREATMENT)
dds.mdv1.l.6<-DESeq(dds.mdv1.l.6)
res.mdv1.l.6<-results(dds.mdv1.l.6)
res.mdv1.l.6
resOrdered.mdv1.l.6 <- res.mdv1.l.6[order(res.mdv1.l.6$pvalue),]
summary(res.mdv1.l.6)
sum(res.mdv1.l.6$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.6 <- results(dds.mdv1.l.6, alpha=0.05)
summary(res05.mdv1.l.6)
sum(res05.mdv1.l.6$padj < 0.05, na.rm=TRUE)

res01.mdv1.l.6 <- results(dds.mdv1.l.6, alpha=0.01)
summary(res01.mdv1.l.6)
sum(res01.mdv1.l.6$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.6 = subset(res05.mdv1.l.6, padj<0.05)
print(resSig05.mdv1.l.6)
write.table(resSig05.mdv1.l.6,file="resSig05.mdv1.l.6.txt")
write.table(res.mdv1.l.6,file="res.mdv1.l.6.txt")

# Volcano plot
vp01<-EnhancedVolcano(res.mdv1.l.6,
                     lab = NA,
                     x = 'log2FoldChange',
                     y = 'padj',
                     #title="(a) FAR3 DEGs",
                     subtitle = NULL,
                     #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                     caption = bquote(" "),
                     pCutoff = 0.05,
                     FCcutoff = 2,
                     pointSize=1,col=c('black', 'black', 'black', 'red3'),
                     boxedLabels = F,
                     legendPosition = "none") +  
#  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 6 hpi") + theme(plot.title=element_text(size=20)) 
vp01
#ggsave("Fig3a.test.tiff", width = 6, height = 6, device='tiff', dpi=300)

##### ##### MDV1  24 hpi       ####

mdv1.l.24<-cbind(local$MDV1.4_C_L_24,  local$MDV1.20_C_L_24, local$MDV1.22_C_L_24, local$MDV1.43_C_L_24,
                 local$MDV1.16_O_L_24, local$MDV1.33_O_L_24, local$MDV1.41_O_L_24, local$MDV1.44_O_L_24)
colnames(mdv1.l.24)<-c("MDV1.4_C_L_24",  "MDV1.20_C_L_24",  "MDV1.22_C_L_24",  "MDV1.43_C_L_24",
                       "MDV1.16_O_L_24", "MDV1.33_O_L_24", "MDV1.41_O_L_24", "MDV1.44_O_L_24")
rownames(mdv1.l.24)<-rownames(local)
coldata.mdv1.l.24<-c("MDV1.4",  "MDV1.20",  "MDV1.22",  "MDV1.43","MDV1.16","MDV1.33","MDV1.41","MDV1.44",
                     "control","control","control","control","infected","infected","infected","infected")
coldata.mdv1.l.24<-matrix(coldata.mdv1.l.24,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv1.l.24)<-colnames(mdv1.l.24)
colnames(coldata.mdv1.l.24)<-c("SAMPLE","TREATMENT")
print(coldata.mdv1.l.24)

dds.mdv1.l.24 <- DESeqDataSetFromMatrix(countData = mdv1.l.24,
                                       colData = coldata.mdv1.l.24,
                                       design = ~ TREATMENT)
dds.mdv1.l.24<-DESeq(dds.mdv1.l.24)
res.mdv1.l.24<-results(dds.mdv1.l.24)
res.mdv1.l.24
resOrdered.mdv1.l.24 <- res.mdv1.l.24[order(res.mdv1.l.24$pvalue),]
summary(res.mdv1.l.24)
sum(res.mdv1.l.24$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.24 <- results(dds.mdv1.l.24, alpha=0.05)
summary(res05.mdv1.l.24)
sum(res05.mdv1.l.24$padj < 0.05, na.rm=TRUE)

res01.mdv1.l.24 <- results(dds.mdv1.l.24, alpha=0.01)
summary(res01.mdv1.l.24)
sum(res01.mdv1.l.24$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.24 = subset(res05.mdv1.l.24, padj<0.05)
print(resSig05.mdv1.l.24)
write.table(resSig05.mdv1.l.24,file="resSig05.mdv1.l.24.txt")
write.table(res.mdv1.l.24,file="res.mdv1.l.24.txt")

# Volcano plot
vp02<-EnhancedVolcano(res.mdv1.l.24,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 24 hpi") + theme(plot.title=element_text(size=20)) 
vp02

##### ##### MDV1  72 hpi       ####

mdv1.l.72<-cbind(local$MDV1.2_C_L_72,  local$MDV1.11_C_L_72,local$MDV1.34_C_L_72,local$MDV1.42_C_L_72,
                 local$MDV1.18_O_L_72, local$MDV1.31_O_L_72,local$MDV1.46_O_L_72,local$MDV1.49_O_L_72)

colnames(mdv1.l.72)<-c("MDV1.2_C_L_72",  "MDV1.11_C_L_72",  "MDV1.34_C_L_72",  "MDV1.42_C_L_72",
                       "MDV1.18_O_L_72", "MDV1.31_O_L_72", "MDV1.46_O_L_72", "MDV1.49_O_L_72")

rownames(mdv1.l.72)<-rownames(local)

coldata.mdv1.l.72<-c("MDV1.2",  "MDV1.11",  "MDV1.34",  "MDV1.42","MDV1.18","MDV1.31","MDV1.46","MDV1.49",
                     "control","control","control","control","infected","infected","infected","infected")

coldata.mdv1.l.72<-matrix(coldata.mdv1.l.72,nrow=8,ncol=2,byrow = FALSE)

rownames(coldata.mdv1.l.72)<-colnames(mdv1.l.72)

colnames(coldata.mdv1.l.72)<-c("SAMPLE","TREATMENT")

print(coldata.mdv1.l.72)

dds.mdv1.l.72 <- DESeqDataSetFromMatrix(countData = mdv1.l.72,
                                       colData = coldata.mdv1.l.72,
                                       design = ~ TREATMENT)
dds.mdv1.l.72<-DESeq(dds.mdv1.l.72)
res.mdv1.l.72<-results(dds.mdv1.l.72)
res.mdv1.l.72
resOrdered.mdv1.l.72 <- res.mdv1.l.72[order(res.mdv1.l.72$pvalue),]
summary(res.mdv1.l.72)
sum(res.mdv1.l.72$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.72 <- results(dds.mdv1.l.72, alpha=0.05)
summary(res05.mdv1.l.72)
sum(res05.mdv1.l.72$padj < 0.05, na.rm=TRUE)

res01.mdv1.l.72 <- results(dds.mdv1.l.72, alpha=0.01)
summary(res01.mdv1.l.72)
sum(res01.mdv1.l.72$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.72 = subset(res05.mdv1.l.72, padj<0.05)
print(resSig05.mdv1.l.72)
write.table(resSig05.mdv1.l.72,file="resSig05.mdv1.l.72.txt")
write.table(res.mdv1.l.72,file="res.mdv1.l.72.txt")

# Volcano plot
vp03<-EnhancedVolcano(res.mdv1.l.72,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 72 hpi") + theme(plot.title=element_text(size=20)) 
vp03

##### ##### MDV1  144 hpi       ####

mdv1.l.144<-cbind(local$MDV1.8_C_L_144,local$MDV1.21_C_L_144,local$MDV1.37_C_L_144,local$MDV1.39_C_L_144,
                  local$MDV1.9_O_L_144,local$MDV1.12_O_L_144,local$MDV1.24_O_L_144,local$MDV1.45_O_L_144)

colnames(mdv1.l.144)<-c("MDV1.8_C_L_144",  "MDV1.21_C_L_144",  "MDV1.37_C_L_144",  "MDV1.39_C_L_144",
                        "MDV1.9_O_L_144", "MDV1.12_O_L_144","MDV1.24_O_L_144","MDV1.45_O_L_144")

rownames(mdv1.l.144)<-rownames(local)

coldata.mdv1.l.144<-c("MDV1.8",  "MDV1.21",  "MDV1.37",  "MDV1.39","MDV1.9","MDV1.12","MDV1.24","MDV1.45",
                      "control","control","control","control","infected","infected","infected","infected")

coldata.mdv1.l.144<-matrix(coldata.mdv1.l.144,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv1.l.144)<-colnames(mdv1.l.144)
colnames(coldata.mdv1.l.144)<-c("SAMPLE","TREATMENT")
print(coldata.mdv1.l.144)

dds.mdv1.l.144 <- DESeqDataSetFromMatrix(countData = mdv1.l.144,
                                         colData = coldata.mdv1.l.144,
                                         design = ~ TREATMENT)

dds.mdv1.l.144<-DESeq(dds.mdv1.l.144)
res.mdv1.l.144<-results(dds.mdv1.l.144)
res.mdv1.l.144
resOrdered.mdv1.l.144 <- res.mdv1.l.144[order(res.mdv1.l.144$pvalue),]
summary(res.mdv1.l.144)
sum(res.mdv1.l.144$padj < 0.1, na.rm=TRUE)

res05.mdv1.l.144 <- results(dds.mdv1.l.144, alpha=0.05)
summary(res05.mdv1.l.144)
sum(res05.mdv1.l.144$padj < 0.05, na.rm=TRUE)

res01.mdv1.l.144 <- results(dds.mdv1.l.144, alpha=0.01)
summary(res01.mdv1.l.144)
sum(res01.mdv1.l.144$padj < 0.01, na.rm=TRUE)

## save significant genes (all and pval<0.05)

resSig05.mdv1.l.144 = subset(res05.mdv1.l.144, padj<0.05)
print(resSig05.mdv1.l.144)
write.table(resSig05.mdv1.l.144,file="resSig05.mdv1.l.144.txt")
write.table(res.mdv1.l.144,file="res.mdv1.l.144.txt")

# Volcano plot
vp04<-EnhancedVolcano(res.mdv1.l.144,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV1 144 hpi") + theme(plot.title=element_text(size=20)) 
vp04


##### ##### MDV2.3  6 hpi       ####

mdv2.3.l.6<-cbind(local$MDV2.3.8_C_L_6, local$MDV2.3.24_C_L_6,local$MDV2.3.29_C_L_6,local$MDV2.3.38_C_L_6,
                  local$MDV2.3.18_O_L_6,local$MDV2.3.21_O_L_6,local$MDV2.3.33_O_L_6,local$MDV2.3.49_O_L_6)
colnames(mdv2.3.l.6)<-c("MDV2.3.8_C_L_6","MDV2.3.24_C_L_6","MDV2.3.29_C_L_6","MDV2.3.38_C_L_6",
                        "MDV2.3.18_O_L_6",  "MDV2.3.21_O_L_6",  "MDV2.3.33_O_L_6",  "MDV2.3.49_O_L_6")
rownames(mdv2.3.l.6)<-rownames(local)
coldata.mdv2.3.l.6<-c("MDV2.3.8","MDV2.3.24","MDV2.3.29","MDV2.3.38","MDV2.3.18","MDV2.3.21","MDV2.3.33","MDV2.3.49",
                      "control","control","control","control","infected","infected","infected","infected")
coldata.mdv2.3.l.6<-matrix(coldata.mdv2.3.l.6,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.6)<-colnames(mdv2.3.l.6)
colnames(coldata.mdv2.3.l.6)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.6)

dds.mdv2.3.l.6 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.6,
                                         colData = coldata.mdv2.3.l.6,
                                         design = ~ TREATMENT)
dds.mdv2.3.l.6<-DESeq(dds.mdv2.3.l.6)
res.mdv2.3.l.6<-results(dds.mdv2.3.l.6)
res.mdv2.3.l.6
resOrdered.mdv2.3.l.6 <- res.mdv2.3.l.6[order(res.mdv2.3.l.6$pvalue),]
summary(res.mdv2.3.l.6)
sum(res.mdv2.3.l.6$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.6 <- results(dds.mdv2.3.l.6, alpha=0.05)
summary(res05.mdv2.3.l.6)
sum(res05.mdv2.3.l.6$padj < 0.05, na.rm=TRUE)

res01.mdv2.3.l.6 <- results(dds.mdv2.3.l.6, alpha=0.01)
summary(res01.mdv2.3.l.6)
sum(res01.mdv2.3.l.6$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.6 = subset(res05.mdv2.3.l.6, padj<0.05)
print(resSig05.mdv2.3.l.6)
write.table(resSig05.mdv2.3.l.6,file="resSig05.mdv2.3.l.6.txt")
write.table(res.mdv2.3.l.6,file="res.mdv2.3.l.6.txt")

# Volcano plot
vp05<-EnhancedVolcano(res.mdv2.3.l.6,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 6 hpi") + theme(plot.title=element_text(size=20)) 
vp05

##### ##### MDV2.3  24 hpi       ####

mdv2.3.l.24<-cbind(local$MDV2.3.5_C_L_24,local$MDV2.3.6_C_L_24,local$MDV2.3.7_C_L_24,local$MDV2.3.43_C_L_24,
                   local$MDV2.3.13_O_L_24, local$MDV2.3.17_O_L_24, local$MDV2.3.30_O_L_24, local$MDV2.3.45_O_L_24)
colnames(mdv2.3.l.24)<-c("MDV2.3.5_C_L_24","MDV2.3.6_C_L_24","MDV2.3.7_C_L_24","MDV2.3.43_C_L_24",
                         "MDV2.3.15_O_L_24", "MDV2.3.17_O_L_24", "MDV2.3.30_O_L_24", "MDV2.3.45_O_L_24")
rownames(mdv2.3.l.24)<-rownames(local)
coldata.mdv2.3.l.24<-c("MDV2.3.5","MDV2.3.6","MDV2.3.7","MDV2.3.43","MDV2.3.13","MDV2.3.17","MDV2.3.30","MDV2.3.45",
                       "control","control","control","control","infected","infected","infected","infected")
coldata.mdv2.3.l.24<-matrix(coldata.mdv2.3.l.24,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.24)<-colnames(mdv2.3.l.24)
colnames(coldata.mdv2.3.l.24)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.24)

dds.mdv2.3.l.24 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.24,
                                          colData = coldata.mdv2.3.l.24,
                                          design = ~ TREATMENT)
dds.mdv2.3.l.24<-DESeq(dds.mdv2.3.l.24)
res.mdv2.3.l.24<-results(dds.mdv2.3.l.24)
res.mdv2.3.l.24
resOrdered.mdv2.3.l.24 <- res.mdv2.3.l.24[order(res.mdv2.3.l.24$pvalue),]
summary(res.mdv2.3.l.24)
sum(res.mdv2.3.l.24$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.24 <- results(dds.mdv2.3.l.24, alpha=0.05)
summary(res05.mdv2.3.l.24)
sum(res05.mdv2.3.l.24$padj < 0.05, na.rm=TRUE)

res01.mdv2.3.l.24 <- results(dds.mdv2.3.l.24, alpha=0.01)
summary(res01.mdv2.3.l.24)
sum(res01.mdv2.3.l.24$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.24 = subset(res05.mdv2.3.l.24, padj<0.05)
print(resSig05.mdv2.3.l.24)
write.table(resSig05.mdv2.3.l.24,file="resSig05.mdv2.3.l.24.txt")
write.table(res.mdv2.3.l.24,file="res.mdv2.3.l.24.txt")

# Volcano plot
vp06<-EnhancedVolcano(res.mdv2.3.l.24,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 24 hpi") + theme(plot.title=element_text(size=20)) 
vp06
##### ##### MDV2.3  72 hpi       ####

mdv2.3.l.72<-cbind(local$MDV2.3.3_C_L_72,local$MDV2.3.11_C_L_72,local$MDV2.3.41_C_L_72,local$MDV2.3.47_C_L_72,
                   local$MDV2.3.15_O_L_72, local$MDV2.3.28_O_L_72, local$MDV2.3.32_O_L_72, local$MDV2.3.44_O_L_72)
colnames(mdv2.3.l.72)<-c("MDV2.3.3_C_L_72","MDV2.3.11_C_L_72","MDV2.3.41_C_L_72","MDV2.3.47_C_L_72",
                         "MDV2.3.15_O_L_72", "MDV2.3.28_O_L_72", "MDV2.3.32_O_L_72", "MDV2.3.44_O_L_72")
rownames(mdv2.3.l.72)<-rownames(local)
coldata.mdv2.3.l.72<-c("MDV2.3.3","MDV2.3.11","MDV2.3.41","MDV2.3.47","MDV2.3.15","MDV2.3.28","MDV2.3.32","MDV2.3.44",
                       "control","control","control","control","infected","infected","infected","infected")
coldata.mdv2.3.l.72<-matrix(coldata.mdv2.3.l.72,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.72)<-colnames(mdv2.3.l.72)
colnames(coldata.mdv2.3.l.72)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.72)

dds.mdv2.3.l.72 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.72,
                                          colData = coldata.mdv2.3.l.72,
                                          design = ~ TREATMENT)
dds.mdv2.3.l.72<-DESeq(dds.mdv2.3.l.72)
res.mdv2.3.l.72<-results(dds.mdv2.3.l.72)
res.mdv2.3.l.72
resOrdered.mdv2.3.l.72 <- res.mdv2.3.l.72[order(res.mdv2.3.l.72$pvalue),]
summary(res.mdv2.3.l.72)
sum(res.mdv2.3.l.72$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.72 <- results(dds.mdv2.3.l.72, alpha=0.05)
summary(res05.mdv2.3.l.72)
sum(res05.mdv2.3.l.72$padj < 0.05, na.rm=TRUE)

res01.mdv2.3.l.72 <- results(dds.mdv2.3.l.72, alpha=0.01)
summary(res01.mdv2.3.l.72)
sum(res01.mdv2.3.l.72$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.72 = subset(res05.mdv2.3.l.72, padj<0.05)
print(resSig05.mdv2.3.l.72)
write.table(resSig05.mdv2.3.l.72,file="resSig05.mdv2.3.l.72.txt")
write.table(res.mdv2.3.l.72,file="res.mdv2.3.l.72.txt")

# Volcano plot
vp07<-EnhancedVolcano(res.mdv2.3.l.72,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 72 hpi") + theme(plot.title=element_text(size=20)) 
vp07
##### ##### MDV2.3  144 hpi       ####

mdv2.3.l.144<-cbind(local$MDV2.3.10_C_L_144,local$MDV2.3.27_C_L_144,local$MDV2.3.40_C_L_144,local$MDV2.3.42_C_L_144,
                    local$MDV2.3.19_O_L_144,local$MDV2.3.26_O_L_144,local$MDV2.3.46_O_L_144,local$MDV2.3.50_O_L_144)
colnames(mdv2.3.l.144)<-c("MDV2.3.10_C_L_144","MDV2.3.27_C_L_144","MDV2.3.40_C_L_144","MDV2.3.42_C_L_144",
                          "MDV2.3.19_O_L_144","MDV2.3.26_O_L_144","MDV2.3.46_O_L_144","MDV2.3.50_O_L_144")
rownames(mdv2.3.l.144)<-rownames(local)
coldata.mdv2.3.l.144<-c("MDV2.3.10","MDV2.3.27","MDV2.3.40","MDV2.3.42","MDV2.3.19","MDV2.3.26","MDV2.3.46","MDV2.3.50",
                        "control","control","control","control","infected","infected","infected","infected")
coldata.mdv2.3.l.144<-matrix(coldata.mdv2.3.l.144,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.mdv2.3.l.144)<-colnames(mdv2.3.l.144)
colnames(coldata.mdv2.3.l.144)<-c("SAMPLE","TREATMENT")
print(coldata.mdv2.3.l.144)

dds.mdv2.3.l.144 <- DESeqDataSetFromMatrix(countData = mdv2.3.l.144,
                                           colData = coldata.mdv2.3.l.144,
                                           design = ~ TREATMENT)
dds.mdv2.3.l.144<-DESeq(dds.mdv2.3.l.144)
res.mdv2.3.l.144<-results(dds.mdv2.3.l.144)
res.mdv2.3.l.144
resOrdered.mdv2.3.l.144 <- res.mdv2.3.l.144[order(res.mdv2.3.l.144$pvalue),]
summary(res.mdv2.3.l.144)
sum(res.mdv2.3.l.144$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l.144 <- results(dds.mdv2.3.l.144, alpha=0.05)
summary(res05.mdv2.3.l.144)
sum(res05.mdv2.3.l.144$padj < 0.05, na.rm=TRUE)

res01.mdv2.3.l.144 <- results(dds.mdv2.3.l.144, alpha=0.01)
summary(res01.mdv2.3.l.144)
sum(res01.mdv2.3.l.144$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l.144 = subset(res05.mdv2.3.l.144, padj<0.05)
print(resSig05.mdv2.3.l.144)
write.table(resSig05.mdv2.3.l.144,file="resSig05.mdv2.3.l.144.txt")
write.table(res.mdv2.3.l.144,file="res.mdv2.3.l.144.txt")

# Volcano plot
vp08<-EnhancedVolcano(res.mdv2.3.l.144,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("MDV2.3 144 hpi") + theme(plot.title=element_text(size=20)) 
vp08
##### ##### VAD2  6 hpi       ####

vad2.l.6<-cbind(local$VAD2.4_C_L_6,local$VAD2.6_C_L_6,local$VAD2.7_C_L_6,local$VAD2.37_C_L_6,
                local$VAD2.14_O_L_6,local$VAD2.15_O_L_6,local$VAD2.29_O_L_6,local$VAD2.44_O_L_6)
colnames(vad2.l.6)<-c("VAD2.4_C_L_6","VAD2.6_C_L_6","VAD2.7_C_L_6","VAD2.37_C_L_6",
                      "VAD2.14_O_L_6",  "VAD2.15_O_L_6",  "VAD2.29_O_L_6",  "VAD2.44_O_L_6")
rownames(vad2.l.6)<-rownames(local)
coldata.vad2.l.6<-c("VAD2.4","VAD2.6","VAD2.7","VAD2.37","VAD2.14","VAD2.15","VAD2.29","VAD2.44",
                    "control","control","control","control","infected","infected","infected","infected")
coldata.vad2.l.6<-matrix(coldata.vad2.l.6,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.6)<-colnames(vad2.l.6)
colnames(coldata.vad2.l.6)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.6)

dds.vad2.l.6 <- DESeqDataSetFromMatrix(countData = vad2.l.6,
                                       colData = coldata.vad2.l.6,
                                       design = ~ TREATMENT)
dds.vad2.l.6<-DESeq(dds.vad2.l.6)
res.vad2.l.6<-results(dds.vad2.l.6)
res.vad2.l.6
resOrdered.vad2.l.6 <- res.vad2.l.6[order(res.vad2.l.6$pvalue),]
summary(res.vad2.l.6)
sum(res.vad2.l.6$padj < 0.1, na.rm=TRUE)

res05.vad2.l.6 <- results(dds.vad2.l.6, alpha=0.05)
summary(res05.vad2.l.6)
sum(res05.vad2.l.6$padj < 0.05, na.rm=TRUE)

res01.vad2.l.6 <- results(dds.vad2.l.6, alpha=0.01)
summary(res01.vad2.l.6)
sum(res01.vad2.l.6$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.6 = subset(res05.vad2.l.6, padj<0.05)
print(resSig05.vad2.l.6)
write.table(resSig05.vad2.l.6,file="resSig05.vad2.l.6.txt")
write.table(res.vad2.l.6,file="res.vad2.l.6.txt")

# Volcano plot
vp09<-EnhancedVolcano(res.vad2.l.6,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 6 hpi") + theme(plot.title=element_text(size=20)) 
vp09

##### ##### VAD2  24 hpi       ####

vad2.l.24<-cbind(local$VAD2.1_C_L_24,local$VAD2.11_C_L_24,local$VAD2.19_C_L_24,local$VAD2.42_C_L_24,
                 local$VAD2.25_O_L_24,local$VAD2.26_O_L_24,local$VAD2.28_O_L_24,local$VAD2.45_O_L_24)
colnames(vad2.l.24)<-c("VAD2.1_C_L_24","VAD2.11_C_L_24","VAD2.19_C_L_24","VAD2.42_C_L_24",
                       "VAD2.25_O_L_24","VAD2.26_O_L_24","VAD2.28_O_L_24","VAD2.45_O_L_24")
rownames(vad2.l.24)<-rownames(local)
coldata.vad2.l.24<-c("VAD2.1","VAD2.11","VAD2.19","VAD2.42","VAD2.25","VAD2.26","VAD2.28","VAD2.45",
                     "control","control","control","control","infected","infected","infected","infected")
coldata.vad2.l.24<-matrix(coldata.vad2.l.24,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.24)<-colnames(vad2.l.24)
colnames(coldata.vad2.l.24)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.24)

dds.vad2.l.24 <- DESeqDataSetFromMatrix(countData = vad2.l.24,
                                        colData = coldata.vad2.l.24,
                                        design = ~ TREATMENT)
dds.vad2.l.24<-DESeq(dds.vad2.l.24)
res.vad2.l.24<-results(dds.vad2.l.24)
res.vad2.l.24
resOrdered.vad2.l.24 <- res.vad2.l.24[order(res.vad2.l.24$pvalue),]
summary(res.vad2.l.24)
sum(res.vad2.l.24$padj < 0.1, na.rm=TRUE)

res05.vad2.l.24 <- results(dds.vad2.l.24, alpha=0.05)
summary(res05.vad2.l.24)
sum(res05.vad2.l.24$padj < 0.05, na.rm=TRUE)

res01.vad2.l.24 <- results(dds.vad2.l.24, alpha=0.01)
summary(res01.vad2.l.24)
sum(res01.vad2.l.24$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.24 = subset(res05.vad2.l.24, padj<0.05)
print(resSig05.vad2.l.24)
write.table(resSig05.vad2.l.24,file="resSig05.vad2.l.24.txt")
write.table(res.vad2.l.24,file="res.vad2.l.24.txt")

# Volcano plot
vp10<-EnhancedVolcano(res.vad2.l.24,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 24 hpi") + theme(plot.title=element_text(size=20)) 
vp10
##### ##### VAD2  72 hpi       ####

vad2.l.72<-cbind(local$VAD2.2_C_L_72,local$VAD2.9_C_L_72,local$VAD2.27_C_L_72,local$VAD2.41_C_L_72,
                 local$VAD2.17_O_L_72,local$VAD2.22_O_L_72,local$VAD2.33_O_L_72, local$VAD2.48_O_L_72)
colnames(vad2.l.72)<-c("VAD2.2_C_L_72","VAD2.9_C_L_72","VAD2.27_C_L_72","VAD2.41_C_L_72",
                       "VAD2.17_O_L_72","VAD2.22_O_L_72","VAD2.33_O_L_72","VAD2.48_O_L_72")
rownames(vad2.l.72)<-rownames(local)
coldata.vad2.l.72<-c("VAD2.2","VAD2.9","VAD2.27","VAD2.41","VAD2.17","VAD2.22","VAD2.33","VAD2.48",
                     "control","control","control","control","infected","infected","infected","infected")
coldata.vad2.l.72<-matrix(coldata.vad2.l.72,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.72)<-colnames(vad2.l.72)
colnames(coldata.vad2.l.72)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.72)

dds.vad2.l.72 <- DESeqDataSetFromMatrix(countData = vad2.l.72,
                                        colData = coldata.vad2.l.72,
                                        design = ~ TREATMENT)
dds.vad2.l.72<-DESeq(dds.vad2.l.72)
res.vad2.l.72<-results(dds.vad2.l.72)
res.vad2.l.72
resOrdered.vad2.l.72 <- res.vad2.l.72[order(res.vad2.l.72$pvalue),]
summary(res.vad2.l.72)
sum(res.vad2.l.72$padj < 0.1, na.rm=TRUE)

res05.vad2.l.72 <- results(dds.vad2.l.72, alpha=0.05)
summary(res05.vad2.l.72)
sum(res05.vad2.l.72$padj < 0.05, na.rm=TRUE)

res01.vad2.l.72 <- results(dds.vad2.l.72, alpha=0.01)
summary(res01.vad2.l.72)
sum(res01.vad2.l.72$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.72 = subset(res05.vad2.l.72, padj<0.05)
print(resSig05.vad2.l.72)
write.table(resSig05.vad2.l.72,file="resSig05.vad2.l.72.txt")
write.table(res.vad2.l.72,file="res.vad2.l.72.txt")

# Volcano plot
vp11<-EnhancedVolcano(res.vad2.l.72,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 72 hpi") + theme(plot.title=element_text(size=20)) 
vp11

##### ##### VAD2  144 hpi       ####

vad2.l.144<-cbind(local$VAD2.8_C_L_144,local$VAD2.36_C_L_144,local$VAD2.40_C_L_144,local$VAD2.43_C_L_144,
                  local$VAD2.18_O_L_144,  local$VAD2.23_O_L_144,  local$VAD2.35_O_L_144,local$VAD2.47_O_L_144)
colnames(vad2.l.144)<-c("VAD2.8_C_L_144","VAD2.36_C_L_144","VAD2.40_C_L_144","VAD2.43_C_L_144",
                        "VAD2.18_O_L_144","VAD2.23_O_L_144","VAD2.35_O_L_144","VAD2.47_O_L_144")
rownames(vad2.l.144)<-rownames(local)
coldata.vad2.l.144<-c("VAD2.8","VAD2.36","VAD2.40","VAD2.43","VAD2.18","VAD2.23", "VAD2.35", "VAD2.47",
                      "control","control","control","control","infected","infected","infected","infected")
coldata.vad2.l.144<-matrix(coldata.vad2.l.144,nrow=8,ncol=2,byrow = FALSE)
rownames(coldata.vad2.l.144)<-colnames(vad2.l.144)
colnames(coldata.vad2.l.144)<-c("SAMPLE","TREATMENT")
print(coldata.vad2.l.144)
dds.vad2.l.144 <- DESeqDataSetFromMatrix(countData = vad2.l.144,
                                         colData = coldata.vad2.l.144,
                                         design = ~ TREATMENT)

dds.vad2.l.144<-DESeq(dds.vad2.l.144)
res.vad2.l.144<-results(dds.vad2.l.144)
res.vad2.l.144
resOrdered.vad2.l.144 <- res.vad2.l.144[order(res.vad2.l.144$pvalue),]
summary(res.vad2.l.144)
sum(res.vad2.l.144$padj < 0.1, na.rm=TRUE)

res05.vad2.l.144 <- results(dds.vad2.l.144, alpha=0.05)
summary(res05.vad2.l.144)
sum(res05.vad2.l.144$padj < 0.05, na.rm=TRUE)

res01.vad2.l.144 <- results(dds.vad2.l.144, alpha=0.01)
summary(res01.vad2.l.144)
sum(res01.vad2.l.144$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l.144 = subset(res05.vad2.l.144, padj<0.05)
print(resSig05.vad2.l.144)
write.table(resSig05.vad2.l.144,file="resSig05.vad2.l.144.txt")
write.table(res.vad2.l.144,file="res.vad2.l.144.txt")

# Volcano plot
vp12<-EnhancedVolcano(res.vad2.l.144,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'padj',
                      #title="(a) FAR3 DEGs",
                      subtitle = NULL,
                      #caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                      caption = bquote(" "),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize=1,col=c('black', 'black', 'black', 'red3'),
                      boxedLabels = F,
                      legendPosition = "none") +  
  #  ggplot2::coord_cartesian(xlim=c(-10, 10),ylim=c(0,15)) +
  ggtitle("VAD2 144 hpi") + theme(plot.title=element_text(size=20)) 
vp12

grid.arrange(vp01, vp05, vp09,vp02, vp06, vp10,vp03, vp07, vp11,vp04, vp08, vp12,
             nrow = 4)
grid.arrange(vp01, vp02, vp03,vp04, vp05, vp06,vp07, vp08, vp09,vp10, vp11, vp12,
             nrow = 3)

list.mdv1.6<-    as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv1.6.txt"))
list.mdv1.24<-   as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv1.24.txt"))
list.mdv1.72<-   as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv1.72.txt"))
list.mdv1.144<-  as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv1.144.txt"))
list.mdv2.3.6<-  as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv2.3.6.txt"))
list.mdv2.3.24<- as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv2.3.24.txt"))
list.mdv2.3.72<- as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv2.3.72.txt"))
list.mdv2.3.144<-as.matrix(read.table("ULMI.LOCAL/LISTS/list.mdv2.3.144.txt"))
list.vad2.6<-    as.matrix(read.table("ULMI.LOCAL/LISTS/list.vad2.6.txt"))
list.vad2.24<-   as.matrix(read.table("ULMI.LOCAL/LISTS/list.vad2.24.txt"))
list.vad2.72<-   as.matrix(read.table("ULMI.LOCAL/LISTS/list.vad2.72.txt"))
list.vad2.144<-  as.matrix(read.table("ULMI.LOCAL/LISTS/list.vad2.144.txt"))

listInput <- list(MDV1_6hpi = list.mdv1.6, MDV1_24hpi = list.mdv1.24, MDV1_72hpi = list.mdv1.72, MDV1_144hpi = list.mdv1.144, 
                  MDV2.3_6hpi = list.mdv2.3.6,MDV2.3_24hpi = list.mdv2.3.24, MDV2.3_72hpi = list.mdv2.3.72,MDV2.3_144hpi = list.mdv2.3.144,
                  VAD2_6hpi = list.vad2.6,VAD2_24hpi = list.vad2.24,VAD2_72hpi = list.vad2.72,VAD2_144hpi = list.vad2.144)


svg("my_upset_plot.svg",width=20, 
    height=10, 
    pointsize=12)
upset(fromList(listInput), 
      nsets = 12,
      nintersects = 120, 
      sets=NULL,
      order.by = "freq", 
      #order.by = "degree",
      #empty.intersections = "on",
      mb.ratio= c(0.6, 0.4), 
      group.by = "sets", 
      cutoff = 10,
      keep.order = T,
      number.angles = -30, 
      point.size = 3, line.size = 1.5, 
      #mainbar.y.label = "Intersection size", sets.x.label = "Set Size", 
      text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 1.2))
# Close the graphics device


## Barplots and Venn's diagrams for each genotype

# BARPLOTS

degs.bp <- data.frame(X = rep(c("MDV1 006hpi","MDV1 024hpi","MDV1 072hpi","MDV1 144hpi",
                                "MDV2.3 006hpi","MDV2.3 024hpi","MDV2.3 072hpi","MDV2.3 144hpi",
                                "VAD2 006hpi","VAD2 024hpi","VAD2 072hpi","VAD2 144hpi"),each=2),
                      Genotype = rep(c("MDV1","MDV2.3","VAD2"),each=8),
                      Time= rep(c("6 hpi","24 hpi","72 hpi", "144 hpi"),each=2,times=3),
                      Expression = rep(c("Induced","Repressed"),times=12),
                      DEGs = c(524,-270,935,-364,789,-210,656,-337,
                               370,-247,577,-716,656,-550,559,-205,
                               636,-759,573,-383,770,-355,714,-477))

barplot<-ggplot(degs.bp, aes(x=fct_inorder(Time), y=DEGs,fill=Expression))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  labs(title="(a)",x = "Time (hours post infection, hpi)", y = "Number of DEGs", fill = "") +
#  theme(axis.text.x = element_text(angle = 45))
  facet_wrap(~Genotype, nrow = 1) +
  geom_bar(position="stack", stat="identity")
barplot

# VENN DIAGRAMS PER GETNOYPE

list1<-list(
  'MDV1 local - 6 hpi' = as.character(list.mdv1.6),
  'MDV1 local - 24 hpi' = as.character(list.mdv1.24),
  'MDV1 local - 72 hpi' = as.character(list.mdv1.72),
  'MDV1 local - 144 hpi' = as.character(list.mdv1.144))

#ggVennDiagram(list[1:4])
#ggVennDiagram(list)

venn1<-ggVennDiagram(list1[1:4], label_alpha = 0,
                   category.names = c("6 hpi (794 DEGs)","24 hpi (1299 DEGs)",
                                      "72 hpi (999 DEGs)","144 hpi (993 DEGs)"),
                   show_intersect = FALSE,
                   set_size=2.5,
                   label_size = 2,
                   label_percent_digit = 1) + 
  scale_x_continuous(expand = expansion(mult=.13)) +
  scale_fill_distiller(palette = 9) + labs(title="(b)") +
  theme_void()
venn1
#v1 + scale_x_continuous(expand = expansion(mult=.13)) +
#  scale_fill_distiller(palette = 9) + labs(title="(b)") +
#  theme_void()



# (c)

list2<-list(
  'MDV2.3 local - 6 hpi' =   as.character(list.mdv2.3.6),
  'MDV2.3 local - 24 hpi' =  as.character(list.mdv2.3.24),
  'MDV2.3 local - 72 hpi' =  as.character(list.mdv2.3.72),
  'MDV2.3 local - 144 hpi' = as.character(list.mdv2.3.144))

#ggVennDiagram(list2[1:4])
#ggVennDiagram(list2)

venn2<-ggVennDiagram(list2[1:4], label_alpha = 0,
                   category.names = c("6 hpi (617 DEGs)","24 hpi (1293 DEGs)",
                                      "72 hpi (1206 DEGs)","144 hpi (764 DEGs)"),
                   show_intersect = FALSE,
                   set_size=2.5,
                   label_size = 2,
                   label_percent_digit = 1) +
  scale_x_continuous(expand = expansion(mult=.13)) +
  scale_fill_distiller(palette = 15) + labs(title="(c)") +
  theme_void()
venn2
#f4c+ scale_x_continuous(expand = expansion(mult=.13)) +
#  scale_fill_distiller(palette = 15) + labs(title="(c)") +
#  theme_void()

# (d)


list3<-list(
  'VAD2 local - 6 hpi' =   as.character(list.vad2.6),
  'VAD2 local - 24 hpi' =  as.character(list.vad2.24),
  'VAD2 local - 72 hpi' =  as.character(list.vad2.72),
  'VAD2 local - 144 hpi' = as.character(list.vad2.144))

#ggVennDiagram(list3[1:4])
#ggVennDiagram(list3)

venn3<-ggVennDiagram(list3[1:4], label_alpha = 0,
                   category.names = c("6 hpi (1395 DEGs)","24 hpi (956 DEGs)",
                                      "72 hpi (1125 DEGs)","144 hpi (1191 DEGs)"),
                   show_intersect = FALSE,
                   set_size=2.5,
                   label_size = 2,
                   label_percent_digit = 0) +
  scale_x_continuous(expand = expansion(mult=.13)) +
  scale_fill_distiller(palette = 17) + labs(title="(d)") +
  theme_void()
venn3
#f4d+ scale_x_continuous(expand = expansion(mult=.13)) +
#  scale_fill_distiller(palette = 17) + labs(title="(d)") +
#  theme_void()


## Grid plots

# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow=2, ncol=3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(f4a, vp=define_region(1:2,1))
print(f4b, vp = define_region(2, 1))
print(f4c, vp = define_region(2, 2))
print(f4d, vp = define_region(2, 3))


#####

################################################################################################
##########  5. A. Time-course of MDV1 LOCAL                              ##########

attach(local)
mdv1.l<-as.data.frame(cbind(MDV1.1_C_L_6, MDV1.3_C_L_6, MDV1.36_C_L_6, MDV1.47_C_L_6,
                            MDV1.4_C_L_24,MDV1.20_C_L_24,MDV1.22_C_L_24,MDV1.43_C_L_24,
                            MDV1.2_C_L_72,MDV1.11_C_L_72,MDV1.34_C_L_72,MDV1.42_C_L_72,
                            MDV1.8_C_L_144,MDV1.21_C_L_144,MDV1.37_C_L_144,MDV1.39_C_L_144,
                            MDV1.13_O_L_6, MDV1.15_O_L_6, MDV1.23_O_L_6, MDV1.50_O_L_6,
                            MDV1.16_O_L_24, MDV1.33_O_L_24,MDV1.41_O_L_24,MDV1.44_O_L_24,
                            MDV1.18_O_L_72,MDV1.31_O_L_72,MDV1.46_O_L_72,MDV1.49_O_L_72,
                            MDV1.9_O_L_144,MDV1.12_O_L_144,MDV1.24_O_L_144,MDV1.45_O_L_144))
detach()
rownames(mdv1.l)<-rownames(local)

coldata.mdv1.l<-c("MDV1.1","MDV1.3","MDV1.36","MDV1.47","MDV1.4","MDV1.20","MDV1.22","MDV1.43",
                  "MDV1.2","MDV1.11","MDV1.34","MDV1.42","MDV1.8","MDV1.21","MDV1.37","MDV1.39",
                  "MDV1.13","MDV1.15","MDV1.23","MDV1.50","MDV1.16","MDV1.33","MDV1.41","MDV1.44",
                  "MDV1.18","MDV1.31","MDV1.46","MDV1.49","MDV1.9","MDV1.12","MDV1.24","MDV1.45",
                  "control","control","control","control","control","control","control","control",
                  "control","control","control","control","control","control","control","control",
                  "infected","infected","infected","infected","infected","infected","infected","infected",
                  "infected","infected","infected","infected","infected","infected","infected","infected",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi",
                  "72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi",
                  "72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi",
                  "control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi",
                  "infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi",
                  "infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi")
coldata.mdv1.l<-matrix(coldata.mdv1.l,nrow=32,ncol=4,byrow=FALSE)
rownames(coldata.mdv1.l)<-colnames(mdv1.l)
colnames(coldata.mdv1.l)<-c("SAMPLE","TREATMENT","TIME","TREAT.TIME")

dds.mdv1.l <- DESeqDataSetFromMatrix(countData = mdv1.l,
                                     colData = coldata.mdv1.l,
                                     design = ~ TIME + TREATMENT + TIME:TREATMENT)

t_data.mdv1.l<-vst(dds.mdv1.l)
head(assay(t_data.mdv1.l))
distances.mdv1.l<-dist(t(assay(t_data.mdv1.l)))
distances_matrix.mdv1.l<-as.matrix(distances.mdv1.l)
rownames(distances_matrix.mdv1.l)<-paste(t_data.mdv1.l$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.mdv1.l<-hclust(distances.mdv1.l)
heatmap.2(distances_matrix.mdv1.l,Rowv = as.dendrogram(hc.mdv1.l),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE,
          key.title = "Color Key and Histogram",keysize = 1.5)
plotPCA(t_data.mdv1.l,intgroup="SAMPLE")
plotPCA(t_data.mdv1.l,intgroup="TREATMENT")
plotPCA(t_data.mdv1.l,intgroup="TIME")
plotPCA(t_data.mdv1.l,intgroup="TREAT.TIME")

pcaData.mdv1.l <- plotPCA(t_data.mdv1.l, intgroup=c("TIME", "TREATMENT"), returnData=TRUE)
percentVar.mdv1.l <- round(100 * attr(pcaData.mdv1.l, "percentVar"))
ggplot(pcaData.mdv1.l, aes(PC1, PC2, color=TREATMENT, shape=TIME)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.mdv1.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.mdv1.l[2],"% variance")) + 
  coord_fixed()

dds.mdv1.l<-DESeq(dds.mdv1.l,test = "LRT", 
                  reduced = ~ TIME + TREATMENT)
res.mdv1.l<-results(dds.mdv1.l)
res.mdv1.l
resOrdered.mdv1.l <- res.mdv1.l[order(res.mdv1.l$pvalue),]
summary(res.mdv1.l)
sum(res.mdv1.l$padj < 0.1, na.rm=TRUE)

res05.mdv1.l <- results(dds.mdv1.l, alpha=0.05)
summary(res05.mdv1.l)
sum(res05.mdv1.l$padj < 0.05, na.rm=TRUE)

res01.mdv1.l <- results(dds.mdv1.l, alpha=0.01)
summary(res01.mdv1.l)
sum(res01.mdv1.l$padj < 0.01, na.rm=TRUE)

resSig05.mdv1.l = subset(res05.mdv1.l, padj<0.05)
print(resSig05.mdv1.l)
write.table(resSig05.mdv1.l,file="resSig05.mdv1.local.txt")



################################################################################################
##########  5. B. Time-course of MDV2.3 LOCAL                              ##########

attach(local)
mdv2.3.l<-as.data.frame(cbind(MDV2.3.8_C_L_6, MDV2.3.24_C_L_6, MDV2.3.29_C_L_6, MDV2.3.38_C_L_6,
                              MDV2.3.5_C_L_24,MDV2.3.6_C_L_24,MDV2.3.7_C_L_24,MDV2.3.43_C_L_24,
                              MDV2.3.3_C_L_72,MDV2.3.11_C_L_72,MDV2.3.41_C_L_72,MDV2.3.47_C_L_72,
                              MDV2.3.10_C_L_144,MDV2.3.27_C_L_144,MDV2.3.40_C_L_144,MDV2.3.42_C_L_144,
                              MDV2.3.18_O_L_6, MDV2.3.21_O_L_6, MDV2.3.33_O_L_6, MDV2.3.49_O_L_6,
                              MDV2.3.13_O_L_24, MDV2.3.17_O_L_24,MDV2.3.30_O_L_24,MDV2.3.45_O_L_24,
                              MDV2.3.15_O_L_72,MDV2.3.28_O_L_72,MDV2.3.32_O_L_72,MDV2.3.44_O_L_72,
                              MDV2.3.19_O_L_144,MDV2.3.26_O_L_144,MDV2.3.46_O_L_144,MDV2.3.50_O_L_144))
detach()
rownames(mdv2.3.l)<-rownames(local)
coldata.mdv2.3.l<-c("MDV2.3.8","MDV2.3.24","MDV2.3.29","MDV2.3.38","MDV2.3.5","MDV2.3.6","MDV2.3.7","MDV2.3.43",
                    "MDV2.3.3","MDV2.3.11","MDV2.3.41","MDV2.3.47","MDV2.3.10","MDV2.3.27","MDV2.3.40","MDV2.3.42",
                    "MDV2.3.18","MDV2.3.21","MDV2.3.33","MDV2.3.49","MDV2.3.13","MDV2.3.17","MDV2.3.30","MDV2.3.45",
                    "MDV2.3.15","MDV2.3.28","MDV2.3.32","MDV2.3.44","MDV2.3.19","MDV2.3.26","MDV2.3.46","MDV2.3.50",
                    "control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",
                    "infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected",
                    "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                    "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                    "control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi",
                    "control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi",
                    "infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi",
                    "infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi")

coldata.mdv2.3.l<-matrix(coldata.mdv2.3.l,nrow=32,ncol=4,byrow=FALSE)
rownames(coldata.mdv2.3.l)<-colnames(mdv2.3.l)
colnames(coldata.mdv2.3.l)<-c("SAMPLE","TREATMENT","TIME","TREAT.TIME")

dds.mdv2.3.l <- DESeqDataSetFromMatrix(countData = mdv2.3.l,
                                       colData = coldata.mdv2.3.l,
                                       design = ~ TIME + TREATMENT + TIME:TREATMENT)
t_data.mdv2.3.l<-vst(dds.mdv2.3.l)
head(assay(t_data.mdv2.3.l))
distances.mdv2.3.l<-dist(t(assay(t_data.mdv2.3.l)))
distances_matrix.mdv2.3.l<-as.matrix(distances.mdv2.3.l)
rownames(distances_matrix.mdv2.3.l)<-paste(t_data.mdv2.3.l$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.mdv2.3.l<-hclust(distances.mdv2.3.l)
heatmap.2(distances_matrix.mdv2.3.l,Rowv = as.dendrogram(hc.mdv2.3.l),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE)
plotPCA(t_data.mdv2.3.l,intgroup="SAMPLE")
plotPCA(t_data.mdv2.3.l,intgroup="TREATMENT")
plotPCA(t_data.mdv2.3.l,intgroup="TIME")
plotPCA(t_data.mdv2.3.l,intgroup="TREAT.TIME")
pcaData.mdv2.3.l <- plotPCA(t_data.mdv2.3.l, intgroup=c("TIME", "TREATMENT"), returnData=TRUE)
percentVar.mdv2.3.l <- round(100 * attr(pcaData.mdv2.3.l, "percentVar"))
ggplot(pcaData.mdv2.3.l, aes(PC1, PC2, color=TREATMENT, shape=TIME)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.mdv2.3.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.mdv2.3.l[2],"% variance")) + 
  coord_fixed()

dds.mdv2.3.l<-DESeq(dds.mdv2.3.l,test = "LRT", 
                    reduced = ~ TIME + TREATMENT)
res.mdv2.3.l<-results(dds.mdv2.3.l)
res.mdv2.3.l
resOrdered.mdv2.3.l <- res.mdv2.3.l[order(res.mdv2.3.l$pvalue),]
summary(res.mdv2.3.l)
sum(res.mdv2.3.l$padj < 0.1, na.rm=TRUE)

res05.mdv2.3.l <- results(dds.mdv2.3.l, alpha=0.05)
summary(res05.mdv2.3.l)
sum(res05.mdv2.3.l$padj < 0.05, na.rm=TRUE)

res01.mdv2.3.l <- results(dds.mdv2.3.l, alpha=0.01)
summary(res01.mdv2.3.l)
sum(res01.mdv2.3.l$padj < 0.01, na.rm=TRUE)

resSig05.mdv2.3.l = subset(res05.mdv2.3.l, padj<0.05)
print(resSig05.mdv2.3.l)
write.table(resSig05.mdv2.3.l,file="resSig05.mdv2.3.local.txt")





################################################################################################
##########  5. C. Time-course of VAD2 LOCAL                             ##########

attach(local)
vad2.l<-as.data.frame(cbind(VAD2.4_C_L_6, VAD2.6_C_L_6, VAD2.7_C_L_6, VAD2.37_C_L_6,
                            VAD2.1_C_L_24,VAD2.11_C_L_24,VAD2.19_C_L_24,VAD2.42_C_L_24,
                            VAD2.2_C_L_72,VAD2.9_C_L_72,VAD2.27_C_L_72,VAD2.41_C_L_72,
                            VAD2.8_C_L_144,VAD2.36_C_L_144,VAD2.40_C_L_144,VAD2.43_C_L_144,
                            VAD2.14_O_L_6, VAD2.15_O_L_6, VAD2.29_O_L_6, VAD2.44_O_L_6,
                            VAD2.25_O_L_24, VAD2.26_O_L_24,VAD2.28_O_L_24,VAD2.45_O_L_24,
                            VAD2.17_O_L_72,VAD2.22_O_L_72,VAD2.33_O_L_72,VAD2.48_O_L_72,
                            VAD2.18_O_L_144,VAD2.23_O_L_144,VAD2.35_O_L_144,VAD2.47_O_L_144))
detach()
rownames(vad2.l)<-rownames(local)
coldata.vad2.l<-c("VAD2.4","VAD2.6","VAD2.7","VAD2.37","VAD2.1","VAD2.11","VAD2.19","VAD2.42",
                  "VAD2.2","VAD2.9","VAD2.27","VAD2.41","VAD2.8","VAD2.36","VAD2.40","VAD2.43",
                  "VAD2.14","VAD2.15","VAD2.29","VAD2.44","VAD2.25","VAD2.26","VAD2.28","VAD2.45",
                  "VAD2.17","VAD2.22","VAD2.33","VAD2.48","VAD2.18","VAD2.23","VAD2.35","VAD2.47",
                  "control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",
                  "infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "6 hpi","6 hpi","6 hpi","6 hpi","24 hpi","24 hpi","24 hpi","24 hpi","72 hpi","72 hpi","72 hpi","72 hpi","144 hpi","144 hpi","144 hpi","144 hpi",
                  "control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 006 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi","control - 024 hpi",
                  "control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 072 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi","control - 144 hpi",
                  "infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 006 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi","infected - 024 hpi",
                  "infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 072 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi","infected - 144 hpi")

coldata.vad2.l<-matrix(coldata.vad2.l,nrow=32,ncol=4,byrow=FALSE)
rownames(coldata.vad2.l)<-colnames(vad2.l)
colnames(coldata.vad2.l)<-c("SAMPLE","TREATMENT","TIME","TREAT.TIME")

dds.vad2.l <- DESeqDataSetFromMatrix(countData = vad2.l,
                                     colData = coldata.vad2.l,
                                     design = ~ TIME + TREATMENT + TIME:TREATMENT)
t_data.vad2.l<-vst(dds.vad2.l)
head(assay(t_data.vad2.l))
distances.vad2.l<-dist(t(assay(t_data.vad2.l)))
distances_matrix.vad2.l<-as.matrix(distances.vad2.l)
rownames(distances_matrix.vad2.l)<-paste(t_data.vad2.l$SAMPLE)
col<-colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
hc.vad2.l<-hclust(distances.vad2.l)
heatmap.2(distances_matrix.vad2.l,Rowv = as.dendrogram(hc.vad2.l),
          symm=TRUE,trace = "none",col = col,
          margins=c(2,10),labCol = FALSE)
plotPCA(t_data.vad2.l,intgroup="SAMPLE")
plotPCA(t_data.vad2.l,intgroup="TREATMENT")
plotPCA(t_data.vad2.l,intgroup="TIME")
plotPCA(t_data.vad2.l,intgroup="TREAT.TIME")
pcaData.vad2.l <- plotPCA(t_data.vad2.l, intgroup=c("TREATMENT", "TIME"), returnData=TRUE)
percentVar.vad2.l <- round(100 * attr(pcaData.vad2.l, "percentVar"))
ggplot(pcaData.vad2.l, aes(PC1, PC2, shape=TIME, color=TREATMENT)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.vad2.l[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vad2.l[2],"% variance")) + 
  coord_fixed()
dds.vad2.l<-DESeq(dds.vad2.l,test = "LRT",
                  reduced = ~ TIME + TREATMENT)
res.vad2.l<-results(dds.vad2.l)
res.vad2.l
resOrdered.vad2.l <- res.vad2.l[order(res.vad2.l$pvalue),]
summary(res.vad2.l)
sum(res.vad2.l$padj < 0.1, na.rm=TRUE)

res05.vad2.l <- results(dds.vad2.l, alpha=0.05)
summary(res05.vad2.l)
sum(res05.vad2.l$padj < 0.05, na.rm=TRUE)

res01.vad2.l <- results(dds.vad2.l, alpha=0.01)
summary(res01.vad2.l)
sum(res01.vad2.l$padj < 0.01, na.rm=TRUE)

resSig05.vad2.l = subset(res05.vad2.l, padj<0.05)
print(resSig05.vad2.l)
write.table(resSig05.vad2.l,file="resSig05.vad2.local.txt",sep="\t",dec=".")



################################################################################################
##########  6. Filtered DEGs from time course by LFC>|2|) in Excel                    ##########


################################################################################################
##########  7. Plot DEGs (also filtered by LFC>|2|) in Venn's, clustering, and plot heatmap ####

##########  7.A. Venn's diagram of DEGs ####
# These lists are DEGs without filetering by LFC 
#mdv1.l.tc.degs<-  read.table("list.mdv1.local.txt",stringsAsFactors = FALSE)
#mdv2.3.l.tc.degs<-read.table("list.mdv2.3.local.txt",stringsAsFactors = FALSE)
#vad2.l.tc.degs<-  read.table("list.vad2.local.txt",stringsAsFactors = FALSE)

# Now, we have the lists after LFC filtering (LFC>|2|), including 617 DEGs for MDV1, 817 DEGs for MDV2.3, 777 DEGs for VAD2  
mdv1.l.tc.degs<-  read.table("ULMI.LOCALISTS/degs.mdv1.time.course.list",stringsAsFactors = FALSE)
mdv2.3.l.tc.degs<-read.table("ULMI.LOCALISTS/degs.mdv2.3.time.course.list",stringsAsFactors = FALSE)
vad2.l.tc.degs<-  read.table("ULMI.LOCALISTS/degs.vad2.time.course.list",stringsAsFactors = FALSE)

mdv1.l.tc.degs<-  as.matrix(mdv1.l.tc.degs)
mdv2.3.l.tc.degs<-as.matrix(mdv2.3.l.tc.degs)
vad2.l.tc.degs<-  as.matrix(vad2.l.tc.degs)

list4<-list(
  'MDV2.3 local - time course'=as.character(mdv2.3.l.tc.degs),
  'MDV1 local - time course' = as.character(mdv1.l.tc.degs),
  'VAD2 local - time course' = as.character(vad2.l.tc.degs))

venn.tc<-ggVennDiagram(list4[1:3], label_alpha = 0,
                  category.names = c("MDV2.3      \n(817 DEGs)      ","MDV1\n(617 DEGs)",
                                     "      VAD2\n      (777 DEGs)"),
                  show_intersect = FALSE,
                  set_size=5,
                  label_size = 5,
                  label_percent_digit = 0) +
  scale_x_continuous(expand = expansion(mult=.13)) +
  scale_fill_distiller(palette = 3) + #labs(title="(d)") +
  theme_void()
venn.tc


##########  7.B. Clustering DEGs and ploting a single heatmap for the 3 genotypes ####
## DEGs are considered when pvalue < 0.05 for each independent genotype and LFC>|2| in at least one time of the specific genotype where is significant

alldegs<-read.delim("degs_time.course.local_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
alldegs<-as.matrix(alldegs)
alldegs.my_hclust_gene <- hclust(dist(alldegs), method = "ward.D2")
as.dendrogram(alldegs.my_hclust_gene) %>% plot(horiz = TRUE)
alldegs.my_gene_col<- cutree(tree=as.dendrogram(alldegs.my_hclust_gene),k=6)
alldegs.my_gene_col<- data.frame(Cluster=alldegs.my_gene_col)
table(alldegs.my_gene_col)

alldegs.my_sample_col<-data.frame(Time=rep(c("6 hpi","24 hpi","72 hpi","144 hpi",
                                             "6 hpi","24 hpi","72 hpi","144 hpi",
                                             "6 hpi","24 hpi","72 hpi","144 hpi"),
                                           c(1,1,1,1,1,1,1,1,1,1,1,1)),
                                  Genotype=rep(c("MDV1","MDV2.3","VAD2"),c(4,4,4)))
alldegs.my_gene_col$Cluster <- factor(alldegs.my_gene_col$Cluster,levels=c("1","2","3","4","5","6")) 

my_color=list(Cluster=c("1"="#8DD3C7","2"="#FB8072",
                        "3"="#80B1D3","4"="#FDB462","5"="#FCCDE5",
                        "6"="#D9D9D9"),
              Genotype=c("MDV1"="#FFED6F","MDV2.3"="#CCEBC5","VAD2"="#BC80BD"))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myColor2<- colorRampPalette(c("blue","white","red"))(300)
myBreaks1 <- c(seq(min(alldegs), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(alldegs)/paletteLength, max(alldegs), length.out=floor(paletteLength/2)))

myBreaks2 <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(10/paletteLength, 10, length.out=floor(paletteLength/2)))

pheatmap(alldegs,show_rownames = FALSE,cluster_cols = F,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         main="DEGs 2",
         annotation_row=alldegs.my_gene_col,
         annotation_col=alldegs.my_sample_col,
         cutree_rows = 6,#cutree_cols = alldegs.my_sample_col$Genotype,
         annotation_colors = my_color,
         cellwidth = 30, 
         #cex = 1,
         color = myColor, 
         #scale="row",
         breaks = myBreaks2,
         #legend_breaks = c(-5,-2,0,2,5),
         ##annotation_legend = TRUE,
         #angle_col = 45,
         #annotation_names_row = F
         )
write.table(mdv1.degs.dm.my_gene_col, file="mdv1.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")

##########  7.C. Clustering MDV1 DEGs and plotting heatmap ####
## DEGs are considered when pvalue < 0.05 for MDV1 genotype and LFC>|2| in at least one time

mdv1.degs<-read.delim("ULMI.LOCAL/degs_mdv1_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
mdv1.degs<-as.matrix(mdv1.degs)
mdv1.degs.my_hclust_gene <- hclust(dist(mdv1.degs), method = "ward.D")
as.dendrogram(mdv1.degs.my_hclust_gene) %>% plot(horiz = TRUE)
mdv1.degs.my_gene_col<- cutree(tree=as.dendrogram(mdv1.degs.my_hclust_gene),k=6)
mdv1.degs.my_gene_col<- data.frame(Cluster=mdv1.degs.my_gene_col)
table(mdv1.degs.my_gene_col)

mdv1.degs.my_sample_col<-data.frame(Time=c("6 hpi","24 hpi","72 hpi","144 hpi"))
mdv1.degs.my_gene_col$Cluster <- factor(mdv1.degs.my_gene_col$Cluster,levels=c("1","2","3","4","5","6")) 

#my_color=list(Cluster=c("1"="#8DD3C7","2"="#FB8072","3"="#80B1D3","4"="#FDB462","5"="#FCCDE5","6"="#D9D9D9"), 
#              Time= c("6 hpi"="#FFED6F","24 hpi"="#CCEBC5","72 hpi"="#BC80BD","144 hpi"="black"))

my_color=list(Cluster=c("1"="#CCEBC5","2"="#FB8072","3"="#80B1D3","4"="#FDB462","5"="#FCCDE5","6"="#FFED6F"))#, 
#              Time= c("6 hpi"="#FFED6F","24 hpi"="#CCEBC5","72 hpi"="#BC80BD","144 hpi"="black"))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
#myColor2<- colorRampPalette(c("blue","white","red"))(300)

#myBreaks1 <- c(seq(min(mdv1.degs), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(mdv1.degs)/paletteLength, max(mdv1.degs), length.out=floor(paletteLength/2))) ## setting max and min in 10 and -10

myBreaks.mdv1 <- c(seq(-15, 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(15/paletteLength, 15, length.out=floor(paletteLength/2)))

p1<-pheatmap(mdv1.degs,show_rownames = FALSE,cluster_cols = F,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         main="DEGs MDV1",
         annotation_row=mdv1.degs.my_gene_col,
#         annotation_col=mdv1.degs.my_sample_col,
         cutree_rows = 6,
         annotation_colors = my_color,
         cellwidth = 30, 
         cex = 1,
         color = myColor, 
         breaks = myBreaks.mdv1,
         legend_breaks = c(-12,-6,0,6,12),
         annotation_legend = TRUE,
         angle_col = 45,
         annotation_names_row = F,annotation_names_col = F)

mdv1.degs.my_gene_col
write.table(mdv1.degs.my_gene_col, file="mdv1.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")

##########  7.D. Clustering MDV2.3 DEGs and plotting heatmap ####
## DEGs are considered when pvalue < 0.05 for mdv2.3 genotype and LFC>|2| in at least one time

mdv2.3.degs<-read.delim("ULMI.LOCAL/degs_mdv2.3_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
mdv2.3.degs<-as.matrix(mdv2.3.degs)
mdv2.3.degs.my_hclust_gene <- hclust(dist(mdv2.3.degs), method = "ward.D")
as.dendrogram(mdv2.3.degs.my_hclust_gene) %>% plot(horiz = TRUE)
mdv2.3.degs.my_gene_col<- cutree(tree=as.dendrogram(mdv2.3.degs.my_hclust_gene),k=6)
mdv2.3.degs.my_gene_col<- data.frame(Cluster=mdv2.3.degs.my_gene_col)
table(mdv2.3.degs.my_gene_col)

#mdv2.3.degs.my_sample_col<-data.frame(Time=rep(c("6 hpi","24 hpi","72 hpi","144 hpi"),c(1,1,1,1)))
mdv2.3.degs.my_gene_col$Cluster <- factor(mdv2.3.degs.my_gene_col$Cluster,levels=c("1","2","3","4","5","6")) 

myBreaks.mdv2.3 <- c(seq(-15, 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(15/paletteLength, 15, length.out=floor(paletteLength/2))) ## setting max and min in 10 and -10

p2<-pheatmap(mdv2.3.degs,show_rownames = FALSE,cluster_cols = F,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D",
             main="DEGs MDV2.3",
             annotation_row=mdv2.3.degs.my_gene_col,
             #         annotation_col=mdv1.degs.my_sample_col,
             cutree_rows = 6,
             annotation_colors = my_color,
             cellwidth = 30, 
             cex = 1,
             color = myColor, 
             breaks = myBreaks.mdv2.3,
             legend_breaks = c(-12,-6,0,6,12),
             annotation_legend = TRUE,
             angle_col = 45,
             annotation_names_row = F,annotation_names_col = F)

mdv2.3.degs.my_gene_col
write.table(mdv2.3.degs.my_gene_col, file="mdv2.3.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")


##########  7.E. Clustering VAD2 DEGs and plotting heatmap ####
## DEGs are considered when pvalue < 0.05 for VAD2 genotype and LFC>|2| in at least one time

vad2.degs<-read.delim("ULMI.LOCAL/degs_vad2_heatmap.txt",header=T,dec=".",sep="\t",check.names = FALSE,row.names = 1)
vad2.degs<-as.matrix(vad2.degs)
vad2.degs.my_hclust_gene <- hclust(dist(vad2.degs), method = "ward.D")
as.dendrogram(vad2.degs.my_hclust_gene) %>% plot(horiz = TRUE)
vad2.degs.my_gene_col<- cutree(tree=as.dendrogram(vad2.degs.my_hclust_gene),k=6)
vad2.degs.my_gene_col<- data.frame(Cluster=vad2.degs.my_gene_col)
table(vad2.degs.my_gene_col)

#vad2.degs.my_sample_col<-data.frame(Time=rep(c("6 hpi","24 hpi","72 hpi","144 hpi"),c(1,1,1,1)))
vad2.degs.my_gene_col$Cluster <- factor(vad2.degs.my_gene_col$Cluster,levels=c("1","2","3","4","5","6")) 

myBreaks.vad2 <- c(seq(-15, 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(15/paletteLength, 15, length.out=floor(paletteLength/2))) ## setting max and min in 10 and -10

p3<-pheatmap(vad2.degs,show_rownames = FALSE,cluster_cols = F,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D",
             main="DEGs VAD2",
             annotation_row=vad2.degs.my_gene_col,
             #         annotation_col=vad2.degs.my_sample_col,
             cutree_rows = 6,
             annotation_colors = my_color,
             cellwidth = 30, 
             cex = 1,
             color = myColor, 
             breaks = myBreaks.vad2,
             legend_breaks = c(-12,-6,0,6,12),
             annotation_legend = TRUE,
             angle_col = 45,
             annotation_names_row = F,annotation_names_col = F)
vad2.degs.my_gene_col
write.table(vad2.degs.my_gene_col, file="vad2.l.tc.tableclust.txt", row.names=TRUE, col.names=TRUE,sep="\t",dec=".")

##########  7.F. Plotting Heatmaps of 3 genotypes together in grid ####
## DEGs are considered when pvalue < 0.05 for each genotype and LFC>|2| in at least one time
## There are three heatmaps together, one per genotype, not one heatmap for all three genotypes together

plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2']]=p2[[4]]
plot_list[['p3']]=p3[[4]]
plot_list
h3<-grid.arrange(grobs=plot_list, ncol=3)

##########  8. Plotting cluster graphs for the 3 genotypes ####
## DEGs for each cluster. Clusters are numbered from 1 to 6, and genotypes as A, B and C for MDV1, MDV2.3 and VAD2 respectively
## There are three heatmaps together, one per genotype, not one heatmap for all three genotypes together

#### Cluster A1 ####
a1<-read.delim("mdv1.cluster1.txt",header=F,dec=".",sep="\t")
rnames.a1<-a1[,1]
cl.a1<-data.frame(a1[,2:ncol(a1)])
rownames(cl.a1)<-rnames.a1
colnames(cl.a1)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A2 ####
a2<-read.delim("mdv1.cluster2.txt",header=F,dec=".",sep="\t")
rnames.a2<-a2[,1]
cl.a2<-data.frame(a2[,2:ncol(a2)])
rownames(cl.a2)<-rnames.a2
colnames(cl.a2)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A3 ####
a3<-read.delim("mdv1.cluster3.txt",header=F,dec=".",sep="\t")
rnames.a3<-a3[,1]
cl.a3<-data.frame(a3[,2:ncol(a3)])
rownames(cl.a3)<-rnames.a3
colnames(cl.a3)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A4 ####
a4<-read.delim("mdv1.cluster4.txt",header=F,dec=".",sep="\t")
rnames.a4<-a4[,1]
cl.a4<-data.frame(a4[,2:ncol(a4)])
rownames(cl.a4)<-rnames.a4
colnames(cl.a4)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A5 ####
a5<-read.delim("mdv1.cluster5.txt",header=F,dec=".",sep="\t")
rnames.a5<-a5[,1]
cl.a5<-data.frame(a5[,2:ncol(a5)])
rownames(cl.a5)<-rnames.a5
colnames(cl.a5)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster A6 ####
a6<-read.delim("mdv1.cluster6.txt",header=F,dec=".",sep="\t")
rnames.a6<-a6[,1]
cl.a6<-data.frame(a6[,2:ncol(a6)])
rownames(cl.a6)<-rnames.a6
colnames(cl.a6)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B1 ####
b1<-read.delim("mdv2.3.cluster1.txt",header=F,dec=".",sep="\t")
rnames.b1<-b1[,1]
cl.b1<-data.frame(b1[,2:ncol(b1)])
rownames(cl.b1)<-rnames.b1
colnames(cl.b1)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B2 ####
b2<-read.delim("mdv2.3.cluster2.txt",header=F,dec=".",sep="\t")
rnames.b2<-b2[,1]
cl.b2<-data.frame(b2[,2:ncol(b2)])
rownames(cl.b2)<-rnames.b2
colnames(cl.b2)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B3 ####
b3<-read.delim("mdv2.3.cluster3.txt",header=F,dec=".",sep="\t")
rnames.b3<-b3[,1]
cl.b3<-data.frame(b3[,2:ncol(b3)])
rownames(cl.b3)<-rnames.b3
colnames(cl.b3)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B4 ####
b4<-read.delim("mdv2.3.cluster4.txt",header=F,dec=".",sep="\t")
rnames.b4<-b4[,1]
cl.b4<-data.frame(b4[,2:ncol(b4)])
rownames(cl.b4)<-rnames.b4
colnames(cl.b4)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B5 ####
b5<-read.delim("mdv2.3.cluster5.txt",header=F,dec=".",sep="\t")
rnames.b5<-b5[,1]
cl.b5<-data.frame(b5[,2:ncol(b5)])
rownames(cl.b5)<-rnames.b5
colnames(cl.b5)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster B6 ####
b6<-read.delim("mdv2.3.cluster6.txt",header=F,dec=".",sep="\t")
rnames.b6<-b6[,1]
cl.b6<-data.frame(b6[,2:ncol(b6)])
rownames(cl.b6)<-rnames.b6
colnames(cl.b6)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C1 ####
c1<-read.delim("vad2.cluster1.txt",header=F,dec=".",sep="\t")
rnames.c1<-c1[,1]
cl.c1<-data.frame(c1[,2:ncol(c1)])
rownames(cl.c1)<-rnames.c1
colnames(cl.c1)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C2 ####
c2<-read.delim("vad2.cluster2.txt",header=F,dec=".",sep="\t")
rnames.c2<-c2[,1]
cl.c2<-data.frame(c2[,2:ncol(c2)])
rownames(cl.c2)<-rnames.c2
colnames(cl.c2)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C3 ####
c3<-read.delim("vad2.cluster3.txt",header=F,dec=".",sep="\t")
rnames.c3<-c3[,1]
cl.c3<-data.frame(c3[,2:ncol(c3)])
rownames(cl.c3)<-rnames.c3
colnames(cl.c3)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C4 ####
c4<-read.delim("vad2.cluster4.txt",header=F,dec=".",sep="\t")
rnames.c4<-c4[,1]
cl.c4<-data.frame(c4[,2:ncol(c4)])
rownames(cl.c4)<-rnames.c4
colnames(cl.c4)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C5 ####
c5<-read.delim("vad2.cluster5.txt",header=F,dec=".",sep="\t")
rnames.c5<-c5[,1]
cl.c5<-data.frame(c5[,2:ncol(c5)])
rownames(cl.c5)<-rnames.c5
colnames(cl.c5)<-c("6 hpi","24 hpi","72 hpi","144 hpi")
#### Cluster C6 ####
c6<-read.delim("vad2.cluster6.txt",header=F,dec=".",sep="\t")
rnames.c6<-c6[,1]
cl.c6<-data.frame(c6[,2:ncol(c6)])
rownames(cl.c6)<-rnames.c6
colnames(cl.c6)<-c("6 hpi","24 hpi","72 hpi","144 hpi")

#### Plot 18 clusters ####

par(mfrow = c(2, 3))
ga1<-matplot(t(cl.a1),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a1)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 1 (100 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a1)))

ga2<-matplot(t(cl.a2),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a2)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 2 (79 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a2)))

ga3<-matplot(t(cl.a3),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a3)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 3 (78 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a3)))

ga4<-matplot(t(cl.a4),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a4)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 4 (149 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a4)))

ga5<-matplot(t(cl.a5),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a5)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 5 (110 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a5)))

ga6<-matplot(t(cl.a6),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.a6)),col="red",lwd=2.5)+
  title("MDV1 - Cluster 6 (101 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.a6)))

gb1<-matplot(t(cl.b1),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b1)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 1 (148 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b1)))

gb2<-matplot(t(cl.b2),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b2)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 2 (267 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b2)))

gb3<-matplot(t(cl.b3),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b3)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 3 (178 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b3)))

gb4<-matplot(t(cl.b4),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b4)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 4 (55 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b4)))

gb5<-matplot(t(cl.b5),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b5)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 5 (79 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b5)))

gb6<-matplot(t(cl.b6),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.b6)),col="red",lwd=2.5)+
  title("MDV2.3 - Cluster 6 (90 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.b6)))

gc1<-matplot(t(cl.c1),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c1)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 1 (171 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c1)))

gc2<-matplot(t(cl.c2),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c2)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 2 (118 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c2)))

gc3<-matplot(t(cl.c3),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c3)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 3 (185 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c3)))

gc4<-matplot(t(cl.c4),type='l',ylab='Fold Change',col="grey",
              lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c4)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 4 (120 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c4)))

gc5<-matplot(t(cl.c5),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c5)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 5 (127 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c5)))

gc6<-matplot(t(cl.c6),type='l',ylab='Fold Change',col="grey",
             lwd=1,xaxt='n',cex.lab=1,cex.axis=1)+
  lines(rowMeans(t(cl.c6)),col="red",lwd=2.5)+
  title("VAD2 - Cluster 6 (56 DEGs)",adj  = 0,cex.main=1.2) +
  axis(1,1:4,labels=rownames(t(cl.c6)))

