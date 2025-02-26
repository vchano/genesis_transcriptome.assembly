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
go.dist$Category <- gsub("CC", "", go.dist$Category)
go.dist$Category <- gsub("MF", "", go.dist$Category)
go.dist$Category <- gsub("BP", "", go.dist$Category)

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

