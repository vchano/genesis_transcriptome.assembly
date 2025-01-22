
# Load necessary libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

# (a)
length <- read.table("1a.seqs.by.length.txt", header = TRUE)
df.length <- as.data.frame(length)
f2a <- ggplot(data = length, aes(x = Length, y = Seqs)) +
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

# (b)
annot.dist <- read.table("1b.data.distrib.txt", header = TRUE)

# Create the pie chart
f2b <- ggplot(annot.dist, aes(x = "", y = Seqs, fill = Annotation_progress)) +
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

# (c)
go.dist <- read.table("1c.GO.level.distrib.txt", header = TRUE)

# Replace Component, Function, and Process with CC, MF, and BP
go.dist$Category <- gsub("Component", "CC", go.dist$Category)
go.dist$Category <- gsub("Function", "MF", go.dist$Category)
go.dist$Category <- gsub("Process", "BP", go.dist$Category)

f2c <- ggplot(go.dist, aes(GO_Level, Annotations, fill = Category)) +
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

## Arrange the plots
combined_plot <- grid.arrange(f2a, f2b, f2c, ncol = 2, layout_matrix = rbind(c(1, 2), c(3, 3)))

# Save the combined figure
ggsave("ulmi_sci.data_figure.1.tiff", plot = combined_plot, width = 18, height = 12, dpi = 300, units = "in")
