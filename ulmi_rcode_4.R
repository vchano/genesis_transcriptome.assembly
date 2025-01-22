# Load necessary libraries
library(ggplot2)
library(forcats)

degs.bp <- data.frame(X = rep(c("MDV1 006hpi","MDV1 024hpi","MDV1 072hpi","MDV1 144hpi",
                           "MDV2.3 006hpi","MDV2.3 024hpi","MDV2.3 072hpi","MDV2.3 144hpi",
                           "VAD2 006hpi","VAD2 024hpi","VAD2 072hpi","VAD2 144hpi"),each=2),
                 Genotype = rep(c("MDV1","MDV2.3","VAD2"),each=8),
                 Time= rep(c("6 hpi","24 hpi","72 hpi", "144 hpi"),each=2,times=3),
                 Expression = rep(c("Induced","Repressed"),times=12),
                 DEGs = c(506,-15,1075,-99,796,-35,715,-56,
                          218,-12,650,-87,489,-22,468,-7,
                          830,-194,745,-20,955,-33,884,-32))

# Create the plot
barplot <- ggplot(degs.bp, aes(x = fct_inorder(Time), y = DEGs, fill = Expression)) + 
  theme_bw() +
  labs(x = "Time (hours post infection, hpi)", y = "Number of DEGs", fill = "") +
  facet_wrap(~Genotype, nrow = 1) +
  geom_bar(position = "stack", stat = "identity") +
  theme(
    plot.title = element_text(size = 1.4 * 20),  # increased font size by 40%
    axis.title.x = element_text(size = 1.4 * 15),
    axis.title.y = element_text(size = 1.4 * 15),
    axis.text.x = element_text(size = 1.4 * 10),
    axis.text.y = element_text(size = 1.4 * 10),
    legend.title = element_text(size = 1.4 * 10),
    legend.text = element_text(size = 1.4 * 10),
    strip.text = element_text(size = 1.4 * 15)  # increased font size for facet titles
  )

# Save the plot as a TIFF file with 300 dpi and specified dimensions
ggsave("ulmi_sci.data_figure.4.tiff", plot = barplot, width = 12, height = 8, dpi = 300, units = "in")

