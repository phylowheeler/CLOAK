library(devtools)
library(patchwork)
library(sem)
library(ggfortify)
library(ggrepel)
library(expm)
library(ca)
library(FactoMineR)

#Figure 1
df_pca <- prcomp(allQ)
df_out <- as.data.frame(df_pca$x)
df_out$group = species
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(df_out,aes(x=PC1,y=PC2, color=group, label = filenames)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=30)) + 
  geom_point() +
  geom_text_repel(size=3) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  xlim(-2.5,1.5) +
  ylim(-2.5,1.5)
ggsave(filename = "archaea.jpg", width=5, height=5, units="in")

#Figure 2
df_pca <- prcomp(allE)
df_out <- as.data.frame(df_pca$x)
df_out$group = species
df_loadings <- as.data.frame(df_pca$rotation)
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(df_out,aes(x=PC1,y=PC2,color=group, label = filenames)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=15)) + 
  geom_point() +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
ggsave(filename = "Transition_Eigen_Decomposition_1.jpg", width=5, height=6, units="in")
ggplot(df_out, aes(x=PC1,y=PC2, color=group, label=filenames)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=15)) + 
  geom_point() +
  geom_point(data = subset(df_loadings, row.names(df_loadings) %in% colnames1), aes(x = PC1, y = PC2), color = "red", inherit.aes=FALSE) +
  geom_text_repel(data = subset(df_loadings, row.names(df_loadings) %in% colnames1), aes(x = PC1, y = PC2), label = colnames1, color = "red", inherit.aes=FALSE, size=5) +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  geom_hline(aes(yintercept = 0),linetype="dashed") +
  geom_vline(aes(xintercept = 0),linetype="dashed") +
  xlim(-0.6,0.5) +
  ylim(-0.3,0.8)
ggsave(filename = "Transition_Eigen_Decomposition_1_biplot_first.jpg", width=5, height=6, units="in")
ggplot(df_out, aes(x=PC1,y=PC2, color=group, label=filenames)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=15)) + 
  geom_point() +
  geom_point(data = subset(df_loadings, row.names(df_loadings) %in% colnames2), aes(x = PC1, y = PC2), color = "red", inherit.aes=FALSE) +
  geom_text_repel(data = subset(df_loadings, row.names(df_loadings) %in% colnames2), aes(x = PC1, y = PC2), label = colnames2, color = "red", inherit.aes=FALSE, size=5) +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  geom_hline(aes(yintercept = 0),linetype="dashed") +
  geom_vline(aes(xintercept = 0),linetype="dashed")
ggsave(filename = "Transition_Eigen_Decomposition_1_biplot_second.jpg", width=5, height=6, units="in")

#Figure 3
df_pca <- prcomp(allE[1:7,])
df_out <- as.data.frame(df_pca$x)
df_loadings <- as.data.frame(df_pca$rotation)
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(df_out,aes(x=PC1,y=PC2, label = filenames[1:7])) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=15)) + 
  geom_point() +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2])
ggsave(filename = "Transition_Eigen_Decomposition_2.jpg", width=5, height=6, units="in")
ggplot(df_out, aes(x=PC1,y=PC2, label=filenames[1:7])) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=15)) + 
  geom_point() +
  geom_point(data = subset(df_loadings, row.names(df_loadings) %in% colnames1), aes(x = PC1, y = PC2), color = "red", inherit.aes=FALSE) +
  geom_text_repel(data = subset(df_loadings, row.names(df_loadings) %in% colnames1), aes(x = PC1, y = PC2), label = colnames1, color = "red", inherit.aes=FALSE, size=5) +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  geom_hline(aes(yintercept = 0),linetype="dashed") +
  geom_vline(aes(xintercept = 0),linetype="dashed")
ggsave(filename = "Transition_Eigen_Decomposition_2_biplot_first.jpg", width=5, height=6, units="in")
ggplot(df_out, aes(x=PC1,y=PC2, label=filenames[1:7])) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=15)) + 
  geom_point() +
  geom_point(data = subset(df_loadings, row.names(df_loadings) %in% colnames2), aes(x = PC1, y = PC2), color = "red", inherit.aes=FALSE) +
  geom_text_repel(data = subset(df_loadings, row.names(df_loadings) %in% colnames2), aes(x = PC1, y = PC2), label = colnames2, color = "red", inherit.aes=FALSE, size=5) +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  geom_hline(aes(yintercept = 0),linetype="dashed") +
  geom_vline(aes(xintercept = 0),linetype="dashed")
ggsave(filename = "Transition_Eigen_Decomposition_2_biplot_second.jpg", width=5, height=6, units="in")



# Define symbols for each point
symbols = ['o', 's', 'D', 'v', 'P', '*', 'X', 'h', '^', '+', 'p']

# Extract x and y coordinates
x_coords = [point[0] for point in points.values()]
y_coords = [point[1] for point in points.values()]

# Plot the points with symbols and labels
for i, (name, point) in enumerate(points.items()):
    if name in ['4_files', '7_files']:
        offset = 0.0025 if name == '4_files' else -0.0025
        plt.plot(point[0] + offset, point[1], marker=symbols[i], markersize=12, label=name)
    elif name == 'Ground Truth':
        plt.plot(point[0], point[1], marker='s', markersize=15, label=name, alpha=0.7)  # Adjusted marker style and size for 'Ground Truth'
        plt.annotate('Ground Truth', xy=(point[0], point[1]), xytext=(point[0]+0.05, point[1]-0.15),
                     arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)
    else:
        plt.plot(point[0], point[1], marker=symbols[i], markersize=12, label=name)

# Set the axis limits
plt.xlim(0, 0.3)
plt.ylim(0.4, 1)

# Set the axis labels
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')

# Create a legend
plt.legend(loc='upper right', fontsize='small')  # Adjusted legend position to 'upper right' and font size to 'small'

# Set the title
plt.title('')

# Display the plot
plt.show()
#end
