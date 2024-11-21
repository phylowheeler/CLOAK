library(expm)
library(factoextra)

# matrices=c("Q.insect","Q.plant","Q.yeast","Q.mammal","Q.bird","Q.insect_cloak","Q.plant_cloak","Q.yeast_cloak","Q.mammal_cloak","Q.bird_cloak")



mitochondrial = c("mtREV", "mtMAM", "mtART", "mtZOA", "mtMet" , "mtVer" , "mtInv")
chloroplast = c("cpREV")
viral = c("HIVb", "HIVw", "FLU", "rtREV")
nuclear = c("WAG", "Dayhoff","JTT", "LG", "VT", "PMB", "Blosum62")
taxa = c("Q.pfam","Q.plant","Q.bird", "Q.mammal", "Q.insect", "Q.yeast","Q.bac","Q.arch")
nonrev = c("NQ.pfam","NQ.plant","NQ.bird", "NQ.mammal", "NQ.insect", "NQ.yeast","NQ.bac","NQ.arch")

matrices=c(nuclear,taxa)

exchangabilities = I(vector('list', 35))
equilibrium_frequencies = I(vector('list', 35))
Q_matrices = I(vector('list', 35))
normalized_Q_matrices = I(vector('list', 35))
Q_eigen_decomposition = I(vector('list', 35))
transition_matrices = I(vector('list', 35))
eigen_decomposition = I(vector('list', 35))
matrix_table = data.frame(matrices,exchangabilities,equilibrium_frequencies,Q_matrices,normalized_Q_matrices,transition_matrices,eigen_decomposition)

allQ = NULL

for(i in 1:35){
	#import exchangabilities and equilibrium frequency vector
	matrix = matrices[i]
	handle = file(matrix,'rt')
	raw = readLines(handle)
	exch = data.matrix(read.table(text=raw,fill=TRUE,col.names=paste("V", 1:20)))
	close(handle)

	if (i <= 27){
		exch = rbind(NA,exch)
		pi = exch[21,]
		matrix_table$equilibrium_frequencies[[i]] = pi
		exch = exch[-21,]
		matrix_table$exchangabilities[[i]] = exch
		# Generate Q matrix
		Q = exch
		Q[upper.tri(Q)] = t(Q)[upper.tri(Q)]
		for (j in 1:20){Q[,j] = Q[,j]*pi[j]}
		for (j in 1:20){Q[j,j] = -sum(Q[j,],na.rm=TRUE)}
		matrix_table$Q_matrices[[i]] = Q
	}

	if (i > 27){
		pi = exch[21,]
		matrix_table$equilibrium_frequencies[[i]] = pi
		exch = exch[-22,]
		exch = exch[-21,]
		# Generate Q matrix
		Q = exch
		matrix_table$Q_matrices[[i]] = Q
	}

	#normalization
	mu=0
	for (j in 1:20){mu=mu-Q[j,j]*pi[j]}
	norm = Q*mu
	matrix_table$normalized_Q_matrices[[i]] = norm
	allQ = rbind(allQ, as.vector(t(norm)))
	Q_eigen = eigen(t(norm))
	matrix_table$Q_eigen_decomposition[[i]] = Q_eigen

	#exponentiation
	trans = expm(norm)
	matrix_table$transition_matrices[[i]] = trans
	eigen = eigen(t(trans))
	matrix_table$eigen_decomposition[[i]] = eigen
}



#----------------------------------------------------------


all = I(vector('list', 8))
first = I(vector('list', 8))
first_second = I(vector('list', 8))
first_second_third = I(vector('list', 8))
Q_first = I(vector('list', 8))
Q_first_second = I(vector('list', 8))
Q_first_second_third = I(vector('list', 8))
PCA_table = data.frame(matrices,all,first,first_second,first_second_third,Q_first,Q_first_second,Q_first_second_third)

for (i in 1:8){
	PCA_table$all[[i]] =c(matrix_table$Q_matrices[[i]])
	first = append(1,matrix_table$eigen_decomposition[[i]]$vectors[,1]/sum(matrix_table$eigen_decomposition[[i]]$vectors[,1]))
	second = append(matrix_table$eigen_decomposition[[i]]$values[2],matrix_table$eigen_decomposition[[i]]$vectors[,2])
	third = append(matrix_table$eigen_decomposition[[i]]$values[3],matrix_table$eigen_decomposition[[i]]$vectors[,3])
	Q_first = append(matrix_table$Q_eigen_decomposition[[i]]$values[1],matrix_table$Q_eigen_decomposition[[i]]$vectors[,1])
	Q_second = append(matrix_table$Q_eigen_decomposition[[i]]$values[2],matrix_table$Q_eigen_decomposition[[i]]$vectors[,2])
	Q_third = append(matrix_table$Q_eigen_decomposition[[i]]$values[3],matrix_table$Q_eigen_decomposition[[i]]$vectors[,3])

	PCA_table$first[[i]] = first
	PCA_table$first_second[[i]] = append(first,second)
	PCA_table$first_second_third[[i]] = append(PCA_table$first_second[[i]],third)
	PCA_table$Q_first[[i]] = Q_first
	PCA_table$Q_first_second[[i]] = append(Q_first,Q_second)
	PCA_table$Q_first_second_third[[i]] = append(PCA_table$Q_first_second[[i]],Q_third)
}



matrix_type = c("transition","Q")
all = I(vector('list', 2))
first = I(vector('list', 2))
first_second = I(vector('list', 2))
first_second_third = I(vector('list', 2))
Q_first = I(vector('list', 2))
Q_first_second = I(vector('list', 2))
Q_first_second_third = I(vector('list', 2))
PCA_results = data.frame(matrix_type,all,first,first_second,first_second_third,Q_first,Q_first_second,Q_first_second_third)

rownames = c("Q.bird","Q.insect","Q.LG","Q.mammal","Q.pfam","Q.plant","Q.yeast" )

PCA_data_all = do.call(rbind,PCA_table$all[19:25])
rownames(PCA_data_all) = rownames
PCA_results$all[[1]] = prcomp(t(PCA_data_all),scale=T)

PCA_data_first = do.call(rbind,PCA_table$first[19:25])
rownames(PCA_data_first) = rownames
PCA_results$first[[1]] = prcomp(t(PCA_data_first),scale=T)

PCA_data_first_second = do.call(rbind,PCA_table$first_second[19:25])
rownames(PCA_data_first_second) = rownames
PCA_results$first_second[[1]] = prcomp(t(PCA_data_first_second),scale=T)

PCA_data_first_second_third = do.call(rbind,PCA_table$first_second_third[19:25])
rownames(PCA_data_first_second_third) = rownames
PCA_results$first_second_third[[1]] = prcomp(t(PCA_data_first_second_third),scale=T)

Q_PCA_data_first = do.call(rbind,PCA_table$Q_first[19:25])
rownames(PCA_data_first) = rownames
PCA_results$Q_first[[1]] = prcomp(t(Q_PCA_data_first),scale=T)

Q_PCA_data_first_second = do.call(rbind,PCA_table$Q_first_second[19:25])
rownames(PCA_data_first_second) = rownames
PCA_results$Q_first_second[[1]] = prcomp(t(Q_PCA_data_first_second),scale=T)

Q_PCA_data_first_second_third = do.call(rbind,PCA_table$Q_first_second_third[19:25])
rownames(PCA_data_first_second_third) = rownames
PCA_results$Q_first_second_third[[1]] = prcomp(t(Q_PCA_data_first_second_third),scale=T)

NQ_rownames = c("NQ.bird","NQ.insect","NQ.mammal","NQ.pfam","NQ.plant","NQ.yeast" )

NQ_PCA_data_first = do.call(rbind,PCA_table$first[29:34])
rownames(NQ_PCA_data_first) = NQ_rownames
PCA_results$first[[2]] = prcomp(t(NQ_PCA_data_first),scale=T)

NQ_PCA_data_first_second = do.call(rbind,PCA_table$first_second[29:34])
rownames(NQ_PCA_data_first_second) = NQ_rownames
PCA_results$first_second[[2]] = prcomp(t(NQ_PCA_data_first_second),scale=T)

NQ_PCA_data_first_second_third = do.call(rbind,PCA_table$first_second_third[29:34])
rownames(NQ_PCA_data_first_second_third) = NQ_rownames
PCA_results$first_second_third[[2]] = prcomp(t(NQ_PCA_data_first_second_third),scale=T)

NQ_Q_PCA_data_first = do.call(rbind,PCA_table$Q_first[29:34])
rownames(NQ_PCA_data_first) = NQ_rownames
PCA_results$Q_first[[2]] = prcomp(t(NQ_Q_PCA_data_first),scale=T)

NQ_Q_PCA_data_first_second = do.call(rbind,PCA_table$Q_first_second[29:34])
rownames(NQ_PCA_data_first_second) = NQ_rownames
PCA_results$Q_first_second[[2]] = prcomp(t(NQ_Q_PCA_data_first_second),scale=T)

NQ_Q_PCA_data_first_second_third = do.call(rbind,PCA_table$Q_first_second_third[29:34])
rownames(NQ_PCA_data_first_second_third) = NQ_rownames
PCA_results$Q_first_second_third[[2]] = prcomp(t(NQ_Q_PCA_data_first_second_third),scale=T)


fviz_pca_ind(PCA_results$first_second[[1]],repel = TRUE)
colnames = c("A1","R1","N1","D1","C1","Q1","E1","G1","H1","I1","L1","K1","M1","F1","P1","S1","T1","W1","Y1","V1","A2","R2","N2","D2","C2","Q2","E2","G2","H2","I2","L2","K2","M2","F2","P2","S2","T2","W2","Y2","V2")






#---------------------------------------------------------
species = c("Q.mammal","Q.bird","Q.insect","Q.plant","Q.yeast","Q.mammal_cloak","Q.bird_cloak","Q.insect_cloak","Q.plant_cloak","Q.yeast_cloak","NQ.mammal","NQ.bird","NQ.insect","NQ.plant","NQ.yeast","NQ.mammal_cloak","NQ.bird_cloak","NQ.insect_cloak","NQ.plant_cloak","NQ.yeast_cloak")
species1 = c("Q.mammal","Q.bird","Q.insect","Q.plant","Q.yeast","Q.mammal_cloak","Q.bird_cloak","Q.insect_cloak","Q.plant_cloak","Q.yeast_cloak")
species2 = c("NQ.mammal","NQ.bird","NQ.insect","NQ.plant","NQ.yeast","NQ.mammal_cloak","NQ.bird_cloak","NQ.insect_cloak","NQ.plant_cloak","NQ.yeast_cloak")


species1 = c("Q.insect","Q.plant","Q.yeast","Q.mammal","Q.bird","Q.insect_cloak","Q.plant_cloak","Q.yeast_cloak","Q.mammal_cloak","Q.bird_cloak")

species = matrices
# species1 = c("Q.mammal","Q.pfam","Q.bac","Q.arch")
# species2 = c("NQ.mammal","NQ.pfam","NQ.bac","NQ.arch")
df_pca <- prcomp(allQ)
df_out <- as.data.frame(df_pca$x)
df_out$group = species
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(df_out,aes(x=PC1,y=PC2, color=group, label = species)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=30)) + 
  geom_point() +
  geom_text_repel(size=5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  xlim(-2,3) +
  ylim(-2,3)
ggsave(filename = "disorder.jpg", width=5, height=5, units="in")





species = c("Q.mammal","Q.pfam","Q.bac","Q.arch","NQ.mammal","NQ.pfam","NQ.bac","NQ.arch")

df_pca <- prcomp(allQ,rank. = 3)
df_out <- as.data.frame(df_pca$x)
df_out$group = species
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
print(percentage)

df_pca <- prcomp(allQ,rank. = 3)
components <- df_pca[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

labx="PC1: 19.63%"
laby="PC2: 12.14%"
labz="PC3: 10.74%"

axx <- list(
  nticks = 8,
  range = c(-2,2),
  title = labx
)

axy <- list(
  nticks = 8,
  range = c(-2,2),
  title = laby
)

axz <- list(
  nticks = 8,
  range = c(-2,2),
  title = labz
)

fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, type="scatter3d", mode="markers", text=species,color=species)
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig <- fig %>% layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=1)))
fig

 %>%
  add_markers(size = 12)
fig <- fig %>%
  layout(
    scene = list(bgcolor = "#e5ecf6")
)
fig


ggplot(df_out,aes(x=PC1,y=PC2, color=group, label = species)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=30)) + 
  geom_point() +
  geom_text_repel(size=7) +
  xlab(percentage[1]) +
  ylab(percentage[2])
ggsave(filename = "allQ.jpg", width=5, height=5, units="in")

df_pca <- prcomp(allQ)
df_out <- as.data.frame(df_pca$x)
df_out$group = species1
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(df_out,aes(x=PC1,y=PC2, color=group, label = species1)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=30)) + 
  geom_point() +
  geom_text_repel(size=7) +
  xlab(percentage[1]) +
  ylab(percentage[2])
ggsave(filename = "Q.jpg", width=5, height=5, units="in")

df_pca <- prcomp(allQ2)
df_out <- as.data.frame(df_pca$x)
df_out$group = species2
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(df_out,aes(x=PC1,y=PC2, color=group, label = species2)) + 
  theme(legend.position="none", legend.title = element_blank(),text = element_text(size=30)) + 
  geom_point() +
  geom_text_repel(size=7) +
  xlab(percentage[1]) +
  ylab(percentage[2])
ggsave(filename = "NQ.jpg", width=5, height=5, units="in")