library(ggplot2)
library(ggrepel)

#Figure 1

figure1 <- read.csv("figure1.csv")

ggplot(figure1,aes(x=1-false_discovery_rate,y=true_positive_rate,label=labels),colour=types)+
	geom_point(size=4,colour=figure1$types)+labs(x="Precision",y="Recall")+
	geom_text_repel(size=7,box.padding=1,min.segment.length = unit(0, 'lines'),segment.size = 1) +
	theme(axis.title=element_text(size=22, face ="bold"),title=element_text(size=22, face ="bold"),axis.text = element_text(size = 18),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))+geom_hline(yintercept=0.7698268362364417, linetype="dashed", color = "blue", linewidth=0.5)+geom_vline(xintercept=1-0.21938855618981626, linetype="dashed", color = "blue", linewidth=0.5)+
	scale_y_continuous(limits = c(0,1), expand = c(0, 0))+scale_x_continuous(limits=c(0.7,1),expand=c(0,0)) +
	theme(plot.margin = margin(0.5,1,0,0, "cm"))

#Figure2









#Figure 3A

lrm1 <- read.csv("figure3a.csv")

ggplot(lrm1, aes(x=-value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("magenta","dark green","gold","cyan2")) +
    geom_density(size = 1) + labs(x="Change in Gene Tree/\n Species Tree Distance") + 
    geom_vline(aes(xintercept=17.3),color="cyan2", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=0),color="black", linetype="solid", size=1) +
	xlim(-40,60) 
	ylim(0,0.1)

#Figure 3B

lrm2 <- read.csv("figure3b.csv")

ggplot(lrm2, aes(x=-value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("magenta","dark green","gold","cyan2","grey")) +
    geom_density(size = 1) + labs(x="Change in Gene Tree/\n Species Tree Distance") + 
    geom_vline(aes(xintercept=17.3),color="black", linetype="solid", size=1) +
	geom_vline(aes(xintercept=20.8),color="magenta", linetype="dashed", size=1) +
    xlim(-40,60) +
	ylim(0,0.1)

#Figure 3C

lrm_gtrpmix <- read.csv("figure3c.csv")

ggplot(lrm_gtrpmix, aes(x=-value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("red")) +
    geom_density(size = 1) + labs(x="Change in Gene Tree/\n Species Tree Distance") + 
    geom_vline(aes(xintercept=22.2),color="black", linetype="solid", size=1) +
	xlim(-40,60) +
	ylim(0,0.1)

#Figure 3D

quartet1 <- read.csv("figure3d.csv")

ggplot(quartet1, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("magenta","dark green","gold","cyan2")) +
    geom_density(size = 1) + labs(x="Change in Quartet Support") + 
    geom_vline(aes(xintercept=0.057),color="cyan2", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=0),color="black", linetype="solid", size=1) +
    scale_x_continuous(breaks=c(-0.1,-0.05,0,0.05,0.1,0.15,0.2),limits=c(-0.15,0.25)) +
    ylim(0,20)

#Figure 3E

quartet2 <- read.csv("figure3e.csv")

ggplot(quartet2, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("magenta","dark green","gold","cyan2","grey")) +
    geom_density(size = 1) + labs(x="Change in Quartet Support") + 
    geom_vline(aes(xintercept=0.057),color="black", linetype="solid", size=1) +
    geom_vline(aes(xintercept=0.076),color="magenta", linetype="dashed", size=1) +
    scale_x_continuous(breaks=c(-0.1,-0.05,0,0.05,0.1,0.15,0.2),limits=c(-0.15,0.25)) +
    ylim(0,20)

#Figure 3F

quartet_gtrpmix <- read.csv("figure3f.csv")

ggplot(quartet_gtrpmix, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("red")) +
    geom_density(size = 1) + labs(x="Change in Gene Tree/\n Species Tree Distance") + 
    geom_vline(aes(xintercept=0.097),color="black", linetype="solid", size=1) +
    scale_x_continuous(breaks=c(-0.1,-0.05,0,0.05,0.1,0.15,0.2),limits=c(-0.15,0.25)) +
    ylim(0,30)

#Figure 3G

dlt1 <- read.csv("figure3g.csv")

ggplot(dlt1, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("magenta","dark green","gold","cyan2")) +
    geom_density(size = 1) + labs(x="Decrease in DLT/#Nodes") + 
    geom_vline(aes(xintercept=0.063),color="cyan2", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=0),color="black", linetype="solid", size=1) +
    xlim(-0.4,0.4) +
    ylim(0,10)

#Figure 3H

dlt2 <- read.csv("figure3h.csv")

ggplot(dlt2, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("magenta","dark green","gold","cyan2","grey")) +
    geom_density(size = 1) + labs(x="Decrease in DLT/#Nodes") + 
    geom_vline(aes(xintercept=0.047),color="black", linetype="solid", size=1) +
    geom_vline(aes(xintercept=0.063),color="magenta", linetype="dashed", size=1)+
    xlim(-0.4,0.4) +
    ylim(0,15)

#Figure 3I

dlt3 <- read.csv("figure3i.csv")

ggplot(dlt_gtrpmix, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=12))+
	scale_color_manual(values=c("red")) +
    geom_density(size = 1) + labs(x="Change in Gene Tree/\n Species Tree Distance") + 
    geom_vline(aes(xintercept=0.063),color="black", linetype="solid", size=1)+
    xlim(-0.4,0.4) +
    ylim(0,15)

#Figure 4A

ggplot(lrm_absolute, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=22))+
	scale_color_manual(values=c("magenta","dark green")) +
    geom_density(size = 1) + labs(x="LRM Tree Distance") +
    scale_x_reverse()

#Figure 4B

ggplot(quartet_absolute, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=22))+
	scale_color_manual(values=c("magenta","dark green")) +
    geom_density(size = 1) + labs(x="Quartet Support") +
    scale_x_continuous(breaks=c(0.4,0.6,0.8,1),limits=c(0.3,1))

#Figure 4C

ggplot(dlt_absolute, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=18))+
	scale_color_manual(values=c("magenta","dark green")) +
    geom_density(size = 1) + labs(x="DLT/#Tips") +
    scale_x_reverse(limits = c(5, 0), breaks = c(4,3,2,1,0)) 

#Figure 4D

ggplot(boot_df, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=22))+
	scale_color_manual(values=c("magenta","dark green")) +
    geom_density(size = 1) + labs(x="Mean Boostrap") + 
	scale_x_continuous(breaks=c(0,25,50,75,100),limits=c(0,100))

#Figure 4E

ggplot(scf, aes(x=value, color=variable))+
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(color = "black", size=22))+
	scale_color_manual(values=c("magenta","dark green")) +
    geom_density(size = 0.75) + labs(x="Mean sCF") + 
	scale_x_continuous(breaks=c(0,25,50,75,100),limits=c(0,100)) +
	ylim(0,0.05)

#Figure 5



#Supplementary Figure 1



#Supplementary Figure 2

