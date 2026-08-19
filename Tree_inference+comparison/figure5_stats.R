figure5 <- read.csv("figure5.csv")

figure5_filtered=filter(figure5,variable=='cleaned')
figure5_unfiltered=filter(figure5,variable=='uncleaned')
unfiltered1=filter(figure5_unfiltered,genes=="1")
unfiltered4=filter(figure5_unfiltered,genes=="4")
unfiltered8=filter(figure5_unfiltered,genes=="8")

cor.test(unfiltered1$distance, unfiltered1$length, method="spearman")
cor.test(unfiltered4$distance, unfiltered4$length, method="pearson")
cor.test(unfiltered8$distance, unfiltered8$length, method="pearson")


t.test(unfiltered1$distance, unfiltered4$distance)
t.test(unfiltered4$distance, unfiltered8$distance)