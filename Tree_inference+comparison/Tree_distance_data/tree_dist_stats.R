


lrm <- read.csv("lrm_distance.csv")
quartet <- read.csv("quartet.csv")
dlt <- read.csv("dlt.csv")


# Lin-Rajan-Moret t tests
t.test(lrm$unfiltered,lrm$cloak,paired=TRUE)
t.test(lrm$unfiltered,lrm$divvier,paired=TRUE)
t.test(lrm$unfiltered,lrm$taper,paired=TRUE)
t.test(lrm$unfiltered,lrm$divvierpf,paired=TRUE)
t.test(lrm$divvierpf,lrm$gene_cloak,paired=TRUE)
t.test(lrm$divvierpf,lrm$gene_divvier,paired=TRUE)
t.test(lrm$divvierpf,lrm$gene_taper,paired=TRUE)
t.test(lrm$divvierpf,lrm$gene_divvierpf,paired=TRUE)
t.test(lrm$divvierpf,lrm$propogate_error,paired=TRUE)
t.test(lrm$gene_cloak,lrm$gtrpmix,paired=TRUE)

# quartet support t tests
t.test(quartet$uncleaned,quartet$cloak,paired=TRUE)
t.test(quartet$uncleaned,quartet$divvier,paired=TRUE)
t.test(quartet$uncleaned,quartet$taper,paired=TRUE)
t.test(quartet$uncleaned,quartet$divvierpf,paired=TRUE)
t.test(quartet$divvierpf,quartet$gene_cloak,paired=TRUE)
t.test(quartet$divvierpf,quartet$gene_divvier,paired=TRUE)
t.test(quartet$divvierpf,quartet$gene_taper,paired=TRUE)
t.test(quartet$divvierpf,quartet$gene_divvierpf,paired=TRUE)
t.test(quartet$divvierpf,quartet$propagate,paired=TRUE)
t.test(quartet$gene_cloak,quartet$gtrpmix,paired=TRUE)

# Duplication/Loss/Transfer t-tests
t.test(dlt$uncleaned,dlt$cloak,paired=TRUE)
t.test(dlt$uncleaned,dlt$divvier,paired=TRUE)
t.test(dlt$uncleaned,dlt$taper,paired=TRUE)
t.test(dlt$uncleaned,dlt$divvierpf,paired=TRUE)
t.test(dlt$divvierpf,dlt$gene_cloak,paired=TRUE)
t.test(dlt$divvierpf,dlt$gene_divvier,paired=TRUE)
t.test(dlt$divvierpf,dlt$gene_taper,paired=TRUE)
t.test(dlt$divvierpf,dlt$gene_divvierpf,paired=TRUE)
t.test(dlt$divvierpf,dlt$propagate,paired=TRUE)
t.test(dlt$gene_cloak,dlt$gtrpmix,paired=TRUE)