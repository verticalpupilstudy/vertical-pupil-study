#this script builds a time-scaled MCC consensus from 1000 sampled trees from Upham posterior


if (!requireNamespace("ape", quietly=TRUE))    install.packages("ape")
if (!requireNamespace("phangorn", quietly=TRUE)) install.packages("phangorn")

library(ape)
library(phangorn)

posterior_all <- read.nexus(
  "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_nexus.trees"
)

idx1000   <- seq(1, length(posterior_all), by = 10)[1:1000]
posterior1k <- posterior_all[idx1000]

mcc1000 <- maxCladeCred(posterior1k, rooted = FALSE)

out_tree_file <- "Upham2019_MCC_1000.tre"
write.tree(mcc1000, file = out_tree_file)


out_pdf <- "Upham2019_MCC_1000.pdf"
pdf(out_pdf, width = 10, height = 6)
plot(mcc1000,
     show.tip.label = FALSE,      
     use.edge.length = TRUE,       
     x.lim = c(0, max(branching.times(mcc1000)))
)
axisPhylo()                     
title("Upham et al. (2019) MCC Chronogram (1 000‑tree subset)")
dev.off()
