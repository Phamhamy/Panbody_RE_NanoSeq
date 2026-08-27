options(stringsAsFactors = F)
library(hdp)
library("RColorBrewer")
lower_threshold=0

chlist <- vector("list", 20)
for (i in 1:20){
  if(file.exists(paste0("hdp_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

mut_example_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi) #This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi,"HDP_multi_chain.Rdata")

pdf("muts_attributed.pdf")
plot_comp_size(mut_example_multi, bty="L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

#dev.new(width=12,height=4)
#par(mfrow=c(3,4))


for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_",i,".pdf"),width=12,height=4)
  
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

#plot_dp_comp_exposure(mut_example_multi,
 #                     dpindices=2:4, incl_numdata_plot=FALSE,
  #                    col=mycolors,
   #                   incl_nonsig=TRUE, cex.names=0.8,
    #                  ylab_exp = 'Signature exposure', leg.title = 'Signature')

#pdf("signature_attribution.pdf",width=10,height=8)
#plot_dp_comp_exposure(mut_example_multi, dpindices=(length(freq)+2):length(mut_example_multi@comp_dp_counts), incl_nonsig = T, ylab_exp = 'Signature exposure', leg.title = 'Signature', col=mycolors, incl_numdata_plot=F)
#dev.off()


dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
#mean_assignment <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),drop=FALSE])

mean_assignment=as.data.frame(comp_dp_distn(mut_example_multi)$mean)
write.table(mean_assignment,"mean_assignment_hdp.txt")
mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
write.table(mean_sigs,"hdp_sigs.txt")
