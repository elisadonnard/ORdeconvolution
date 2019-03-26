library(ggplot2)
# multimapper or gapped reads that align to more than one OR
multi = read.table("multimapped_and_gap_readnames_ORnames", 
                 col.names = c("id","read"))
# reference table with OR genename and transcript ID
ref = read.table("reftable_OR_name", 
                 col.names = c("id","OR"))
###
multi = merge(multi,ref,by="id")
multi = multi[,-1]
multi = multi[!duplicated(multi),]
multiwide = reshape2::dcast(splitstackshape::getanID(multi, 'read'), read~.id, value.var='OR')
multiwide = multiwide[!is.na(multiwide$`2`),]
multiwide = multiwide[,-1]
pairs = multiwide[is.na(multiwide$`3`),1:2]
colnames(pairs) = c("OR1","OR2")
groups = multiwide[!is.na(multiwide$`3`),]
for (i in 1:nrow(groups)) {
  combo = data.frame(t(combn(groups[i,], 2)))
  newpairs = as.data.frame(cbind(unlist(combo[,1]), unlist(combo[,2])))
  colnames(newpairs) = c("OR1","OR2")
  newpairs = newpairs[!is.na(newpairs$OR2),]
  pairs = rbind(pairs, newpairs)
}
pairs = data.frame(t(apply(pairs, 1, sort)))
colnames(pairs) = c("OR1","OR2")
counts = as.data.frame(table(pairs))
counts = counts[counts$Freq>1,]
ggplot(counts, aes(Freq)) +
  geom_density() +
  theme_bw() +
  scale_x_log10()
ggsave("ORpair_frequency_distribution.pdf")
over1 = paste(counts$OR1,counts$OR2, sep = "_")
pairs$pair = paste(pairs$OR1,pairs$OR2,sep = "_")
pairs = pairs[pairs$pair%in%over1,]
write.table(pairs, file = "empirical_OR_pairs", row.names = F, col.names = F, quote = F)
write.table(over1, file = "empirical_OR_pairs_filtered", row.names = F, col.names = F, quote = F)
