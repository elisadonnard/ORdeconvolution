#!/usr/bin/env Rscript
# Removing reads from secondary ORs based on known similar ORs from multimappers and gapped reads
###
knownpairs = scan("empirical_OR_pairs_filtered",what = character())
#####################################################################
library(ggplot2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
file = args[1]
cut = args[2]
sample = gsub("clean_ORs_|.umi.distributions.txt","",file)
cleanfile = paste0(sample,"_umiClean.txt")
input = read.table(file, col.names = c("cell","gene","umi","reads"))
#####################################################################
if (cut==1) {
	# Add umi and read sum per OR, and number of ORs expressed for each cell
        cleancells = input %>%
          add_count(cell, gene) %>%
          rename(umisum=n) %>%
	  select(c("gene","cell","umisum")) %>%
	  distinct()
        # prep final OR expression matrix for merging
        orexprs = as.data.frame(tidyr::spread(cleancells, cell, umisum, fill=0))
        rownames(orexprs) = orexprs$gene
        orexprs = t(orexprs[,-1])
        orexprs = tibble::rownames_to_column(as.data.frame(orexprs),var="cell")
} else {
	###
	# Add umi and read sum per OR, and number of ORs expressed for each cell
	cells = input %>% 
	  add_count(cell, gene) %>% 
	  rename(umisum=n) %>%
	  add_count(cell) %>% 
	  rename(umitotal=n) %>%
	  group_by(cell) %>%
	  mutate(ornum = n_distinct(gene), umifrtotal = umisum/umitotal) %>%
	  group_by(cell,gene) %>%
	  mutate(readsum = sum(reads))
	ggplot(cells[cells$ornum>1,], aes(umifrtotal)) +
	  geom_density() +
	  theme_bw()
	ggsave(paste0(sample,"_UMI_fraction_per_OR_distribution.pdf"))

	###
	# Select cells with >1 OR and add max umi sum and read sum per cell to calculate ratios, 
	# collapse table to one line per OR
	multior = cells %>%
	  filter(ornum>1) %>%
	  group_by(cell) %>%
	  mutate(umimax = max(umisum), umiratio = umisum/umimax, readmax = max(readsum), readratio = readsum/readmax) %>%
	  select(-c(umi,reads)) %>%
	  distinct()

	###
	# make a list of pairs ocurring in these cells
	multiwide = reshape2::dcast(splitstackshape::getanID(multior, 'cell'), cell~.id, value.var='gene')
	pairs = multiwide[is.na(multiwide$`3`),1:3]
	colnames(pairs) = c("cell","OR1","OR2")
	groups = multiwide[!is.na(multiwide$`3`),]
	for (i in 1:nrow(groups)) {
	  cellname = groups$cell[i]
	  combo = data.frame(t(combn(groups[i,-1], 2)))
	  newpairs = as.data.frame(cbind(unlist(combo[,1]), unlist(combo[,2])))
	  colnames(newpairs) = c("OR1","OR2")
	  newpairs = newpairs[!is.na(newpairs$OR2),]
	  newpairs$cell = cellname
	  newpairs = newpairs[,c(3,1,2)]
	  pairs = rbind(pairs, newpairs)
	}
	rowsorted = data.frame(t(apply(pairs[,-1], 1, sort)))
	pairs$OR1 = rowsorted$X1
	pairs$OR2 = rowsorted$X2
	pairs$pair = paste(pairs$OR1,pairs$OR2,sep = "_")
	pairs$known = "no"
	pairs$known[pairs$pair%in%knownpairs] = "yes"
	tomerge = multior[,1:3]
	colnames(tomerge) = c("cell","OR1","OR1_umisum")
	pairs = merge(pairs, tomerge, by=c("cell","OR1"))
	colnames(tomerge) = c("cell","OR2","OR2_umisum")
	pairs = merge(pairs, tomerge, by=c("cell","OR2"))
	rm1 = pairs[pairs$known=="yes" & pairs$OR1_umisum<pairs$OR2_umisum,c("cell","OR1")]
	colnames(rm1) = c("cell","gene")
	rm2 = pairs[pairs$known=="yes" & pairs$OR2_umisum<pairs$OR1_umisum,c("cell","OR2")]
	colnames(rm2) = c("cell","gene")
	toremove = rbind(rm1,rm2)
	toremove$remove = "yes"
	cleaninput = anti_join(input,toremove,by=c("cell","gene"))
	# which cells will be affected by this filter?
	altered = as.character(toremove$cell)

	###
	# Add umi and read sum per OR, and number of ORs expressed for each cell
	cleancells = cleaninput %>% 
	  add_count(cell, gene) %>% 
	  rename(umisum=n) %>%
	  add_count(cell) %>% 
	  rename(umitotal=n) %>%
	  group_by(cell) %>%
	  mutate(umifrtotal = umisum/umitotal, maxfrcell = max(umifrtotal)) %>%
	  mutate(ornum = n_distinct(gene)) %>%
	  group_by(cell,gene) %>%
	  mutate(readsum = sum(reads))
	# keep one OR for cells with umifrtotal>=cut and all ors for cells where maxfrcell<cut
	cleanoutput = cleancells[cleancells$umifrtotal >= cut | cleancells$maxfrcell < cut,
	                         c("gene","cell","umisum")] %>% distinct()
	# which cells will be affected by this filter?
	altered = c(altered, unique(as.character(cleancells$cell[cleancells$umifrtotal >= cut & cleancells$ornum>1])))
	data.table::fwrite(as.data.frame(altered), paste0("cellnames_altered_",sample), sep = "\t")
	# prep final OR expression matrix for merging
	orexprs = as.data.frame(tidyr::spread(cleanoutput, cell, umisum, fill=0))
	rownames(orexprs) = orexprs$gene
	orexprs = t(orexprs[,-1])
	orexprs = tibble::rownames_to_column(as.data.frame(orexprs),var="cell")
}
###
# Read _umiClean.txt files and replace OR counts with new counts
exprsinput = read.table(cleanfile, stringsAsFactors = F, header = T, row.names = 1)
exprsinput = t(exprsinput[grep("Olfr", rownames(exprsinput), invert = T),])
exprsinput = tibble::rownames_to_column(as.data.frame(exprsinput),var="cell")
final = full_join(exprsinput,orexprs,by="cell")
rownames(final) = final$cell
final = final[,-1]
final[is.na(final)] = 0
final = t(final)
data.table::fwrite(as.data.frame(final), file=paste0("fixedOR_",cleanfile), sep = "\t", row.names = T)
