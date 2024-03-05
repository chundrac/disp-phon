require(phytools)
require(ggtree)
require(ggplot2)
require(Cairo)

load('data_for_branch_vis.Rdata')

p <- ggtree(tree) + geom_tiplab(size=2,color='slategray',alpha=1)
ordered.taxa <- get_taxa_name(p)

ordered.inds <- c(1:length(tree$tip.label))
names(ordered.inds) <- ordered.taxa

tree$tip.label <- paste('  ',ordered.inds[tree$tip.label],'. ',tree$tip.label,sep='')

freqs$language <- tree$tip.label

#CairoPDF('branch_level_trends.pdf',width=7,height=7)
#pdf('branch_level_trends.pdf')
tikzDevice::tikz('branch_level_trends_final.tex',width=7,height=7)
ggtree(tree,aes(col=exp(rhos.med)),size=1) %<+% freqs + geom_tippoint(aes(size=inventory.size,alpha=inventory.size),color='black') +scale_size(name='Inv.\\ size',range=c(.01,1))  +
  coord_cartesian(clip='off') + 
  theme(plot.margin=margin(1, 260, 1, 1),legend.text.align = 0,legend.key.size = unit(.4, 'cm'),legend.key.height = unit(.4, 'cm'),legend.key.width = unit(.4, 'cm'),legend.position=c(.15,.775),legend.text=element_text(size=5))+
  scale_colour_gradient(name='Dist.\\ from origin', low='blue', high = 'yellow') + 
  geom_tiplab(size=2,color='slategray',alpha=1)

dev.off()

lines <- readLines(con="branch_level_trends.tex")
lines <- lines[-which(grepl("\\path\\[clip\\]*", lines,perl=F))]
lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines,perl=F))]
writeLines(lines,con="branch_level_trends.tex")

