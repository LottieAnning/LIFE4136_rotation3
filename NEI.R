#Create a matrix of pairwise distances for **individuals**:
aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE)

#Create a matrix of pairwise distances for **populations**:
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)

#Create the distance objects:
colnames(aa.D.ind) <- rownames(aa.D.ind)
aa.D.ind.dist <-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels") <-rownames(aa.D.ind)

colnames(aa.D.pop) <- rownames(aa.D.pop) 
aa.D.pop.dist <-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels") <-rownames(aa.D.pop)

#Plot and save the Neighbour Joining (NJ) tree:
plot(nj(aa.D.ind), typ="unrooted", cex=0.7)
title(expression("Neighbour-joining tree of distance-based analysis of "*italic(Arabidposis)*" "))
write.tree(nj(aa.D.pop),file="NJ.distance_tree_outgroups.tre")
