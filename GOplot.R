library(GOplot)

data(EC)
str(EC)
head(EC$david)

circ <- circle_dat(EC$david, EC$genelist)
head(circ)
#Generate a simple barplot
par(mfrow = c(1,2))
GOBar(subset(circ, category == 'BP'))

# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')

# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

# Generate the bubble plot with a label threshold of 3
GOBubble(circ, labels = 3)

# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)

# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3) 

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)

#Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)

chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)

# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = EC$genes)
# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = EC$process)

# Create the plot
chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(3, 0), gene.order = 'logFC')

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-8], nlfc = 0)

GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)

l1 <- subset(circ, term == 'heart development', c(genes,logFC))
l2 <- subset(circ, term == 'plasma membrane', c(genes,logFC))
l3 <- subset(circ, term == 'tissue morphogenesis', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('heart development', 'plasma membrane', 'tissue morphogenesis'))
