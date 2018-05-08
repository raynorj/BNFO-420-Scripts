GetAllGenes <- function(filepath) {
	geneSet = vector()
	
	#for each GO line
	for (fline in readLines(filepath)) {
		fline = unlist(strsplit(fline, "\t"))
		lastIndex = length(fline)

		if (length(fline) > 0) {
			if (startsWith(fline[1], "GO:") == TRUE) {
				print("hello")
				#split annotated element into list of genes
				geneVec = unlist(strsplit(fline[lastIndex], ","))
				
				#trim whitespace
				for (i in 1:length(geneVec)) {
					geneVec[i] = trimws(geneVec[i])
				}
			
				#add GO term to vector
				geneSet = union(geneSet, geneVec)
			}
		}
	}
	
	return(geneSet)
}

GetGOGenes <- function(filepath, go_term) {
	goVec = vector()
	
	#for each GO line
	for (fline in readLines(filepath)) {
		fline = unlist(strsplit(fline, "\t"))
		lastIndex = length(fline)
		
		#if this is the GO term we're looking for, store gene vector
		if (length(fline) > 0) {
			#if (fline[1] == go_term) {
			if (startsWith(fline[1], go_term) == TRUE) {
				#split annotated element into list of genes
				geneVec = unlist(strsplit(fline[lastIndex], ","))
				
				#trim whitespace
				for (i in 1:length(geneVec)) {
					geneVec[i] = trimws(geneVec[i])
				}
				
				#add genes to GO term vector
				goVec = c(goVec, geneVec)
				break
			}
		}
	}
	
	return(goVec)
}

GetRankings <- function(evcTable) {
	#sort names by eigenvector centrality values in increasing order
	evcTable = evcTable[order(evcTable$V1),]
	#return sorted names
	return(as.character(evcTable$V2))
}

ConstructTable <- function(geneSet, goSet, rankedNames, evcTable) {
	membership = vector()
	ranks = vector()
	evcTable = evcTable[order(evcTable$V1),]
	
	for (i in 1:length(geneSet)) {
		#add gene membership to membership vector
		if (geneSet[i] %in% goSet) {
			membership = c(membership, 1)
		} else {
			membership = c(membership, 0)
		}
		
		#add ranked eigenvector centrality to evc
		ranks = c(ranks, match(geneSet[i], rankedNames))
	}
	
	geneTable = data.frame(geneSet, membership, ranks)
	colnames(geneTable) <- c("Gene", "GOTerm", "Rank")
	geneTable = geneTable[order(geneTable$Rank),]
	geneTable$Centrality = evcTable$V1

	return(geneTable)
}

RankBiserialCorrelation <- function(geneTable) {
	#filter out go term members
	members = subset(geneTable, GOTerm == 1)
	#filter out non-members
	nonmembers = subset(geneTable, GOTerm == 0)
	
	#calculate mean of member ranks
	Y1 = mean(members$Rank)
	#calculate mean of non-member ranks
	Y2 = mean(nonmembers$Rank)
	#get number of genes in table
	n = length(members$Rank) + length(nonmembers$Rank)
	
	#calculate 2 * (member mean - non-member mean) / # of genes
	rb = 2 * (Y1 - Y2) / n

	return(rb)
}

###From https://sebastiansauer.github.io/convert_logit2prob/
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
############################################################

month = "m2"

#specify GO term of interest
GOTerms = c("GO:0023052","GO:0022610","GO:0008283","GO:0050896","GO:0002376","GO:XXXXXXX")
#get all genes in gene enrichment
geneList = as.character(read.table(paste(month, "_evc_values.txt", sep=""), sep="\t")$V2)
#get gene names in rank-order from eigenvector centrality file
rankedNames = GetRankings(read.table(paste(month, "_evc_values.txt", sep=""), sep="\t"))
evcTable = read.table(paste(month, "_evc_values.txt", sep=""), sep="\t")
counter = 0

for (GOTerm in GOTerms) {
	GOList = GetGOGenes("query_tab_amended.txt", GOTerm)
	geneTable = ConstructTable(geneList, GOList, rankedNames, evcTable)
	
	GOModel <- glm(GOTerm ~ Centrality, family=binomial(link="logit"), data=geneTable)
	geneTable$Probability = predict(GOModel, newdata=geneTable, terms="response")

	for (i in 1:length(geneTable$Probability)) {
	  geneTable$Probability[i] = logit2prob(geneTable$Probability[i])
	}
	
	correlation = cor(geneTable$Centrality, geneTable$Probability, method="spearman")
  GeneLine <- lm(Probability ~ Centrality, data=geneTable)
	
	png(filename=paste(month, "_", counter, ".png", sep=""))
	plot(geneTable$Centrality, geneTable$Probability, xlab="Eigenvector Centrality", ylab="Probability", xlim = c(0.0, 1.0), ylim = c(0.0, 1.0), main = paste(GOTerm, " - Month 2", sep=""))
  text(x=0.3, y=0.9, paste("rho = ", correlation, sep=""))
  lines(geneTable$Centrality, fitted(GeneLine))
  dev.off()
  
	print(summary(GOModel))
	
	counter = counter + 1
}