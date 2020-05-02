library(data.table)
library(stringr)
library(magrittr)
library(progress)

#First load flower genotype to phenotype map.
#AKA maps genes to colors for all flowers
#This data grabbed from https://aiterusawato.github.io/guides/acnh/flowers.html
#Downloaded and tranformation is in file TransformPhenotypeTable.R and saved to ACNH_PHENOTYPES_TIDY.CSV
# RUN TransformPhenotypeTable.R FIRST to get latest data. File also in repo for ease of use.
print("Loading phenotype file.")
allColors <- fread("ACNH_PHENOTYPES_TIDY.CSV")
setkey(allColors, G1,G2,G3,G4,Flower)

#Helper functions for handling genes.
#Genes basically have two representations, list and character
#Character each gene is 4 characters in trinary notation, each character being either 0, 1, or 2.
#list representation is a list of four two element character vectors, with the character vectors holding one of the below sub genes

#subgenes (gene value)
GV_0 <- c(0,0)
GV_1 <- c(0,1)
GV_2 <- c(1,1)

#genes of seed colors for roses
SEED_ROSE_WHITE <- list(GV_0, GV_0, GV_1, GV_0)
SEED_ROSE_RED   <- list(GV_2, GV_0, GV_0, GV_1)
SEED_ROSE_YELLOW <- list(GV_0, GV_2, GV_0, GV_0)

geneListToVec <- function(gene) {
  sapply(gene, sum)
}

geneListToStr <- function(gene) {
  gv <- geneListToVec(gene)
  paste0(gv, collapse='')
}

geneStrToList <- function(geneString) {
  lapply(str_split_fixed(geneString, "", 4), function(x) {
    get(paste0("GV_", x))
  })
}

#Helper functions to lookup the color of a gene
geneListToColor <- function(gene, flower = "Rose") {
  geneSum <- lapply(gene, sum)
  allColors[geneSum][Flower==flower, Color]
}
geneStrToColor <- function(gene, flower = "Rose") {
  sapply(lapply(gene, geneStrToList), geneListToColor, flower=flower)
}
geneToColor <- function(gene, flower = "Rose") {
  if(is.character(gene)) geneStrToColor(gene, flower) else geneListToColor(gene, flower)
}


#Function to breed two sets of flowers together and return resulting combinations
#parent 1 parent 2 - gene lists (list of two character vectors in binary notation valid char vectors are 00, 01, 11)
#rawOut -- return list of genes rathe than summary table
## Some examples of using breedFlowers
# breedFlowers(SEED_ROSE_WHITE)
# breedFlowers(SEED_ROSE_RED)
# 
# breedFlowers(SEED_ROSE_RED, SEED_ROSE_YELLOW)
# breedFlowers(SEED_ROSE_RED, SEED_ROSE_WHITE)
# 
# breedFlowers(SEED_ROSE_WHITE, SEED_ROSE_YELLOW)
breedFlowers <- function(parent1, parent2 = parent1, rawOut = FALSE) {
  #parent1 <- list(c(0,1), c(0,1), c(0,1), c(0,1))
  #parent2 <- list(c(0,1), c(0,1), c(0,1), c(0,1))
  
  if(is.character(parent1)) {
    parent1 <- geneStrToList(parent1)
  }
  
  if(is.character(parent2)) {
    parent2 <- geneStrToList(parent2)
  }
  
  if(geneListToStr(parent2) < geneListToStr(parent1)) {
    #swap
    tmp <- parent1
    parent1 <- parent2
    parent2 <- tmp
  }
  
  p1Gives <- do.call(expand.grid, parent1)
  p2Gives <- do.call(expand.grid, parent2)
  
  parentPicker <- expand.grid(p1 = seq(nrow(p2Gives)), p2 = seq(nrow(p2Gives)))
  
  finalGenes <- lapply(seq(nrow(parentPicker)), function(rn) {
    gene1 <- c(p1Gives[parentPicker[rn, "p1"],])
    gene2 <- c(p2Gives[parentPicker[rn, "p2"],])
    
    finalGene <- lapply(seq_along(gene1), function(gn) {
      sort(c(gene1[[gn]], gene2[[gn]]))
    })
    finalGene
  })
  finalGenes
  
  if(rawOut) return(finalGenes)
  
  geneOutcomes <- as.data.table(sapply(finalGenes, geneListToStr) %>% table() %>% divide_by(length(finalGenes)))
  geneOutcomes[,Parent1 := geneListToStr(parent1)]
  geneOutcomes[,Parent2 := geneListToStr(parent2)]
  setnames(geneOutcomes, c('N', '.'), c("Probability", "ChildGene"))
  setcolorder(geneOutcomes, c("Parent1", "Parent2", "ChildGene", "Probability"))
  
  geneOutcomes[, Color := geneStrToColor(ChildGene)]
  geneOutcomes[,ColorIsUnique := .N==1, by=Color]
  
  geneOutcomes[, NumOutputOfColor := .N, by=Color]
  geneOutcomes[, TotalColorProb := sum(Probability), by=Color]
  geneOutcomes[, ColorPickProb := Probability/TotalColorProb, by=Color]
  geneOutcomes[, AmbiguityWeightedProb := Probability * ColorPickProb, by=Color]
  
  geneOutcomes[]
}



##### PreCalculating Combination Results

#Get al possible gene strings, sorted
allPossibleGenes <- apply(expand.grid(rep(list(0:2), 4)), 1, paste0, collapse="") %>% sort
#Get all possible combinations
allPossibleBreeds <- data.table(expand.grid(Parent1 = allPossibleGenes, Parent2 = allPossibleGenes, stringsAsFactors = F))
#Now remove the mirror duplicates and standardize parent1 < parent2
allPossibleBreeds[,ParentL := ifelse(Parent1 < Parent2, Parent1, Parent2)]
allPossibleBreeds[,ParentH := ifelse(Parent1 < Parent2, Parent2, Parent1)]
allPossibleBreeds[,Parent1 := ParentL]
allPossibleBreeds[,Parent2 := ParentH]
allPossibleBreeds[,ParentL := NULL]
allPossibleBreeds[,ParentH := NULL]
allPossibleBreeds <- allPossibleBreeds[!duplicated(allPossibleBreeds)]

#Master list of all possible breed results
knownBreedList <- NULL

#Represent seeds in this list
addSeedGeneToBreedList <- function(gene) {
  geneStr <- geneListToStr(gene)
  geneRow <- data.table(Parent1 = "*SEED", Parent2 = "*SEED", ChildGene = geneStr, Probability = 1, Color = geneListToColor(gene), ColorIsUnique = TRUE,
                        NumOutputOfColor = 1, TotalColorProb = 1, ColorPickProb = 1, AmbiguityWeightedProb = 1)
  knownBreedList <<-  rbind(knownBreedList, geneRow)
  invisible(knownBreedList[])
}

knownBreedList <- NULL
addSeedGeneToBreedList(SEED_ROSE_WHITE)
addSeedGeneToBreedList(SEED_ROSE_RED)
addSeedGeneToBreedList(SEED_ROSE_YELLOW)
setkey(knownBreedList, Parent1, Parent2)

knownBreedList

#create a progress bar because this step takes a while.
pb <- progress_bar$new(
  format = " Calculating Breeds [:bar] :current/:total (:percent) eta: :eta",
  total = nrow(allPossibleBreeds), clear = FALSE)

#Loop over all possible combos and add to known breed list
for(breedRow in 1:nrow(allPossibleBreeds)) {
  pb$tick()
  #pick next pair
  breedPair <-   allPossibleBreeds[breedRow]
  #2. breed
  breedResults <- breedFlowers(breedPair[[1]], breedPair[[2]])
  #3. store
  knownBreedList <- rbind(knownBreedList, breedResults)
  if(any(duplicated(knownBreedList))) stop("fail dupes in breed list")
}

knownBreedList
setkey(knownBreedList, Parent1, Parent2)


##### Pathfinding (modified dijkstra) through flower breeds
#Step 1 is to create a shortest path tree through the breed graph.

pathfindTable <- data.table(gene = knownBreedList$ChildGene %>% unique() %>% sort())
setkey(pathfindTable, gene)

pathfindTable[,visited := FALSE] #not really needed
pathfindTable[, incomingCombo := rep(list(NULL), nrow(pathfindTable))] #two parents incoming to make this
pathfindTable[, minDist := Inf] #two parents incoming to make this
pathfindTable[, totalProbability := 0] #what is the multiplied probability of getting this? used as tiebreaker after minDist
pathfindTable[, expectedSpawns := Inf]
pathfindTable[, inSpt := FALSE] #in shortest path tree?
pathfindTable[, color := geneStrToColor(gene)]

pathfindPairTable <- knownBreedList[,1:2, with=F]
pathfindPairTable <- pathfindPairTable[!duplicated(pathfindPairTable)]
pathfindPairTable[,visited := FALSE]

#initialize starting points
addSeedToPathfindTable <- function(geneName) {
  pathfindTable[geneName, incomingCombo := .(list(c('*SEED', '*SEED')))]
  pathfindTable[geneName, minDist := 0]
  pathfindTable[geneName, visited := TRUE]
  pathfindTable[geneName, totalProbability := 1]
  pathfindTable[geneName, expectedSpawns := 1]
}

addSeedToPathfindTable("0010")
addSeedToPathfindTable("2001")
addSeedToPathfindTable("0200")


#Loopy Bits
pb <- progress_bar$new(
  format = " Calculating Shortest Path Tree [:bar] :current/:total (:percent) eta: :eta",
  total = length(allPossibleGenes), clear = FALSE)
for(i in seq_along(allPossibleGenes)) { 
  #This means some breeds aren't possible without unambiguous paths
  
  #grab next gene to expand from visited but not finished
  openList <- pathfindTable[visited==T & inSpt == FALSE][order(minDist, -totalProbability)]
  if(nrow(openList) == 0) break;
  
  nextGene <- openList[1,gene] #the full list without subset is the "open list"
  pathfindTable[nextGene, inSpt := TRUE]
  # print("Expanding Gene:")
  # print(pathfindTable[nextGene,])
  
  #all other breedable combos with unique colors
  #note: this filters non unique combos.
  #XXX CHANGE BETWEEN THIS AND THE COMMENTED LINE TO SWICH TO AMBIGUOUS VS UNAMBIGUOUS SEARCH
  allCombosWithGene <- knownBreedList[(Parent1 == nextGene | Parent2 == nextGene) & ChildGene != nextGene & ColorIsUnique == T] #unambig
  # allCombosWithGene <- knownBreedList[(Parent1 == nextGene | Parent2 == nextGene) & ChildGene != nextGene] #ambig
  
  
  allCombosWithGene[, otherGene := ifelse(Parent1 == nextGene, Parent2, Parent1)]
  #filter down to unique other genes we know how to make
  uniqueKnownOthers <- pathfindTable[unique(allCombosWithGene$otherGene)][visited == TRUE, gene]
  
  newCombos <- allCombosWithGene[otherGene %in% uniqueKnownOthers]
  #now update all of these new combos we just made
  for(ncr in seq(nrow(newCombos))) {
    comboData <- newCombos[ncr,]
    updateGene <- comboData$ChildGene
    pathfindTable[updateGene, visited := TRUE]
    
    nextDist <- max(pathfindTable[nextGene, minDist], pathfindTable[comboData$otherGene, minDist]) + 1
    nextProb <- min(pathfindTable[nextGene, totalProbability], pathfindTable[comboData$otherGene, totalProbability]) * comboData$AmbiguityWeightedProb
    nextExpectSpawns <- max(pathfindTable[nextGene, expectedSpawns], pathfindTable[comboData$otherGene, expectedSpawns]) + (1 / comboData$AmbiguityWeightedProb)^1
    #if(nextDist < pathfindTable[updateGene, minDist] || nextDist == pathfindTable[updateGene, minDist] && nextProb > pathfindTable[updateGene, totalProbability]) {
    if(nextExpectSpawns < pathfindTable[updateGene, expectedSpawns]) {
      pathfindTable[updateGene, minDist := nextDist]
      pathfindTable[updateGene, totalProbability := nextProb]
      pathfindTable[updateGene, expectedSpawns := nextExpectSpawns]
      pathfindTable[updateGene, incomingCombo := .(list(c(comboData$Parent1, comboData$Parent2)))]
    }
  }
  pb$tick()
}
cat("Finished building SPT.\nPaths found for", pathfindTable[inSpt == T, .N], "of", nrow(pathfindTable),"genes.") #awk syntax since cant just write !inSpt


# pathTableAmbiguous <- copy(pathfindTable)
# pathTableUnambiguous <- copy(pathfindTable)
# pathfindTable <- copy(pathTableAmbiguous)

# step 2 of pathfinding is to work backwards up the shortest path tree from your destination to the root(s)
#work backwards to produce path
#destinationGene <- '2220' #blue is 2220
findGene <- function(destinationGene, includeDetails = F) {
  #Find all involved genes
  if(is.numeric(destinationGene)) {
    destinationGene <- str_pad(destinationGene, 4, "left", "0")
  }
  closedParentSet <- c('*SEED', destinationGene)
  openParentSet <- destinationGene
  while(length(openParentSet) > 0) {
    openParentSet <- pathfindTable[openParentSet, incomingCombo] %>% unlist %>% unique
    openParentSet <- setdiff(openParentSet, closedParentSet)
    closedParentSet <- c(closedParentSet, openParentSet)  
  }
  involvedGenes <- setdiff(closedParentSet, "*SEED")
  rm(openParentSet); rm(closedParentSet)
  
  #Look up breed info from involved genes
  pathfindTable[involvedGenes, ]
  p1 <- sapply(pathfindTable[involvedGenes,incomingCombo], head, 1)
  p2 <- sapply(pathfindTable[involvedGenes,incomingCombo], tail, 1)
  
  setkey(knownBreedList, Parent1, Parent2, ChildGene)
  breedingInfo <- knownBreedList[J(p1, p2, involvedGenes)]
  breedingInfo <- merge(breedingInfo, pathfindTable[involvedGenes, .(gene, Order = minDist, totalProbability, expectedSpawns)], by.x = 'ChildGene', by.y='gene')
  breedingInfo[, P1Color := NA_character_]
  breedingInfo[, P2Color := NA_character_]
  breedingInfo[str_sub(Parent1,1,1) != '*', P1Color := geneStrToColor(Parent1)]
  breedingInfo[str_sub(Parent2,1,1) != '*',P2Color := geneStrToColor(Parent2)]
  breedingInfo <- breedingInfo[order(Order), .(Order, Parent1, Parent2, ChildGene, Probability = round(Probability, 3), WeightedProb = round(AmbiguityWeightedProb,3), ExpectedSpawns = expectedSpawns, P1Color, P2Color, ChildColor = Color)]
  
  stringOutput <- paste0("Breeding target: ", destinationGene, " - ", breedingInfo[ChildGene == destinationGene, ChildColor],"\n")
  #iterate over each breed and get table for it
  for(breedingRow in seq(nrow(breedingInfo))) {
    brInfo <- breedingInfo[breedingRow, ]
    if(brInfo$Parent1 == "*SEED") {
      stringOutput <- paste0(stringOutput,"  Phase ", brInfo$Order, " :  seeds -> ", brInfo$ChildGene, " ", str_pad(brInfo$ChildColor, 6, "right"),  "\n")
    } else {
      outcomeTable <- knownBreedList[Parent1 == brInfo$Parent1 & Parent2 == brInfo$Parent2]
      isUnique <- outcomeTable[Color == brInfo$ChildColor, .N] == 1
      stringOutput <- paste0(stringOutput,"  Phase ", brInfo$Order, " :  ", 
                             brInfo$Parent1, " ", str_pad(brInfo$P1Color, 6, "right"), " X ", 
                             brInfo$Parent2, " ", str_pad(brInfo$P2Color, 6, "right"),  " -> ", 
                             brInfo$ChildGene, " ", str_pad(brInfo$ChildColor, 6, "right"), 
                             str_pad(paste0("(", brInfo$Probability*100, "%) "), 10, "left") ,
                             "Total Avg. Spawns: ", brInfo$ExpectedSpawns, " ",
                             ifelse(isUnique, "Unique", paste0("AMBIGUOUS (Penalized Prob: ", brInfo$WeightedProb, "%)")), " \n")  
      
      if(includeDetails) {
        stringOutput <- paste0(stringOutput,"    Detailed Outcomes:\n")
        detailOutcome <- capture.output(print(outcomeTable))
        stringOutput <- paste0(stringOutput, paste0(sapply(detailOutcome, function(x) paste0("    ", x)), collapse="\n"), "\n")
      }
      # print(knownBreedList[Parent1 == brInfo$Parent1 & Parent2 == brInfo$Parent2])
      # cat("\n")
    }
  }
  stringOutput
}

cat(findGene("2220", F)) #blue


cat(findGene("1120"))
cat(findGene("1210"))
cat(findGene("2110"))

#### Trying to make igraph of stuff
library(igraph)
geneGraphData <- pathfindTable[order(expectedSpawns),]
setkey(geneGraphData, gene)
geneGraphData <- geneGraphData[visited==TRUE,.(gene, incomingCombo, minDist, totalProbability, expectedSpawns, color)]
geneGraphData[, leftParent := sapply(geneGraphData$incomingCombo, head, 1)]
geneGraphData[, rightParent := sapply(geneGraphData$incomingCombo, tail, 1)]

geneGraphData[, hexColor := rgb(t(col2rgb(color))/255)]

#create a subset for graphing. Here I'm picking the min numberr for all colors
minRowsForAllColors <- geneGraphData[,.(i=.I,c=color)][!duplicated(c),max(i)]
geneGraphData <- geneGraphData[1:minRowsForAllColors,]

# geneGraphData <- geneGraphData[J(involvedGenes)]

geneGraphData[, leftParentIdx := match(leftParent, gene)]
geneGraphData[, rightParentIdx := match(rightParent, gene)]

leftParentEdges <- c(t(geneGraphData[,cbind(leftParentIdx,.I)]))
leftParentEdges[is.na(leftParentEdges)] <- nrow(geneGraphData) + 1

rightParentEdges <- c(t(geneGraphData[,cbind(rightParentIdx,.I)]))
rightParentEdges[is.na(rightParentEdges)] <- nrow(geneGraphData) + 1

geneGraphData[color=='black', color:='grey25']
geneGraphData[color=='blue', color:='skyblue']

g <- make_empty_graph(n = nrow(geneGraphData) + 1)
g %<>% set_vertex_attr("gene", value = c(geneGraphData$gene, "SEED"))
g %<>% set_vertex_attr("color", value = c(geneGraphData$color, 'lightgrey'))


g %<>% add_edges(leftParentEdges)
g %<>% add_edges(rightParentEdges)


lay <- layout_as_tree(g)
lay <- layout_with_dh(g)
plot(g, layout=lay, vertex.label=c(geneGraphData$gene, "SEED"), vertex.size=30)

id <- tkplot(g, layout=lay, vertex.label=c(geneGraphData$gene, "SEED"), vertex.size=30)

#### #print out all breeds to a file.
bigString <- "All ACNH rose gene unambiguous paths starting from seeds.\nOutput from script authored by Michael Walker (mwalk10@gmail.com)\n\n"
for(geneNum in 1:nrow(pathfindTable)) {
  genePathInfo <- pathfindTable[geneNum, ]
  print(genePathInfo)
  if(genePathInfo$inSpt) {
    bigString <- paste0(bigString, findGene(genePathInfo$gene))
  } else {
    bigString <- paste0(bigString, "Breeding target: ",  genePathInfo$gene," ", genePathInfo$color,"\n  No unambiguous path found :(\n")
  }
  #cat(bigString)
}
write(bigString, "Output/ACNH All Flower Paths Unambiguous.txt")
