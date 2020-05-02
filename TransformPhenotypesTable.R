#install.packages("xml2")

#Scrape this table:
library(rvest)
library(dplyr)
library(xml2)
library(magrittr)

library(RCurl)
library(stringr)
library(data.table)

#First load flower genotype to phenotype map.
#AKA maps genes to colors for all flowers
#This data grabbed from https://aiterusawato.github.io/guides/acnh/flowers.html
#Downloaded and tranformation is in file TransformPhenotypeTable.R and saved to ACNH_PHENOTYPES_TIDY.CSV
# RUN TransformPhenotypeTable.R FIRST to get latest data. File also in repo for ease of use.
genotypePage <- "https://aiterusawato.github.io/guides/acnh/genotypes.html"

pageLines <- download.file(genotypePage, "genotypes.html")

# # match rows like this:
# table tr:nth-child(7) td:nth-child(8),
# table tr:nth-child(8) td:nth-child(9),
# table tr:nth-child(7) td:nth-child(10),
# table tr:nth-child(7) td:nth-child(11),
#
colorCellMatches <- str_match(pageLines, "child\\((\\d+).*child\\((\\d+)\\)")
naturalFlowerCoords <-  colorCellMatches[!is.na(colorCellMatches[,1]),][,-1]
matchesDim <- dim(naturalFlowerCoords)
naturalFlowerCoords <- as.numeric(naturalFlowerCoords)
dim(naturalFlowerCoords) <- matchesDim
rm(matchesDim)
#[row, column]

pageLines <- readLines(genotypePage)

webpage <- getURL(genotypePage)
webpage <- readLines(tc <- textConnection("genotypes.html")); close(tc)

genotypeHtml <- read_html(genotypePage)

genotypeHtml
str(genotypeHtml)

header_row <- genotypeHtml %>% 
  html_node('thead tr') %>% 
  html_children()

columnNames <- header_row %>% html_node("span") %>% html_text()
#first three should be g1-g3
columnNames[4:6] <- paste0(columnNames[4:6], 0:2)

#then rose rose rose, and the other colors
data_rows <- genotypeHtml %>% 
  html_node('tbody') %>%
  html_children()

combinedDataFrame <- lapply(data_rows, function(datarow) {
  dataRow1Spans <- datarow %>% html_children() %>% html_children()  
  columnTypes <- dataRow1Spans %>% html_attr("class")
  columnText <- dataRow1Spans %>% html_text()
  geneColumns <- as.numeric(columnText[1:3])
  
  titleColumns <- dataRow1Spans[!is.na(columnTypes)] %>% html_children() %>% html_attr("title")
  
  rowdt <- do.call(data.frame, c(as.list(as.numeric(geneColumns)), as.list(titleColumns)))
  names(rowdt) <- columnNames
  rowdt
}) %>% rbindlist()

naturalGeneTable <- data.table(naturalFlowerCoords)
setnames(naturalGeneTable, c("row", "column"))
naturalGeneTable[,flowerName := as.matrix(combinedDataFrame)[naturalFlowerCoords]]
naturalGeneTable[,c("G1", "G2", "G3") := combinedDataFrame[row,.(G1,G2,G3)]]
tempG4 <- naturalGeneTable[, column-4]
tempG4[tempG4 > 2] <- 0
naturalGeneTable[, G4 := tempG4]
setkey(naturalGeneTable, G1, G2, G3, G4, flowerName)
#setkey(naturalGeneTable, flowerName)
#TODO: need to set row numbers here and compare against css code above.

meltedData <- melt(combinedDataFrame, id.vars=1:3)
meltedData[, G4 := 0]
meltedData[str_sub(variable,1,4) == "Rose", G4 := as.numeric(str_sub(variable,-1,-1))]
meltedData[str_sub(variable,1,4) == "Rose", variable := "Rose"]
meltedData[, GeneString := paste0(G1,G2,G3,G4)]
splitCol <- meltedData$value %>% str_split_fixed(" ", 2)
meltedData[, Color := tolower(splitCol[,1])]
meltedData[, IsNatural := FALSE]
setnames(meltedData, "variable", "Flower")
setnames(meltedData, "value", "flowerName")

setkey(meltedData, G1, G2, G3, G4, flowerName)
meltedData[naturalGeneTable, IsNatural := T]
meltedData[,flowerName := NULL]

setcolorder(meltedData, c("G1", "G2", "G3", "G4", "GeneString", "Flower", "Color", "IsNatural"))

setkey(meltedData, GeneString, Flower)
write.csv(meltedData, "ACNH_PHENOTYPES_TIDY.CSV", row.names = F)
