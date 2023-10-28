
library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)
library(compiler)
library(gsubfn)

# htn 'C0020538', 'C0085580'

#########
# DB setup
#########
drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "causalehr",
    host = "localhost",
    port = 5432,
    user = "scottalexandermalec", # "postgres",
    password = "mandarin"
  )

#########
# manually defined exclusions
#########
manuallyDefinedExclusions <- " ('C0450442', 'C0243077', 'C0021083', 'C0013216', 'C0009429', 'C0199176', 'C0085104',
                               'C0017296', 'C0597357', 'C1519033', 'C1280903', 'C0003289', 'C0279025', 'C0279025',
                               'C0040808', 'C0009244', 'C0175723', 'C0011581', 'C0026868', 'C0683465', 'C1519595',
                               'C0042890', 'C0150321', 'C0581602', 'C1874451', 'C0012854', 'C0033972', 'C0279493',
                               'C0014935', 'C0086132', 'C0344315', 'C0016452', 'C1137094', 'C2718059', 'C0013806',
                               'C0183185', 'C0678908', 'C2827774', 'C0031765', 'C0003827', 'C0034991', 'C0581601',
                               'C0282402', 'C0150593', 'C1254359', 'C3858690', 'C1096024', 'C0020431', 'C1293131',
                               'C1964029', 'C0204727', 'C0005893', 'C0185027', 'C0439775', 'C0418981', 'C0152060',
                               'C0183210', 'C0870230', 'C0677850', 'C0376626', 'C0677850', 'C0949766', 'C0150133',
                               'C1515119', 'C0079613', 'C0455142', 'C0585941', 'C0279494', 'C0920425', 'C0033968', 
                               'C0935576', 'C1527374', 'C1269683', 'C0034618', 'C0376547', 'C0020431', 'C0185117',
                               'C0204514', 'C0439857', 'C0680038', 'C0846672', 'C0015618', 'C0182537', 'C0183115', 
                               'C0220908', 'C0279494', 'C0935576', 'C1268930', 'C0013218', 'C0034764', 'C0183210', 
                               'C0232164', 'C0281481', 'C1139730', 'C1511300', 'C1512629', 'C2986605', 'C3177188',
                               'C0030567', 'C0497327', 'C0474169', 'C0011269', 'C0599917', 'C0599917',
                               'C0290627', 'C0243072', 'C0042890', 'C0026336') "

length(manuallyDefinedExclusions)

#########
# Covariates that are measured for subjects in sample from electronic health database
#########
getGoodCovars <- function(n) {
  fNames <- list.files(path = "data/r3dataset/data")
  i = 0
  goodData <- c()
  for (f in fNames) {
    if (sum(read.csv(file = paste("data/r3dataset/data/", f, sep = ""), header = FALSE)) >= n) 
    { #print(f) 
      i = i + 1 
      #print(i)
      goodData <- c(goodData, f) 
    }
  }
  #print(goodData)
  goodData <- gsubfn(goodData, pattern = "\\.csv", replacement = "")
  return(goodData)
}

predicates <- " ('CAUSES') " #, 'PREDISPOSES') " #, 'PREVENTS', 'PREDISPOSES') " 
# NOTE: 'TREATS' -- excluding because of the ambiguous causal meaning of the verb 'to treat'
#predicates <- " ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'AFFECTS', 'DISRUPTS', 'AUGMENTS') "

#########
# getConfounders for each dementia subtype
#########
#
getConfounders <- function(predicates, DEMsqldef) {
  confounders <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.subject_cui AS covar, COUNT(*) AS scnt, cp.subject_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND cp.subject_cui NOT IN ", manuallyDefinedExclusions, "
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY subject_cui, subject_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND cp.subject_cui NOT IN ", manuallyDefinedExclusions, "
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ", DEMsqldef, 
"                       AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, lower(zx.subject_name) as covarname FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  #write.table(file = "covariates/rawConfounders_for_depression2AD.txt", x = confounders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(confounders)
}


#########
# getColliders for each dementia subtype
#########
#
getColliders <- function(predicates, DEMsqldef) {
  colliders <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY object_cui, object_name),
                        ZY AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ", DEMsqldef, 
"                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.object_cui, cp.object_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, lower(zx.object_name) as covarname FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  #write.table(file = "covariates/rawColliders_for_depression2AD.txt", x = colliders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(colliders)
}

#########
# Define each dementia subtype by corresponding UMLS CUIs
#########
#

#Dementia type	Definition in terms of CUIs in the UMLS
Dsql <-	" ('C0497327', 'C0011268') " # hypernym of the subtypes
ADsql <- " ('C0002395') " # Alzheimer's disease
EOADsql <- " ('C0750901', 'C0276496', 'C0751071', 'C0276496') " # Early Onset AD
LOADsql <- " ('C0494463') " # Late Onset AD
FTDsql	<- " ('C0338451') " # Frontotemporal Dementia
ADCsql <-	" ('C0001849') " # AIDS Dementia Complex 
VaDsql	<- " ('C0011269', 'C0393561') " # Vascular Dementia

#########
# Get confounders and colliders for each dementia subtype and remove any colliders 
#########
#
Dconfounders <- getConfounders(predicates, Dsql)
ADconfounders <- getConfounders(predicates, ADsql)
EOADconfounders <- getConfounders(predicates, EOADsql)
LOADconfounders <- getConfounders(predicates, LOADsql)
FTDconfounders <- getConfounders(predicates, FTDsql)
ADCconfounders <- getConfounders(predicates, ADCsql)
VaDconfounders <- getConfounders(predicates, VaDsql)

#Dcolliders <- getColliders(predicates, Dsql)
#ADcolliders <- getColliders(predicates, ADsql)
#EOADcolliders <- getColliders(predicates, EOADsql)
#LOADcolliders <- getColliders(predicates, LOADsql)
#FTDcolliders <- getColliders(predicates, FTDsql)
#ADCcolliders <- getColliders(predicates, ADCsql)
#VaDcolliders <- getColliders(predicates, VaDsql)

# confounders minus colliders
#setdiff(Dconfounders[,4], Dcolliders[,4])


#colliderConfounders <- intersect(confounders, colliders)
#colliderConfounders <- setdiff(colliderConfounders, mediators)

#########
# Set operations to compare confounders
#########
# intersection of all subtypes
intersect(Dconfounders[,4], intersect(ADconfounders[,4], intersect(EOADconfounders[,4], intersect(LOADconfounders[,4], intersect(FTDconfounders[,4], intersect(ADCconfounders[,4], VaDconfounders[,4]))))))


setdiff(Dconfounders[,4], ADconfounders[,4])
setdiff(ADconfounders[,4], Dconfounders[,4])
intersect(ADconfounders[,4], Dconfounders[,4])

intersect(VaDconfounders[,4], ADconfounders[,4])
setdiff(ADconfounders[,4], VaDconfounders[,4])
setdiff(VaDconfounders[,4], ADconfounders[,4])


Dconfounders 
ADconfounders  
EOADconfounders  
LOADconfounders  
FTDconfounders  
ADCconfounders  
VaDconfounders 



#########
# Create Venn diagrams
#########
#
library(VennDiagram)

# your data
Dconfs <- sort(Dconfounders[,4]) #c('a','b','c','d')
ADconfs <- sort(ADconfounders[,4]) #c('a','e','f','g')
EOADconfs <- sort(EOADconfounders[,4])
FTDconfs <- sort(FTDconfounders[,4])
#ADC <- ADCpvars[,3]
VaDconfs <- sort(VaDconfounders[,4])

# Generate plot
v <- venn.diagram(main = "Distribution of confounders for depression and dementia subtypes", main.cex = 2.5, #main.fontface =  "italics",
                  #sub = "for depression and dementia subtypes", sub.cex = 2.5,
                  list(AD=ADconfs, VaD=VaDconfs, dementia=Dconfs, EOAD=EOADconfs, FTD=FTDconfs),
                  fill = c("red", "yellow", "orange", "blue", "cyan"),
                  alpha = c(0.5, 0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5, cat.fontface = "bold", resolution = 300,
                  filename=NULL)

# have a look at the default plot
grid.newpage()
grid.draw(v)
