#
#' create a SummarizedExperiment dataset from metadata and expression data
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param metadata Path to the metadata file
#' @param dataFile Path to the data file
#' @param species name of the experimental species e.g. human, mouse, etc
#' @param featureIdformat format for featureIDs in your data e.g. SYMBOL, ENSEMBLID, e.t.c
#' @return A summarizedExperiment object of the infile
#' @export
#'
dataPackaging<-function(metadataFile,dataFile,species,featureIdformat){
exp_desc<-data.table::fread(metadataFile,sep = ",",stringsAsFactors = F,header = T)
metaData<-exp_desc
#
colnames(metaData)[1]<-"label"
#
xdata<-data.table::fread(dataFile,sep = ",",stringsAsFactors = F,header = T)
#
#

# Choose what variables you want to model
# --------------------------------------------------------------
ConditionVariable<-readline("Enter condition name or column number you want to visualize: ")
#
exp_des<-dplyr::select(metaData, sample=label,all_of(ConditionVariable))

colnames(exp_des)<-c("label","condition")
exp_des<- dplyr::mutate(exp_des,
                        condition = factor(condition, levels = sort(unique(condition)),
                                           labels =sort(unique(condition))))


#
#----------------------- determine replication levels for each condition---------------------------
#
replicate<-rnorm(dim(exp_des)[1])
levs<-levels(exp_des$condition)
for (lev in 1:length(levs)){
  ii<-which(exp_des$condition %in% levs[lev])
  replicate[ii]<-seq(c(1:length(ii)))
}
exp_des<-dplyr::mutate(exp_des,replicates = replicate)
#
#---------------------------------------------------------------------------------------------------
#
sampNams<-colnames(xdata)[2:dim(xdata)[2]]
corresp<-all(sampNams==exp_des$label)
if(corresp)
{ print("SUCCESS: All sample names in the metadata match those in the data file.")
}else{
  cat("ERROR: One or more sample names in the metadata do not match those in the data file. Check the following list for FALSE cases and try again...\n\n")
  corresp
}
ColsToAnalyze <- exp_des$label
colData<-metaData
colnames(colData)[1]<-"label"
colData<-dplyr::mutate(colData,exp_des)#merge(colData,exp_des)
#
newColIDs<-paste0("counts.",ColsToAnalyze)
colnames(xdata)[2:dim(xdata)[2]]<-newColIDs
colnames(xdata)[1]<-"Gene.names"
#
#
# ---(1)--------------------------GENE ANNOTATIONS-------------------------------------------------------------
#
#spec<-c("human","mouse","rat","fruitfly","zebrafish")
if(tolower(species)=="human"){spec="hsapiens"}
if(tolower(species)=="mouse"){spec="musculus"}
if(tolower(species)=="rat"){spec="norvegicus"}
if(tolower(species)=="zebrafish"){spec="rerio"}
if(tolower(species)=="fruitFly"){spec="melanogaster"}
if(tolower(species)=="yeast"){spec="cerevisiae"}
if(tolower(species)=="cow"){spec="taurus"}
if(tolower(species)=="c.elegans"){spec="elegans"}
if(tolower(species)=="frog"){spec="tropicalis"}
#
#ensembl <- biomaRt::useMart("ensembl")
ensembl<-biomaRt::useEnsembl(biomart = "ensembl", mirror = "uswest")
ensembl <-biomaRt::useDataset(biomaRt::searchDatasets(mart = ensembl, pattern = spec)$dataset[1],mart=ensembl)
#
filterValues <- as.character(xdata$Gene.names)
filterType<-{}
if(tolower(featureIdformat)=="symbol"){filterType="external_gene_name"}
if(tolower(featureIdformat)=="entrezid"){filterType="entrezgene_id"}
if(tolower(featureIdformat)=="ensemblid"){filterType="ensembl_gene_id"}
#
attributeNames <- c('ensembl_gene_id','entrezgene_id' ,'external_gene_name','chromosome_name','start_position','end_position','strand')
annot <- biomaRt::getBM(attributes=attributeNames,
               filters = filterType,
               values = filterValues,
               mart = ensembl)

colnames(annot)<-c("ENSEMBLID","ENTREZID","SYMBOL","Chr","start","end","strand")
colnames(xdata)[1]<-toupper(featureIdformat)
#Annot<-annot[which(annot$ENSEMBLID %in% xdata$Gene.names),]
# Determine the indices for the non-NA genes
xdata<-as.data.frame(xdata)
ii<-match(xdata[,toupper(featureIdformat)],annot[,toupper(featureIdformat)])
annot<-annot[ii,]
non_na_idx <- which(is.na(annot$SYMBOL) == FALSE)

# Return only the genes with annotations using indices
annot <- annot[non_na_idx, ]

xdata2<-merge(annot,xdata)
#rownames(xdat2)<-xdata2$ENSEMBLID
#
#########------------------------------------------------------------
make_unique <- function(rawdata, names, ids, delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(rawdata),
                          is.character(names),
                          length(names) == 1,
                          is.character(ids),
                          length(ids) == 1,
                          is.character(delim),
                          length(delim) == 1)

  col_names <- colnames(rawdata)
  # Show error if inputs do not contain required columns
  if(!names %in% col_names) {
    stop("'", names, "' is not a column in '",
         deparse(substitute(rawdata)), "'",
         call. = FALSE)
  }
  if(!ids %in% col_names) {
    stop("'", ids, "' is not a column in '",
         deparse(substitute(rawdata)), "'",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(rawdata))
    rawdata <- as.data.frame(rawdata)

  # Select the name and id columns, and check for NAs
  double_NAs <- apply(rawdata[,c(names, ids)], 1, function(x) all(is.na(x)))
  if(any(double_NAs)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }

  # Take the first identifier per row and make unique names.
  # If there is no name, the ID will be taken.
 # unique_features <- rawdata %>%
   unique_features<- dplyr::mutate(rawdata,name = gsub(paste0(delim, ".*"), "", get(names)),
           ID = gsub(paste0(delim, ".*"), "", get(ids)),
           name = make.unique(ifelse(name == "" | is.na(name), ID, name)))

  return(unique_features)
}
########### ------------------------------------------------------------

data_unique <- make_unique(xdata2, "ENSEMBLID", "ENTREZID",  delim = ";")

#####
#########--------------------------------------------------------------
#
count_columns <- grep("counts.", colnames(data_unique))
# Select the assay data
rownames(data_unique) <- data_unique$name
raw <- data_unique[, count_columns]
raw[raw == 0] <- NA
raw <- log2(raw)
# Generate the colData from the experimental design
# and match these with the assay data
#import::from(magrittr, "%>%")
colnames(raw)<-sampNams
expdesign<-exp_des
rownames(expdesign) <- expdesign$label
#
matched <- match(make.names((expdesign$label)),
                 make.names((colnames(raw))))
colnames(raw)[matched] <- expdesign$label
raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
# Select the rowData
row_data <- data_unique[, -count_columns]
rownames(row_data) <- row_data$name

# Generate the SummarizedExperiment object
se_data <- SummarizedExperiment(assays = as.matrix(raw),
                                colData = expdesign,
                                rowData = row_data)

return(se_data)
}

