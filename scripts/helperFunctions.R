# helper functions
# v 0.9
#get CWD data
CWDverify <- function(){
  require(data.table)
  
  ## Pull down the CWD catalogue
  CWD <- list()
  CWD$data <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",skip = 1,stringsAsFactors = FALSE,select = c(2,3),showProgress = FALSE)
  CWD$version <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",nrows = 1,stringsAsFactors = FALSE,select=1,showProgress = FALSE)
  
  ## Pull down the hla_nom.txt, Deleted_alleles.txt and allelelist.txt files to create a table of v3.0.0+ deleted alleles, their ACCs,their replacements, and their ACCs
  deletedHLA <- list()
  # Temporarily store the entire hla_nom.txt in $version
  deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",skip=6, stringsAsFactors = FALSE,sep = ";", col.names = c("Locus","AlleleName","NewName","Event"),select = c(1,2,5,6),showProgress = FALSE)
  ## Exclude entries without allele name changes
  deletedHLA$data <- deletedHLA$version[deletedHLA$version$NewName !="",]
  # Exclude pre-db release 3.0.0 alleles
  deletedHLA$data <- deletedHLA$data[grep(":",deletedHLA$data$AlleleName,fixed=TRUE),]
  
  ## Process and extract the accession numbers from the Deleted_alleles.txt file, temporarily stored in $version
  deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,skip = 7,sep=",",header=TRUE,fill=TRUE,showProgress = FALSE)
  ## Below to account for one extra comma in line 106 (hopefully, can be deleted in a future release)
  if(ncol(deletedHLA$version)==4) {deletedHLA$version$Description[98] <- paste(deletedHLA$version$Description[98],deletedHLA$version$V4[98],sep=" ")
  deletedHLA$version <- deletedHLA$version[,1:3] }
  # Store the pertinent accession numbers in the data element
  deletedHLA$data$origAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$AlleleName,sep=""),deletedHLA$version$Allele)]
  # Temporarily store the allelelist.txt file in $version
  deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",skip=6, stringsAsFactors = FALSE,sep = ",", header=TRUE,showProgress = FALSE)
  deletedHLA$data$newAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep=""),deletedHLA$version$Allele)]
  # overwrite the Deleted_alelles.txt files with the version information
  deletedHLA$version <- cbind(fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE,showProgress = FALSE),fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE,showProgress = FALSE), fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",nrows=5, stringsAsFactors = FALSE,sep = "?", header=TRUE,showProgress = FALSE))
  
  ## Match accession numbers in CWD to the Accession numbers in the deleted alleles.
  changeCWD <- match(CWD$data$`IMGT/HLA Accession Number`,deletedHLA$data$origAccession)
  # Create full allele names for the new names
  deletedHLA$data$NewName <- paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep="")
  CWD$data[!is.na(changeCWD),] <- cbind(deletedHLA$data[changeCWD[!is.na(changeCWD)],6],deletedHLA$data[changeCWD[!is.na(changeCWD)],3])
  
  # Rename the columns of the verified CWD table
  colnames(CWD$data) <- c("Accession","AlleleName")
  
  CWD$data
}

dupdiff <- function(x,y) x[-match(
  make.unique(as.character(y)),
  make.unique(as.character(x)),
  nomatch=0
)]

countUniqueMotifs <- function(BCdata) {
  
  motifSize <- length(unlist(strsplit(colnames(BCdata)[3],":",fixed=TRUE)))
  
  nPairs <- (ncol(BCdata)-2)/2
  unique <- 0
  
  for(i in 1:(nPairs)) {
    unique <- unique + length(unique(unique(BCdata[,(i*2)+1]),unique(BCdata[,(i*2)+2])))
  }
  
  c(motifSize,nPairs,unique)
  
}