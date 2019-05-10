  #variantAA_extractor function
  #V 0.5
  #By: L Tran
  #05/09/19
  

  ##Source Files
  
  #Allelelist.3310.txt
  #obtained from https://github.com/ANHIG/IMGTHLA/blob/Latest/Allelelist.3310.txt
  #Various versions of Allelelist.___.txt are available depending on the version of HLA alleles used
  #This file is a list of documented HLA alleles with their allele IDs
  #needed for the alleleListfile parameter in cwdID
  
  #MS_Eur.txt
  #Patient HLA haplotypes in BIGDAWG format
  
  #AA_atlas.rda
  #an .rda file saved from AA_atlas, where AA_atlas is a set of exon boundaries for a given locus
  
  ########
  
  ##REQUIRED PACKAGES 
  require(data.table)
  require(stringr)
  require(BIGDAWG)
  
  
###REQUIRED FUNCTIONS:
#function to count spaces in between regions of interest
#to determine where  start for the alignment sequence 
  countSpaces <- function(x){
    counter <- 0
    coll <- numeric()
    vec <- strsplit(x," ")[[1]]
    for(i in 1:length(vec)){
      if (vec[i]==""){
        counter <- counter+1
      }
      else{
        if (counter!=0) coll <- c(coll,counter)
        counter <- 1
      }
    }
    coll
  }
  
#CWDverify, by SJ MACK -- compares CWD catalogue to changed and deleted files
  #outputs accurate CWD files 
  CWDverify <- function(){
    require(data.table)
    
    ## Pull down the CWD catalogue
    CWD <- list()
    CWD$data <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",skip = 1,stringsAsFactors = FALSE,select = c(2,3))
    CWD$version <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",nrows = 1,stringsAsFactors = FALSE,select=1)
    
    ## Pull down the hla_nom.txt, Deleted_alleles.txt and allelelist.txt files to create a table of v3.0.0+ deleted alleles, their ACCs,their replacements, and their ACCs
    deletedHLA <- list()
    # Temporarily store the entire hla_nom.txt in $version
    deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",skip=6, stringsAsFactors = FALSE,sep = ";", col.names = c("Locus","AlleleName","NewName","Event"),select = c(1,2,5,6))
    ## Exclude entries without allele name changes
    deletedHLA$data <- deletedHLA$version[deletedHLA$version$NewName !="",]
    # Exclude pre-db release 3.0.0 alleles
    deletedHLA$data <- deletedHLA$data[grep(":",deletedHLA$data$AlleleName,fixed=TRUE),]
    
    ## Process and extract the accession numbers from the Deleted_alleles.txt file, temporarily stored in $version
    deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,skip = 7,sep=",",header=TRUE,fill=TRUE)
    ## Below to account for one extra comma in line 106 (hopefully, can be deleted in a future release)
    if(ncol(deletedHLA$version)==4) {deletedHLA$version$Description[98] <- paste(deletedHLA$version$Description[98],deletedHLA$version$V4[98],sep=" ")
    deletedHLA$version <- deletedHLA$version[,1:3] }
    # Store the pertinent accession numbers in the data element
    deletedHLA$data$origAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$AlleleName,sep=""),deletedHLA$version$Allele)]
    # Temporarily store the allelelist.txt file in $version 
    deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",skip=6, stringsAsFactors = FALSE,sep = ",", header=TRUE)
    deletedHLA$data$newAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep=""),deletedHLA$version$Allele)]
    # overwrite the Deleted_alelles.txt files with the version information
    deletedHLA$version <- cbind(fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE),fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE), fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",nrows=5, stringsAsFactors = FALSE,sep = "?", header=TRUE))
    
    ## Match accession numbers in CWD to the Accession numbers in the deleted alleles. 
    changeCWD <- match(CWD$data$`IMGT/HLA Accession Number`,deletedHLA$data$origAccession)
    # Create full allele names for the new names
    deletedHLA$data$NewName <- paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep="")
    CWD$data[!is.na(changeCWD),] <- cbind(deletedHLA$data[changeCWD[!is.na(changeCWD)],6],deletedHLA$data[changeCWD[!is.na(changeCWD)],3])
    
    # Rename the columns of the verified CWD table
    colnames(CWD$data) <- c("Accession","AlleleName")
    
    CWD$data
  }
  
#loads AA_atlas
load("AA_atlas.rda")

###BEGIN SCRIPT FOR variantAAextractor()
  
  variantAAextractor<-function(loci,genotypefiles){
    
  #reads in genotype data  
  gdata<-read.table(genotypefiles, sep="\t", header=T, check.names = F, stringsAsFactors = F)
  
  #sets blank cells to NA 
  #if cells do not contain NA, locus names are pasted to the allele in the MS_file
  for (i in 3:ncol(gdata)){
    gdata[gdata==""]<-NA
    gdata[[i]]<-ifelse(is.na(gdata[[i]])==FALSE, paste(colnames(gdata[i]),gdata[,i],sep="*"), NA)}
  
  #removes rows with only ALL NA data 
  gdata<-gdata[!(rowSums(is.na(gdata))==NCOL(gdata)-2),]
  
  #creates empty variables for future for loops
  start<-end<-alignment<-list()
  
  #creates empty variables where each element is named after the used loci 
  
  #empty variables for correspondence table 
  inDels<-corr_table<-cols<-downloaded_segments<-w<-alignment_positions<-alignment_length<-alignment_start<-prot_extractions<-refblock_number<-end_char<-space_diff<-
    
    #empty variables for exon_extractor function   
    variantAApositions<-geno_exonlist<-missing_geno_output<-missing_geno<-rep_variantAA<-mastertablecols<-mastertable<-position_parsed<-nonCWD_checked<-nonCWDtrunc<-singleAA_exon<-singleAA_alleles<-pastedAAseq<-columns<-all_gdata<-genotype_variants<-geno_alleles<-AA_segments<-AA_aligned <-refexon<-pepsplit<-alignment<-exonlist<- sapply(loci, function(x) NULL)
  
  #begin for loop   
  for(i in 1:length(loci)){
    #downloads relevant locus alignment file -- readLines allows for space preservation, which is important in
    #finding where the alignment sequence starts 
    alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(ifelse(loci[[i]]=="DRB1","DRB",loci[[i]]),"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
    
    #alters alignment file by cutting out non-pertinent information in beginning
    #and endind of alignment file 
    alignment[[loci[i]]] <- head(alignment[[loci[i]]],-3)
    alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-7)
    
    #see countSpaces function at beginning of script 
    #Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
    #the alignment sequence actually starts 
    space_diff[[loci[i]]]<-(nchar(strsplit(alignment[[loci[i]]][3], " ")[[1]][2])+countSpaces(alignment[[loci[i]]][3])[2]+1)-countSpaces(alignment[[loci[i]]][2])[1]
    
    #reduces repeated whitespace in alignment file and removes rows with empty values for proper
    #start and stop subsetting 
    alignment[[loci[i]]] <-str_squish(alignment[[loci[i]]])
    alignment[[loci[i]]] <-alignment[[loci[i]]][-which(alignment[[loci[i]]] == "")]
    
    #determines positions of "Prot" and the end of that reference block segment
    start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
    end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,length(alignment[[loci[i]]])))
    
    #counts number of characters in the very last allele to add onto the last Prot enumeration block
    #to obtain end length 
    end_char[[loci[i]]]<-nchar(sapply(strsplit(gsub(" ", "", sub(" ", "~", str_squish(tail(alignment[[loci[i]]], 1)))), "~"), "[", 2))-1
    
    #extracts rows with "Prot" and reference sequence position information 
    #extracts only relevant reference sequence positions
    #NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
    #as the actual sequence start will always be 1 
    for (j in 1:length(start[[loci[i]]])){
      
      prot_extractions[[loci[i]]][j]<-strsplit(alignment[[loci[i]]][start[[loci[i]]][j]], " ")
      
      refblock_number[[loci[i]]][j]<-as.numeric(sapply(prot_extractions[[loci[i]]][j], "[", 2))
      
      
      #determines the alignment start by adding -30 to the difference between white spaces found above 
      alignment_start[[loci[i]]]<-refblock_number[[loci[i]]][1]+space_diff[[loci[i]]]
    }
    
    #closes all white space in the alignment file, except for the white space separating the allele and peptide sequence
    alignment[[loci[i]]] <-paste(substr(alignment[[loci[i]]],1,regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(alignment[[loci[i]]],regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE),nchar(alignment[[loci[i]]]))),sep = "")
    
    #string splits at white spaces to yield allele and peptide sequences
    alignment[[loci[i]]]  <- strsplit(alignment[[loci[i]]]," ", fixed=T)
    
    #binds the previously split strings by row 
    alignment[[loci[i]]] <- do.call(rbind,alignment[[loci[i]]])
    
    #if the pepseq column is equal to the allele column due to premature peptide termination, 
    #insert a blank in place of the allele in the pepseq column 
    alignment[[loci[i]]][which(alignment[[loci[i]]][,1]==alignment[[loci[i]]][,2]),2] <- ""
    
    #renames columns to "alleles" and "pepseq"
    colnames(alignment[[loci[i]]])<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")
    
    #due to ANHIG formatting, cases where an allele contains newly reference peptide sequences will not 
    #contain the same number of rows as previous reference peptide blocks
    #this for loop is invoked to add "."for all other alleles for each character in the newly reference peptide
    #to preserve structural integrity 
    for(k in 1:length(start[[loci[i]]])){
      if(nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
        x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
        colnames(x)<-paste(loci[[i]], "alleles", sep="_")
        x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(alignment[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
        y<-data.frame(tail(alignment[[loci[i]]],1), stringsAsFactors = F)
        x$pepseq[match(y[,1], x[,1])]<-y$pepseq
        alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -1), x))
        start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
        end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))}}
    
    #if a locus has extra formatting, resulting in unqeual rows, start and end will be updated to reflect subsetting
    #if a locus has no extra formatting, start and end will remain the same, as procured by earlier code
    for(e in 1:length(start[[loci[i]]])){
      AA_segments[[loci[i]]]<-cbind(AA_segments[[loci[i]]], alignment[[loci[i]]][start[[loci[i]]][e]:end[[loci[i]]][e],])}
    
    #removes first two rows containing AA position and "Prot"
    AA_segments[[loci[i]]] <- AA_segments[[loci[i]]][-c(1,2),]
    
    #designates columns to be combined as every other so allele names are not included
    #in pasting all the amino acid sequences together 
    cols<-seq(0, ncol(AA_segments[[loci[i]]]), by=2)
    AA_segments[[loci[i]]]<-cbind(AA_segments[[loci[i]]][,1], apply(AA_segments[[loci[i]]][,cols], 1 ,paste, collapse = ""))
    
    #creates a new matrix with the number of columns equal to the number of characters in the reference sequence 
    corr_table[[loci[i]]]<-matrix(, nrow = 2, ncol = as.numeric(nchar(AA_segments[[loci[i]]][,2][1])))
    
    #determines alignment length based on the total number of characters plus the alignment start (which is negative ) 
    alignment_length[[loci[i]]]<-as.numeric(nchar(AA_segments[[loci[i]]][,2][1]))+alignment_start[[loci[[i]]]]
    
    #pastes alignment_start to alignment_length together in sequential order, with inDels accounted for 
    #captures output as "w"
    w[[i]] <- capture.output(cat(alignment_start[[loci[i]]]:alignment_length[[loci[i]]]))
    
    #splits string formed by cat for separate character variables
    alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(w[[loci[i]]], " ")))
    
    #eliminates "0", as the alignment sequence from ANHIG does not contain 0
    alignment_positions[[loci[i]]]<-alignment_positions[[loci[i]]][-which(alignment_positions[[loci[i]]] == 0)]
    
    #contains alignment sequence information 
    corr_table[[loci[i]]][2,]<-alignment_positions[[loci[i]]]
    
    #string splits to extract locus in the allele name
    #assigns to new variable "AA_aligned"
    AA_aligned[[loci[i]]]<- as.matrix(do.call(rbind,strsplit(AA_segments[[loci[i]]][,1],"[*]")))
    
    #adds a new column of pasted locus and trimmed two field alleles to AA_aligned
    AA_aligned[[loci[i]]]<- cbind(AA_aligned[[loci[i]]], paste(AA_aligned[[loci[i]]][,1], apply(AA_aligned[[loci[i]]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2], sep="*"))
    
    #binds AA_aligned and AA_segments -- renames columns 
    AA_segments[[loci[i]]] <- cbind(AA_aligned[[loci[i]]], AA_segments[[loci[i]]])
    colnames(AA_segments[[loci[i]]]) <- c("locus", "full_allele", "trimmed_allele", "allele_name", "AAsequence")
    
    #sets refexon to a reference peptide for each HLA locus based on the reference sequences in AA_segments 
    refexon[[loci[i]]] <- rbind(AA_segments[[loci[i]]][1,])[which(rbind(AA_segments[[loci[i]]][1,])[,"locus"]==loci[[i]]),'AAsequence']
    
    #splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
    pepsplit[[loci[i]]] <- sapply(AA_segments[[loci[i]]][,"AAsequence"],strsplit,split="*")
    
    #fills in space with NA for alleles with premature termination to make it the same number of characters
    #as the reference sequence 
    pepsplit[[loci[i]]]<- lapply(pepsplit[[loci[i]]],function(x) c(x,rep("NA",nchar(refexon[[loci[i]]])-length(x))))
    
    #binds pep_split together by element in its previous list form by row
    pepsplit[[loci[i]]]<- do.call(rbind,pepsplit[[loci[i]]])
    
    #nullifies row names 
    rownames(pepsplit[[loci[i]]]) <- NULL
    
    #binds all columns together to form desired ouput, as described above
    AA_segments[[loci[i]]] <- cbind.data.frame(AA_segments[[loci[i]]][,1:4],pepsplit[[loci[i]]], stringsAsFactors=FALSE)
    
    #finds positions in AA_segments that have ".", indicating an inDel 
    inDels[[loci[[i]]]]<-colnames(AA_segments[[loci[[i]]]][1, 5:ncol(AA_segments[[loci[[i]]]])][AA_segments[[loci[[i]]]][1, 5:ncol(AA_segments[[loci[[i]]]])] %in% "."])
    
    #inputs AA_segments alignment sequence into the corr_table with "InDel" still present
    corr_table[[loci[[i]]]][1,]<-names(AA_segments[[loci[[i]]]][5:ncol(AA_segments[[loci[[i]]]])])
    
    for(b in 1:length(inDels[[loci[[i]]]])){
      corr_table[[loci[[i]]]][2,][inDels[[loci[[i]]]][[b]]==corr_table[[loci[[i]]]][1,]]<-paste("InDel", b, sep="_")
    }
    
    #fixes enumerations following "InDel"
    corr_table[[loci[[i]]]][2,][!grepl("InDel", corr_table[[loci[[i]]]][2,])]<-(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("InDel", corr_table[[loci[[i]]]][2,])]))+alignment_start[[loci[[i]]]]))[!(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("InDel", corr_table[[loci[[i]]]][2,])]))+alignment_start[[loci[[i]]]]))==0]
    
    #renames columns in AA_segments
    colnames(AA_segments[[loci[i]]]) <- c("locus","allele","trimmed_allele","allele_name", 1:ncol(corr_table[[loci[[i]]]]))
    
    #distributes  reference sequence from row 1
    #into all other rows, if they contain a "-"
    #amino acids with changes will not be impacted
    for(k in 5:ncol(AA_segments[[loci[i]]])) {
      AA_segments[[loci[i]]][,k][which(AA_segments[[loci[i]]][,k]=="-")] <- AA_segments[[loci[i]]][,k][1]}  
    
  
    #for loop for subsetting AA_segments by matching exon start and end cells from AA_atlas
    #column names of AA_segments, which are AA positions
    #subsets relevant amino acids, inputting them into a list
    #binds previous columns with locus, allele, trimmed allele, and allele name information
    
    #subsets first exon for all loci
    #HLA-A, B, and C's first exons end at -1 (i.e exon 2 begins at position 1), so 
    #the matching end atlas coordinate must be substracted by 2, since there is 
    #no position zero in the alignment
    
    #HLA-DQB1, DRB1, and DPB1's first exon ends at a number other than -1 
    #(i.e. exon 2 begins at position #2<, the matching end atlas coordinate is 
    #only subtracted by 1, since we do not need to
    #account for there being no position zero in the alignment)
    if((loci[[i]]=="A") || (loci[[i]]=="B") || (loci[[i]]=="C")){
      exonlist[[i]][[1]]<-cbind(AA_segments[[loci[i]]][,1:4], AA_segments[[loci[i]]][,5:match(as.numeric(AA_atlas[match(loci[[i]],names(AA_atlas))][[loci[i]]][[2]][[1]]-1), colnames(AA_segments[[loci[i]]]))])}
  
      if((loci[[i]]=="DRB1") || (loci[[i]]=="DQB1") || (loci[[i]]=="DPB1")){
      exonlist[[i]][[1]]<-cbind(AA_segments[[loci[i]]][,1:4], AA_segments[[loci[i]]][,5:match(as.numeric(AA_atlas[match(loci[[i]],names(AA_atlas))][[loci[i]]][[2]][[1]]), colnames(AA_segments[[loci[i]]]))])}
  
    #subsets last exon for loci 
    exonlist[[loci[i]]][[nrow(AA_atlas[[match(loci[[i]],names(AA_atlas))]])+1]]<-cbind(AA_segments[[loci[i]]][,1:4], AA_segments[[loci[i]]][match(AA_atlas[[match(loci[[i]],names(AA_atlas))]][[2]][[length(AA_atlas[match(loci[[i]],names(AA_atlas))][[loci[i]]][[2]])]]:names(AA_segments[[loci[i]]][ncol(AA_segments[[loci[i]]])]), colnames(AA_segments[[loci[i]]]))])
    
    #subsets N-1 exons 
    for(j in 1:(nrow(AA_atlas[[match(loci[i],names(AA_atlas))]])-1)){
      exonlist[[loci[i]]][[j+1]]<-cbind(AA_segments[[loci[i]]][,1:4], AA_segments[[loci[i]]][,match(AA_atlas[match(loci[i],names(AA_atlas))][[loci[i]]][[2]][[j]], colnames(AA_segments[[loci[i]]])):match(as.numeric(AA_atlas[match(loci[i],names(AA_atlas))][[loci[i]]][[2]][[j+1]]),colnames(AA_segments[[loci[i]]]))])}
    
    #for loop for subsetting exonlist alleles to only those found in genotype data
    #focuses on subsetting via the third column in exonlist, which consists of trimmed_allele data 
    #variable e in for loop represents number of columns per locus, which is how BIGDAWG input data is formatted
    for(d in 1:length(exonlist[[loci[i]]])){
      for(e in 1:2){
        
        #finds which exonlist alleles are present in genotype data alleles 
        geno_alleles[[loci[i]]][[e]]<-exonlist[[loci[i]]][[d]][,3][which(exonlist[[loci[i]]][[d]][,3] %in% gdata[which(colnames(gdata)%in%loci[[i]]==TRUE)][,e]==TRUE)]
      }}
    
    #merges both sets of unique alleles found in exonlist and gets rid of duplicates 
    geno_alleles[[loci[i]]]<-unique(append(geno_alleles[[loci[i]]][[1]], geno_alleles[[loci[i]]][[2]]))
    
    #creates a variable geno_exonlist, with the number of elements equal to how many exons there are for an allele
    geno_exonlist[[loci[i]]]<-sapply(exonlist[[loci[i]]], function(x) NULL)
    
    #reads in 3310 HLA alleles 
    HLA_alleles<-read.csv("Allelelist.3340.txt", header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")
    
    #compiles a list of CWD alleles and inserts them into a new variable
    CWDalleles<-CWDverify()
    
    #makes a list of lists based on the number of exons for a given locus 
    nonCWD_checked[[loci[[i]]]]<-singleAA_exon[[loci[[i]]]]<-singleAA_alleles[[loci[[i]]]]<-pastedAAseq[[loci[[i]]]]<-columns[[loci[[i]]]]<-all_gdata[[loci[[i]]]]<-nonCWDtrunc[[loci[[i]]]]<-genotype_variants[[loci[[i]]]]<-sapply(exonlist[[loci[[i]]]], function(x) NULL)
    
    #subsets exonlist alleles to those found in genotype data and inserts them into a new list
    #geno_exonlist
    for(d in 1:length(exonlist[[loci[i]]])){
      geno_exonlist[[loci[i]]][[d]]<-subset(exonlist[[loci[i]]][[d]], exonlist[[loci[i]]][[d]][,3]%in%geno_alleles[[loci[i]]])
      geno_exonlist[[loci[i]]][[d]]<-cbind.data.frame("accessions"=HLA_alleles$AlleleID[match(geno_exonlist[[loci[i]]][[d]]$allele_name, HLA_alleles$Allele)], geno_exonlist[[loci[i]]][[d]], stringsAsFactors=FALSE)
      geno_exonlist[[loci[i]]][[d]]<-cbind.data.frame("CWD"=ifelse(geno_exonlist[[loci[i]]][[d]]$accessions %in% CWDalleles$Accession, "CWD", "NON-CWD"), geno_exonlist[[loci[i]]][[d]], stringsAsFactors=FALSE)
      
      
      #subsets geno_exonlist to only containing CWD alleles via accession number
      #and stores it to a new variable, all_gdata
      #NOTE: all g_data will be a master copy of all variants of genotype data alleles
      if(any(geno_exonlist[[loci[i]]][[d]]$CWD=="CWD")){
        all_gdata[[loci[i]]][[d]]<-na.omit(geno_exonlist[[loci[i]]][[d]][geno_exonlist[[loci[i]]][[d]]$accessions%in%CWDalleles$Accession,])}
      
      #compares whether all truncated alleles in all_gdata are in geno_alleles
      #returns truncated alleles that are not CWD, but that are present in geno_alleles
      nonCWDtrunc[[loci[i]]]<-cbind(geno_alleles[[loci[i]]]%in%all_gdata[[loci[i]]][[d]]$trimmed_allele, geno_alleles[[loci[i]]])[which(cbind(geno_alleles[[loci[i]]], geno_alleles[[loci[i]]]%in%all_gdata[[loci[i]]][[d]]$trimmed_allele)==FALSE)]
      
      #obtains non-CWD genotype variants in the genotype dataset
      for(b in 1:length(nonCWDtrunc[[loci[i]]])){
        genotype_variants[[loci[i]]][[d]][[b]]<-subset(geno_exonlist[[loci[i]]][[d]], geno_exonlist[[loci[i]]][[d]]$trimmed_allele==nonCWDtrunc[[loci[i]]][[b]])
        
        #if the non-CWD allele only has one variant, bind it to all_gdata
        if(nrow(genotype_variants[[loci[i]]][[d]][[b]])==1){all_gdata[[loci[[i]]]][[d]]<-rbind(all_gdata[[loci[[i]]]][[d]],genotype_variants[[loci[[i]]]][[d]][[b]])}
        
        #if the non-CWD allele has more than one variant, extract number of amino acid columns
        #present for a given exon 
        if(nrow(genotype_variants[[loci[i]]][[d]][[b]])>1){
          columns[[loci[i]]][[d]]<-7:length(genotype_variants[[loci[i]]][[d]][[b]])
          
          #if an exon for a non-CWD allele has more than one amino acid column, paste all the columns together to obtain
          #the amino acid sequence which is stored in pastedAAseq
          #pastedAAseq is evaluated to find which allele variant has the most complete sequence by counting the number of
          #character, omitting * (notation for unknown amino acid)
          #the allele with the most compelte sequence is bound to all_gdata
          if(length(columns[[loci[i]]][[d]])>1){
            pastedAAseq[[loci[i]]][[d]]<-apply(genotype_variants[[loci[i]]][[d]][[b]][ , columns[[loci[i]]][[d]]] , 1 , paste , collapse = "" )
            all_gdata[[loci[i]]][[d]]<-rbind(all_gdata[[loci[i]]][[d]], genotype_variants[[loci[i]]][[d]][[b]][names(pastedAAseq[[loci[i]]][[d]][which.max(nchar(gsub("[*^]","",pastedAAseq[[loci[i]]][[d]])))]),])}
          
          
          #if an exon for a non-CWD allele has one amino acid column (i.e. exon 8 for HLA-A), store it into a separate
          #variable, singleAA_alleles
          if(length(columns[[loci[i]]][[d]])==1){
            singleAA_exon[[loci[i]]][[b]]<-genotype_variants[[loci[i]]][[d]][[b]][ncol(genotype_variants[[loci[i]]][[d]][[b]])==7]
            singleAA_alleles[[loci[i]]]<-singleAA_exon[[loci[i]]][lapply(singleAA_exon[[loci[i]]], length)>0]}}}
      
      
      #evaluates whether a variant amino acid is present and subsets it to nonCWD_checked if there is one
      #otherwise, if nonCWDchecked only contains *, use *
      for(c in 1:length(singleAA_alleles[[loci[i]]])){
        if(any(singleAA_alleles[[loci[i]]][[c]][7:length(singleAA_alleles[[loci[i]]][[c]])]!="*")==TRUE) {nonCWD_checked[[loci[i]]][[c]]<-subset(singleAA_alleles[[loci[i]]][[c]], singleAA_alleles[[loci[i]]][[c]][7:length(singleAA_alleles[[loci[i]]][[c]])]!="*")[1,]}
        if(any(singleAA_alleles[[loci[i]]][[c]][7:length(singleAA_alleles[[loci[i]]][[c]])]!="*")==FALSE){nonCWD_checked[[loci[i]]][[c]]<-subset(singleAA_alleles[[loci[i]]][[c]], singleAA_alleles[[loci[i]]][[c]][7:length(singleAA_alleles[[loci[i]]][[c]])]=="*")[1,]}
      }
      
      #binds narrowed down non-CWD alleles for one amino acid exons and inputs it back IF there is a one columned amino acid
      #if not, nothing happens 
      if(length(columns[[loci[i]]][[d]])==1)
      {all_gdata[[loci[i]]][[d]]<-rbind(all_gdata[[loci[i]]][[d]][ncol(all_gdata[[loci[i]]][[d]])==7], rbind(nonCWD_checked[[loci[i]]][[1]], nonCWD_checked[[loci[i]]][[2]]))}}
    
    
    #creates a new variable, position_parsed, with pre-defined elements based on
    #column names in AA_segments (i.e. position in the peptide sequence)
    position_parsed[[loci[i]]]<-sapply(colnames(AA_segments[[loci[i]]][,5:ncol(AA_segments[[loci[i]]])]), function(x) NULL)
    
    #for loop to extract only variant amino acids and input them into their respective element positions
    #in position_parsed 
    #extracts only variant amino acids, discounting NA and unknown alleles (*)
    for(a in 1:length(all_gdata[[loci[i]]])){
      for(b in 1:length(7:ncol(all_gdata[[loci[i]]][[a]]))){
        position_parsed[[loci[i]]][match(colnames(all_gdata[[loci[i]]][[a]][7:ncol(all_gdata[[loci[i]]][[a]])]), names(position_parsed[[loci[i]]]))][[b]]<-unique(subset(all_gdata[[loci[i]]][[a]][c(5,b+6)], (all_gdata[[loci[i]]][[a]][b+6]!=all_gdata[[loci[i]]][[a]][,b+6][1]) & (all_gdata[[loci[i]]][[a]][b+6] != "*") & (all_gdata[[loci[i]]][[a]][b+6] != "NA")))}}
    
    #removes invariant positions (i.e elements with no rows )
    #inDels will be filtered out via a is.null application
    position_parsed[[loci[i]]]<-position_parsed[[loci[i]]][sapply(position_parsed[[loci[[i]]]][which(lapply(position_parsed[[loci[[i]]]], is.null)==FALSE)], nrow)>0]
    
    #further subsets position_parsed to only variant positions with polymorphic amino acids 
    for(g in 1:length(position_parsed[[loci[i]]])){
      position_parsed[[loci[i]]][[g]]<-subset(position_parsed[[loci[i]]][[g]], length(unique(position_parsed[[loci[i]]][[g]][,2]))!=1)}
    
    #removes elements without polymorphic amino acids 
    position_parsed[[loci[i]]]<-position_parsed[[loci[i]]][sapply(position_parsed[[loci[i]]], nrow)>0]
    
    
    variantAApositions[[loci[[i]]]]<-sapply(position_parsed[[loci[[i]]]], function(x) NULL)
    
    for(j in 1:length(all_gdata[[loci[[i]]]])){
      for(k in 1:length(names(variantAApositions[[loci[[i]]]]))){
        if(any(colnames(all_gdata[[loci[[i]]]][[j]])==names(variantAApositions[[loci[[i]]]])[[k]])){variantAApositions[[loci[[i]]]][names(variantAApositions[[loci[[i]]]])==names(variantAApositions[[loci[[i]]]])][[k]]<-cbind.data.frame(trimmed_allele=all_gdata[[loci[[i]]]][[1]][,5], all_gdata[[loci[[i]]]][[j]][colnames(all_gdata[[loci[[i]]]][[j]])==names(variantAApositions[[loci[[i]]]])[[k]]], stringsAsFactors=FALSE)}}}
    
    
    #creates a dataframe that will go into BIGDAWG,     #where each variant position has 2 columns to match each locus specific
    #column in genotype data
    #columns 1 and 2 of this dataframe are adapted from genotype data columns
    #patientID and disease status 
    mastertable[[loci[[i]]]]<- data.frame(gdata[,c(1,2)], matrix("", ncol =   length(variantAApositions[[loci[[i]]]])*2), stringsAsFactors = F)
    mastertablecols[[loci[[i]]]]<-names(position_parsed[[loci[[i]]]])
    
    #repeats variant amino acid positions twice and stores them for future naming of
    #master table column 
    for(t in 1:length(mastertablecols[[loci[[i]]]])){
      rep_variantAA[[loci[[i]]]][[t]]<-rep(mastertablecols[[loci[[i]]]][[t]],2)}
    
    #renames column names 
    colnames(mastertable[[loci[[i]]]])<-c("SampleID", "Disease", unlist(rep_variantAA[[loci[[i]]]]))
    
    
    
    for(u in 1:length(gdata[loci[[i]]==colnames(gdata)])){
      for(s in 1:length(variantAApositions[[loci[[i]]]])){
        mastertable[[loci[[i]]]][names(variantAApositions[[loci[[i]]]][[s]][2]) == names(mastertable[[loci[[i]]]])][[u]]<-variantAApositions[[loci[[i]]]][[s]][,2][match(gdata[loci[[i]]==colnames(gdata)][[u]], variantAApositions[[loci[[i]]]][[s]][,1])]
      }}
  }}


#example run of the function saved to variantAAtable 
variantAAtable<-variantAAextractor(loci="A", genotypefiles = "MS_EUR.txt")

