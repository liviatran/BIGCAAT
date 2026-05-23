##### VARIANT AMINO ACID EXTRACTION #####
# v 0.9.1
#creates a dataframe with variant amino acid positions for a given locus

variantAAextractor<-function(loci, dataset, exon_analyze){
  
  AA_atlas<-HLAtools::HLAatlas$prot[c(loci)]
  
  #checks if exons entered for exon specific analysis exist
  if(!is.null(exon_analyze)){
    num_exons <- seq(1:ncol(AA_atlas[[loci]])+1)
    if(any(!exon_analyze %in% num_exons)){
      return(paste("The following exons do not exist for HLA-", loci, ": ", paste(exon_analyze[!exon_analyze %in% num_exons], collapse=","), sep =""))
    }
  }
  
  for (a in 3:ncol(dataset)){
    dataset[[a]]<-ifelse(is.na(dataset[[a]])==FALSE, paste(colnames(dataset[a]),dataset[,a],sep="*"), NA)
  }
  
  #removes rows with only ALL NA data
  dataset<-dataset[!(rowSums(is.na(dataset))==ncol(dataset)-2),]
  
  #initalize variables
  variantAApositions<-genoExonList<-missingGenoOutput<-missingGeno<-repVariantAA<-mastertablecols<-mastertable<-positionParsed<-nonCWDChecked<-nonCWDtrunc<-singleAA_exon<-singleAAAlleles<-pastedAAseq<-columns<-genoCWD<-genotypeVariants<-genotypeAlleles<-locusAAsegments<-AAAligned <-refexon<-pepsplit<-alignment<-exonList<- sapply(loci, function(x) NULL)
  
  for(i in 1:length(loci)){
    
    locusAASegments<-HLAtools::buildAlignments(loci[[i]], 'AA')[[1]]$AA
    
    columnNames<- colnames(locusAASegments)
    locusAtlas<-AA_atlas[[loci[[i]]]]
    atlasCols<-ncol(AA_atlas[[loci[[i]]]])
    
    # the reference sequence for DQB1 is an insertion for the exon 4-5 boundary
    # these positions are 226.1-227 for the reference, but are actual amino
    # acid positions for DQB1*05:03 alleles
    # use 226.1 as the exon 4-5 boundary and 227 for the exon 5-6 boundary
    if(loci[[i]] == "DQB1"){
      locusAtlas[['E.4-5']] <- '226.1'
      locusAtlas[['E.5-6']] <- '227'
    }
    
    #for loop for subsetting locusAAsegments by matching exon start and end cells from AA_atlas
    #column names of locusAAsegments, which are AA positions
    #subsets relevant amino acids, inputting them into a list
    #binds previous columns with locus, allele, trimmed allele, and allele name information
    
    #HLA-A, B, and C's first exons end at -1 (i.e exon 2 begins at position 1), so 
    #the matching end atlas coordinate must be substracted by 2, since there is 
    #no position zero in the alignment
    
    #HLA-DQB1, DRB1, and DPB1's first exon ends at a number other than -1 
    #(i.e. exon 2 begins at position #2<, the matching end atlas coordinate is 
    #only subtracted by 1, since we do not need to
    #account for there being no position zero in the alignment)
    
    ##### EXON 1 COLUMN EXTRACTION #####

    exonList[[loci[[i]]]][[1]] <- cbind(locusAASegments[, 1:4], locusAASegments[,5:(match(as.numeric(locusAtlas[,1]), columnNames)-1)])

    ##### LAST EXON COLUMN EXTRACTION #####
    
    exonList[[loci[i]]][[atlasCols+1]]<-cbind(locusAASegments[,1:4], locusAASegments[match(as.numeric(locusAtlas[[atlasCols]]), columnNames):ncol(locusAASegments)])  
    
    ##### REMAINING EXON EXTRACTION ####
    #subsets N-1 exons 
    for(j in 1:(atlasCols-1)){
      exonList[[loci[i]]][[j+1]]<-cbind(locusAASegments[,1:4], locusAASegments[,match(locusAtlas[[j]], colnames(locusAASegments)):(match(as.numeric(locusAtlas[[j+1]]),columnNames)-1)])
    }
    
    #subset exonList to input exons if exon specific analysis is desired
    if(!is.null(exon_analyze)){
      exonList[[loci[[i]]]]<-exonList[[loci[[i]]]][c(exon_analyze)]
      
      #check to see if entered exon # exists for the specified locus
      exon_check <- lapply(exonList[[loci[[i]]]], is.null) == TRUE
      
      if(any(exon_check) == TRUE){
        return(paste("Exon ", exon_analyze[exon_check], " does not exist for HLA-", loci[[i]], ".", sep = ""))
      }
    }
  
    #for loop for subsetting exonList alleles to only those found in genotype data
    #focuses on subsetting via the third column in exonList, which consists of trimmed_allele data 
    #variable e in for loop represents number of columns per locus, which is how BIGDAWG input data is formatted
    for(d in 1:length(exonList[[loci[i]]])){
      for(e in 1:2){
        
        #finds which exonList alleles are present in genotype data alleles 
        genotypeAlleles[[loci[i]]][[e]]<-exonList[[loci[i]]][[d]][,3][which(exonList[[loci[i]]][[d]][,3] %in% dataset[which(colnames(dataset)%in%loci[[i]]==TRUE)][,e])]
        
      }
    }
    
    #merges both sets of unique alleles found in exonList and gets rid of duplicates
    genotypeAlleles[[loci[i]]]<-unique(append(genotypeAlleles[[loci[i]]][[1]], genotypeAlleles[[loci[i]]][[2]]))
    
    #creates a variable genoExonList, with the number of elements equal to how many exons there are for an allele
    genoExonList[[loci[i]]]<-sapply(exonList[[loci[i]]], function(x) NULL)
    
    #reads in text file of of latest, full allele history -- chooses most recent allele release to set as HLAAlleles
    #LT
    HLAAlleles<-read.csv("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist_history.txt", header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")[,c(1,2)]
    
    #compiles a list of CWD alleles and inserts them into a new variable
    CWDalleles<-CWDverify()
    
    #makes a list of lists based on the number of exons for a given locus
    nonCWDChecked[[loci[[i]]]]<-singleAA_exon[[loci[[i]]]]<-singleAAAlleles[[loci[[i]]]]<-pastedAAseq[[loci[[i]]]]<-columns[[loci[[i]]]]<-nonCWDtrunc[[loci[[i]]]]<-genotypeVariants[[loci[[i]]]]<-genoCWD[[loci[[i]]]]<-sapply(exonList[[loci[[i]]]], function(x) NULL)
    
    #subsets exonList alleles to those found in genotype data and inserts them into a new list
    #genoExonList
    for(d in 1:length(exonList[[loci[i]]])){
      
      genoExonList[[loci[i]]][[d]]<-exonList[[loci[[i]]]][[d]] %>%
        filter(trimmed_allele %in%genotypeAlleles[[loci[i]]])
      
      #associate accession # from HLA allele history with matching alleles
      genoExonList[[loci[[i]]]][[d]]<-genoExonList[[loci[[i]]]][[d]] %>%
        add_column("accessions"=HLAAlleles[,1][match(genoExonList[[loci[i]]][[d]]$allele_name, HLAAlleles[,2])], .before = 1)
      
      #based on accession number, determine if allele is CWD or not
      genoExonList[[loci[[i]]]][[d]]<-genoExonList[[loci[[i]]]][[d]] %>%
        add_column('CWD'=ifelse(genoExonList[[loci[i]]][[d]]$accessions %in% CWDalleles$Accession, "CWD", "NON-CWD"), .before = 1)
      
      #subsets genoExonList to only containing CWD alleles
      #NOTE: all g_data will be a master copy of all variants of genotype data alleles
      if(any(genoExonList[[loci[i]]][[d]]$CWD=="CWD")){
        
        genoCWD[[loci[i]]][[d]]<-genoExonList[[loci[i]]][[d]] %>%
          filter(CWD == 'CWD')
      }
      
      #compares whether all truncated alleles in genoCWD are in genotypeAlleles
      #returns truncated alleles that are not CWD, but that are present in genotypeAlleles
      nonCWDbool<-cbind.data.frame(bool=genotypeAlleles[[loci[i]]]%in%genoCWD[[loci[i]]][[d]]$trimmed_allele, allele=genotypeAlleles[[loci[i]]])
      
      nonCWDtrunc[[loci[i]]]<- nonCWDbool %>%
        filter(bool == FALSE) %>%
        pull(allele)
      
      if (length(nonCWDtrunc[[loci[i]]]) != 0) {
        
        #obtains non-CWD genotype variants in the genotype dataset
        for(b in 1:length(nonCWDtrunc[[loci[i]]])){
          
          genotypeVariants[[loci[i]]][[d]][[b]]<-genoExonList[[loci[i]]][[d]] %>%
            filter(trimmed_allele == nonCWDtrunc[[loci[i]]][[b]])
          
          #if the non-CWD allele only has one variant, bind it to genoCWD
          if(nrow(genotypeVariants[[loci[i]]][[d]][[b]])==1){
            genoCWD[[loci[[i]]]][[d]]<-rbind(genoCWD[[loci[[i]]]][[d]],genotypeVariants[[loci[[i]]]][[d]][[b]])
          }
          
          #if the non-CWD allele has more than one variant, extract number of amino acid columns
          #present for a given exon
          else{
            
            columns[[loci[i]]][[d]]<-7:length(genotypeVariants[[loci[i]]][[d]][[b]])
            
            #if an exon for a non-CWD allele has more than one amino acid column, paste all the columns together to obtain
            #the amino acid sequence which is stored in pastedAAseq
            #pastedAAseq is evaluated to find which allele variant has the most complete sequence by counting the number of
            #character, omitting * (notation for unknown amino acid)
            #the allele with the most compelte sequence is bound to genoCWD
            
            if(length(columns[[loci[i]]][[d]])>1){
              
              pastedAAseq[[loci[i]]][[d]]<-apply(genotypeVariants[[loci[i]]][[d]][[b]][ , columns[[loci[i]]][[d]]] , 1 , paste , collapse = "" )
              genotypeVariants[[loci[i]]][[d]][[b]][pastedAAseq[[loci[i]]][[d]][which.max(nchar(gsub("[*^]","",pastedAAseq[[loci[i]]][[d]])))],]
              
              genoCWD[[loci[i]]][[d]]<-rbind(genoCWD[[loci[i]]][[d]], genotypeVariants[[loci[i]]][[d]][[b]][which.max(nchar(gsub("[*^]","",pastedAAseq[[loci[i]]][[d]]))),])
            }
            
            #if an exon for a non-CWD allele has one amino acid column (i.e. exon 8 for HLA-A), store it into a separate
            #variable, singleAAAlleles
            if(length(columns[[loci[i]]][[d]])==1){
              singleAA_exon[[loci[i]]][[b]]<-genotypeVariants[[loci[i]]][[d]][[b]][ncol(genotypeVariants[[loci[i]]][[d]][[b]])==7]
              singleAAAlleles[[loci[i]]]<-singleAA_exon[[loci[i]]][lapply(singleAA_exon[[loci[i]]], length)>0]}}}
        
      }
      
      #evaluates whether a variant amino acid is present and subsets it to nonCWDChecked if there is one
      #otherwise, if nonCWDchecked only contains *, use *
      if(!is.null(unlist(singleAAAlleles[[loci[i]]]))){
        for(c in 1:length(singleAAAlleles[[loci[i]]])){
          if(any(singleAAAlleles[[loci[i]]][[c]][7:length(singleAAAlleles[[loci[i]]][[c]])]!="*")==TRUE) {
            nonCWDChecked[[loci[i]]][[c]]<-subset(singleAAAlleles[[loci[i]]][[c]], singleAAAlleles[[loci[i]]][[c]][7:length(singleAAAlleles[[loci[i]]][[c]])]!="*")[1,]
          } else{
            nonCWDChecked[[loci[i]]][[c]]<-subset(singleAAAlleles[[loci[i]]][[c]], singleAAAlleles[[loci[i]]][[c]][7:length(singleAAAlleles[[loci[i]]][[c]])]=="*")[1,]}
        }
      }
      
      #binds narrowed down non-CWD alleles for one amino acid exons and inputs it back IF there is a one columned amino acid
      #if not, nothing happens
      if(length(columns[[loci[i]]][[d]])==1){
        genoCWD[[loci[i]]][[d]]<-rbind(genoCWD[[loci[i]]][[d]][ncol(genoCWD[[loci[i]]][[d]])==7], rbind(nonCWDChecked[[loci[i]]][[1]], nonCWDChecked[[loci[i]]][[2]]))
      }
    }
    
    #creates a new variable, positionParsed, with pre-defined elements based on
    #column names in locusAAsegments (i.e. position in the peptide sequence)
    positionParsed[[loci[i]]]<-sapply(colnames(locusAASegments[,5:ncol(locusAASegments)]), function(x) NULL)
    
    #for loop to extract only variant amino acids and input them into their respective element positions
    #in positionParsed
    #extracts only amino acids that are not equal to the reference, discounting . and unknown alleles (*)
    #changed NA to "." 6/6/22
    for(a in 1:length(genoCWD[[loci[i]]])){
      for(b in 1:length(7:ncol(genoCWD[[loci[i]]][[a]]))){

        positionParsed[[loci[i]]][match(colnames(genoCWD[[loci[i]]][[a]][7:ncol(genoCWD[[loci[i]]][[a]])]), names(positionParsed[[loci[i]]]))][[b]]<-unique(subset(genoCWD[[loci[i]]][[a]][c(5,b+6)], (genoCWD[[loci[i]]][[a]][b+6]!=genoCWD[[loci[i]]][[a]][,b+6][1]) & (genoCWD[[loci[i]]][[a]][b+6] != "*") & (genoCWD[[loci[i]]][[a]][b+6] != ".")))
      }
    }
    
    #turn any length 0 elements into nrow = 0
    #length 0 elements occur if the user desires exon specific analysis
    if (!is.null(exon_analyze)) {
      positionParsed[[loci[[i]]]] <- compact(positionParsed[[loci[[i]]]])
    }
    
    #removes invariant positions (i.e elements with no rows )
    positionParsed[[loci[i]]]<-positionParsed[[loci[i]]][sapply(positionParsed[[loci[[i]]]], nrow)>0]
    
    
    for(d in 1:length(genoCWD[[loci[[i]]]])){
      for(z in 1:length(positionParsed[[loci[[i]]]])){
        if(names(positionParsed[[loci[[i]]]])[[z]] %in% colnames(genoCWD[[loci[[i]]]][[d]])==TRUE){
          positionParsed[[loci[[i]]]][[z]]<-genoCWD[[loci[[i]]]][[d]] %>% 
            select(trimmed_allele, names(positionParsed[[loci[[i]]]])[[z]])
        }
      }
    }
    
    variantAApositions[[loci[[i]]]]<-sapply(positionParsed[[loci[[i]]]], function(x) NULL)
    
    for (j in 1:length(genoCWD[[loci[[i]]]])) {
      for (k in 1:length(names(variantAApositions[[loci[[i]]]]))) {
        if (any(colnames(genoCWD[[loci[[i]]]][[j]]) == names(variantAApositions[[loci[[i]]]])[[k]])) {
          variantAApositions[[loci[[i]]]][names(variantAApositions[[loci[[i]]]]) ==
                                            names(variantAApositions[[loci[[i]]]])][[k]] <-
            cbind.data.frame(trimmed_allele = genoCWD[[loci[[i]]]][[1]][, 5],
                             genoCWD[[loci[[i]]]][[j]][colnames(genoCWD[[loci[[i]]]][[j]]) == names(variantAApositions[[loci[[i]]]])[[k]]],
                             stringsAsFactors = FALSE)
        }
      }
    }
    
    #creates a dataframe that will go into BIGDAWG,     
    #where each variant position has 2 columns to match each locus specific
    #column in genotype data
    #columns 1 and 2 of this dataframe are adapted from genotype data columns
    #patientID and disease status
    mastertable[[loci[[i]]]]<- data.frame(dataset[,c(1,2)], matrix("", ncol = length(variantAApositions[[loci[[i]]]])*2), stringsAsFactors = F)
    mastertablecols[[loci[[i]]]]<-names(positionParsed[[loci[[i]]]])
    
    #repeats variant amino acid positions twice and stores them for future naming of
    #master table column
    for(t in 1:length(mastertablecols[[loci[[i]]]])){
      repVariantAA[[loci[[i]]]][[t]]<-rep(mastertablecols[[loci[[i]]]][[t]],2)
    }
    
    #renames column names
    colnames(mastertable[[loci[[i]]]])<-c("SampleID", "Disease", unlist(repVariantAA[[loci[[i]]]]))
    
    #add variant aa positions toa mastertable
    for(u in 1:length(dataset[loci[[i]]==colnames(dataset)])){
      for(s in 1:length(variantAApositions[[loci[[i]]]])){
        mastertable[[loci[[i]]]][names(variantAApositions[[loci[[i]]]][[s]][2]) == names(mastertable[[loci[[i]]]])][[u]]<-variantAApositions[[loci[[i]]]][[s]][,2][match(dataset[loci[[i]]==colnames(dataset)][[u]], variantAApositions[[loci[[i]]]][[s]][,1])]
      }
    }
  
  ########check if all positions in mastertable contain polymorphic amino acids
  mt_length<-sapply(mastertablecols[[loci[[i]]]], function(x) NULL )
  
  #find length of unique amino acids for each pair of columns in mastertable - input into vector 
  for(y in 1:length(mastertable[[loci[[i]]]][c(TRUE, FALSE)][2:length(mastertable[[loci[[i]]]][c(TRUE, FALSE)])])){
    mt_length[[y]]<-c(length(unique(mastertable[[loci[[i]]]][c(FALSE, TRUE)][2:length(mastertable[[loci[[i]]]][c(FALSE, TRUE)])][[y]])),length(unique(mastertable[[loci[[i]]]][c(TRUE, FALSE)][2:length(mastertable[[loci[[i]]]][c(FALSE, TRUE)])][[y]])))}
  
  mt_remove<-NULL
  
  #finds amino acid positions that only have one amino acid variant
  for(z in 1:length(mt_length)){
    if(mt_length[[z]][1]==1 & mt_length[[z]][2]==1){
      mt_remove[[z]]<-names(mt_length)[[z]]
    }}
  
  
  mt_remove<-mt_remove[!is.na(mt_remove)]
  
  #finds those positions in mastertable and removes them 
  for(w in 1:length(mt_remove)){
    mastertable[[loci[[i]]]][colnames(mastertable[[loci[[i]]]]) %in% mt_remove[[w]]]<-NULL
  }
  
  cat(sprintf("Variant AA table extraction complete for HLA-%s! \n", loci[[i]]))
  
  }

  mastertable
}
