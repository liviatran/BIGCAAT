# wrapper function for combiAnalyzer
# v 0.9
runCombiAnalyzer <- function(loci, variantAAtable, loop) {

  BOLO_list<-KDLO_list<-UMLO_list<-list()
  
  #sets motif_list to NULL
  motif_list<-NULL
  
  myData<-variantAAtable[[loci]]
  
  #initiates recursion with stop=FALSE and begins the counter at 0
  stop<-FALSE
  
  counter<-0
  
  totalStart<-proc.time()
  
  ###BEGIN RECURSION -- as long as stop==FALSE, combiAnalyzer will be run until the maximum OR
  #is reached, or the end of the motif_list is reached
  #the recursive program receives input from combiAnalyzer, where stop=TRUE once the maximum OR
  #is reached, either because the BOLO is empty, the KDLO is empty, or new combination names
  #can be made
  while(stop==FALSE){
    
    Start<-proc.time()
    
    cat("Evaluating",ifelse(counter==0,"initial comparison of 1-mers to null hypothesis \n",paste(counter+1,"-mers \n",sep="")))
    
    interim<-combiAnalyzer(loci, myData, BOLO ,KDLO, UMLO, counter, motif_list, KDLO_list, UMLO_list, variantAAtable, loop)
    
    counter<-counter+1
    
    myData<-interim$combidf
    
    KDLO<-KDLO_list[[counter]]<-interim$KDLO
    BOLO<-BOLO_list[[counter]]<-interim$BOLO
    UMLO<-UMLO_list[[counter]]<-interim$UMLO
    
    if(is.null(nrow(KDLO))==TRUE){
      cat("Maximal significant OR values identified. End of analysis of the",loci,"locus.\n\n") ### SJM cosmetic & informative changes
      Results <- (list(KDLO = KDLO_list, BOLO = BOLO_list, UMLO = UMLO_list))
      return (Results)
    }
    
    if((is.null(nrow(KDLO))==FALSE) & (length(motif_list)!=counter)){
      ##    cat("BIGCAAT: Dataset is able to be further analyzed - moving on to next iteration.\n") ### SJM added break, and removed message
    }
    
    if((is.null(nrow(KDLO))==FALSE) & length(motif_list)==counter){
      cat("BIGCAAT: WARNING: end of motif_list analysis, but further analysis is possible.\n") ### SJM added break
      stop=TRUE
    }
    
    if((is.null(nrow(KDLO))==TRUE) & length(motif_list)==counter){
      cat("BIGCAAT: End of motif_list analysis - maximal OR has been reached.\n") ### SJM added break
    }
    
    print(sprintf('elapsed time for combiAnalyzer %s seconds', round((proc.time()-Start)['elapsed'])))
  }
}

# combiAnalyzer 
combiAnalyzer<-function(loci, myData, KDLO, BOLO, UMLO, counter, motif_list, KDLO_list, UMLO_list, variantAAtable, loop){

  #specifies a default motif list if one is not provided
  if((is.null(motif_list)==TRUE)&(counter==0)){
    motif_list<-c(0,2,3,4,5,6,7)
    #  cat("BIGCAAT: A motif list has not been provided - BIGCAAT will run until maximal OR is reached. \n") ### SJM Currently no way to provide a motif list
  }
  
  ######################################################## BIGDAWG ######################################################## 
  
  #set output as T for statistical outputs
  silenceBD <- capture.output(BOLO<-BIGDAWG(myData, HLA=F, Run.Tests="L", Missing = 2, Cores.Lim = 3L, Return=T, Output = F, Verbose = F, Data.Type = 'Motif')) ### SJM Verbose OFF, and BIGDAWG output captured to silenceBD
  
  BOLO<-data.frame(lapply(as.data.frame(BOLO$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)
  
  #filter out NCalc entries before calculating OR differences
  BOLO<-subset(BOLO, BOLO$OR != "NCalc")
  #######################################################################################################################
  
  ######################################################## MAORI ######################################################## 
  # MODULE ADJUSTMENT ORI
  # calculates the OR difference for amino acid motifs by subtracting the
  # OR difference of the current iteration and the OR of the motif subset of the
  # previous iteration. this is known as OR.diff.A.
  # EX: if the current iteration is 1:2:3 K~L~Y, the motif subset for the previous
  # iteration would be 1:2 K~L. The OR difference would be 1:2:3 K~L~Y's OR - 
  # 1:2 K~L's OR.
  # for 3+mers, an additional OR difference is calculated for the current iteration OR
  # and the OR of the singleton residue added. this is known as OR.diff.B.
  # EX: if the current iteration is 1:2:3 K~L~Y, the OR difference would be 1:2:3
  # K~L~Y OR - 3 ~Y's OR.
  
  #creates dummy_KDLO for comparison to first BOLO ONLY on the 0th iteration
  if(counter==0){
    #makes dummy KDLO based on previous BOLO
    dummy_KDLO<-as.data.frame(t(c(Locus="TBA-loc",Allele="TBA-allele",OR=1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F)[rep(seq_len(nrow(as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F))), each=nrow(BOLO)),]
    dummy_KDLO[,1]<-BOLO$Locus
    dummy_KDLO[,2]<-BOLO$Allele
    
    BOLO<-getORdiff(BOLO, dummy_KDLO, counter)
  } 
  
  #subsets out binned alleles and any alleles with NA combinations
  if(counter>0){
    BOLO<-subset(BOLO, (BOLO$Allele!="binned") & (!grepl("NA", BOLO$Allele)) & BOLO$OR != "NCalc")
  }
  
  if(nrow(BOLO)==0){
    return(list(KDLO, BOLO, UMLO))
  }
  
  #get OR differences for iterations 1+
  #counter = 2+ will have OR.diff.B calculations
  if(counter != 0){
    BOLO<-getORdiff(BOLO, KDLO, counter, KDLO_list)
  }

  #subsets out NS values
  KDLO<-subset(BOLO,BOLO[,7]=="*")
  
  #### LOOP SPECIFICATIONS ####
  
  #filters out predisposing ORs for analysis
  if(loop==1){
    KDLO<-KDLO %>% filter(OR > 1.0)
  }
  
  #filters out protective ORs for analysis
  if(loop==2){
    KDLO<-KDLO %>% filter(OR <1.0)
  }
  
  if(nrow(KDLO)==1|nrow(KDLO)==0){
    return(list(KDLO, BOLO, UMLO="none"))
  }
  
  #statement for returning BOLO if KDLO=0
  #if((counter>0) & (nrow(KDLO)==0)){
  #  return(list(KDLO, BOLO, UMLO))
  #}
  
  #subsets out variants that have not shown an improvement based on the provided
  #OR difference threshold from their previous variants and
  #singular amino acids
  if(counter>1) {
    
    KDLO<-KDLO %>%
      mutate(OR.diff.B = round(OR.diff.B, 3)) %>%
      filter(OR.diff.B >0.1)
    
  }

  KDLO<-KDLO %>%
    mutate(OR.diff.A = round(OR.diff.A, 3)) %>%
    filter(OR.diff.A >0.1)
  
  #check for exon specific analysis for counter = 0
  #if no entries in KDLO make it past the OR difference 
  #threshold and counter is 0, return UMLO as none
  if(nrow(KDLO)==0){
    if(counter == 0){
      return(list(KDLO, BOLO, UMLO = "none"))
    } else{
      return(list(KDLO, BOLO, UMLO))
    }
  }
  
  #adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant
  KDLO<-BOLO %>%
    mutate(OR.diff.A = round(OR.diff.A, 3)) %>%
    filter(Locus %in% KDLO$Locus) %>%
    rbind.data.frame(KDLO) %>%
    distinct()

  if(!counter %in% c(0,1)){
     KDLO<-KDLO %>%
        mutate(OR.diff.B = round(OR.diff.B, 3))
  }
  
  #finds unassociated positions from current iteration
  unassociated_posi<-unique(BOLO$Locus[!BOLO$Locus %in% KDLO$Locus])
  
  #if length(unassociated_posi==0), return KDLO -- this means KDLO and BOLO are the same
  #and max improvement has been reached
  if(length(unassociated_posi)==0){
    return(list(KDLO, BOLO, UMLO))
  }
  
######################################################## NAME GENERATION ######################################################## 
  
  combinames<-NULL
  
  #pair name generation
  if(counter==0){
    start1<-unique(KDLO$Locus)
    
    #if nothing is in the KDLO, return KDLO and BOLO
    if((length(start1))==0){
      return(list(KDLO, BOLO))
    }
    
    for(i in 1:(length(start1)-1)){
      for(j in (i+1):length(start1)){
        combinames<-append(combinames, paste(start1[[i]], start1[[j]], sep = ':'))
      } 
    }
  }
  
  if(counter>0){
    start1<-unique(KDLO_list[[1]]$Locus)
  }
  
  if(counter>0){

    combinames <- foreach(i = 1:length(unique(KDLO$Locus)), .packages = c("gtools"), .combine = c) %dopar% {
      possibleCombiNameSplit <- strsplit(unique(KDLO$Locus)[[i]], ':')[[1]]
      residueAppend <- setdiff(start1, possibleCombiNameSplit)
      sapply(residueAppend, function(k) 
        paste(mixedsort(c(possibleCombiNameSplit, k)), collapse = ':')
      )
    }
    
    combinames<-unique(mixedsort(combinames))
  }
  
  ######################################################## SUBSET OUT UNASSOCIATED POSITIONS ######################################################## 
  
  #^ in grepl pattern indicates the beginning of a string
  #need to filter out unassociated positions where the new motif STARTS with the
  #unassociated position, or immediately after a colon
  
  # no unassociated positions for counter = 0 
  if(counter != 0){
    
    boolVector<-filterUnassociatedPositions(unassociated_posi, combinames)
    combinames<-combinames[!boolVector]
    
    if(counter==2) {
      
      boolVector_UMLO<-filterUnassociatedPositions(UMLO_list[[counter]], combinames)
      combinames<-combinames[!boolVector_UMLO]
    }
    
    if (counter > 2) {

      boolVector_UMLO_2<-filterUnassociatedPositions(unlist(UMLO_list[2:counter]), combinames)
      combinames<-combinames[!boolVector_UMLO_2]
      
    }
    
    #end if no combo names are generated after filtering out unassociated positions
    if(length(combinames)==0) {
      return(list(KDLO, BOLO, UMLO))
    }
  }
  
  ######################################################## CREATE NEW DF FOR NEXT ITERATION EVAL ######################################################## 
  
  #df for pairs -- length is number of unique pairs * 2,
  combidf<-data.frame(variantAAtable[[loci]][,c(1,2)], matrix("", ncol =length(rep(combinames, 2))), stringsAsFactors = F)
  
  temp<-NULL
  for(i in 1:length(combinames)){
    temp[[i]]<-rep(combinames[i], 2)  
  }
  
  #fills in column names
  colnames(combidf)<-c("SampleID", "Disease", unlist(temp))

  #observes number of columns for those needed to be pasted together
  cols=c(1:length(strsplit(combinames[[1]], ":")[[1]]))
  
  #[[1]] to contain amino acid combos of TRUE/FALSE
  #[[2]] to contain amino acid combos of FALSE/TRUE
  dfAA<-sapply(1:2, function(x) NULL)
  
  #fills in element names in the lists formed in the above lists
  for(j in 1:length(dfAA)){
    dfAA[[j]]<-sapply(combinames, function(x) NULL)
  }
  
  combiCols <- generateAAcols(combinames, variantAAtable[[loci]], cols)

  #fills into pair_df
  combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-  sapply(combiCols, '[[', 1)
  combidf[,3:length(combidf)][,c(FALSE,TRUE)]<- sapply(combiCols, '[[', 2)
  
  #saves each iteration into specified elements in a list in a variable "myData"
  #returns myData
  myDataFinal<-list("KDLO"=KDLO, "BOLO"=BOLO, "combidf"=combidf, "UMLO"=unassociated_posi, "combinames"=combinames)
  
  return(myDataFinal)
}


#wrapper function for creating new motifs
generateAAcols<-function(combinations, aaTable, columns){
  
  aaCols <- foreach(i = 1:length(combinations), .packages = c("dplyr", "tidyr")) %dopar% {
    
    left_vec <- aaTable[c(TRUE, FALSE)][strsplit(combinations, ":")[[i]]][, columns] %>%
      unite(newNames, sep = '~') %>% pull(newNames)
    
    right_vec <- aaTable[c(FALSE, TRUE)][strsplit(combinations, ":")[[i]]][, columns] %>%
      unite(newNames, sep = '~') %>% pull(newNames)
    
    list(left_vec = left_vec, right_vec = right_vec)
  }
  
  return(aaCols)
}

