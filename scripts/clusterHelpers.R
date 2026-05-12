# cluster initiation for parallelization
# v 0.9
initCluster<-function(){
  
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  return(cl)
  
}

endCluster<-function(cluster){
  stopCluster(cluster)
}

# generates a collapsed report of unique motifs and unique alleles from the BIDS
# $final data frame in the BIDS return object 
# grouped by p-value, OR, and confidence intervals
finalReport <-function(finalFrame) {
  
  p_val_summary<-groupUniquePVals(finalFrame)
  
  # Generates a collapsed report of unique motifs and unique alleles from the BIDS $final data frame in the BIDS return object 
  newFrame <- data.frame(matrix(ncol=4,nrow=nrow(p_val_summary)))
  
  for(i in 1:nrow(newFrame)) {
    
    newFrame[i,1] <-paste(unique(unlist(p_val_summary$motif[i])),collapse=",")
    newFrame[i,2] <- length(unique(unlist(p_val_summary$motif[i])))
    newFrame[i,3] <-paste(unique(unlist(p_val_summary$DSalleles[i])),collapse=",")
    if(substr(newFrame[i,3],1,1) == ",") {newFrame[i,3] <- substr(newFrame[i,3],2,nchar(newFrame[i,3]))}
    newFrame[i,4] <- length(unique(unlist(p_val_summary$DSalleles[i])))
  }
  
  collapsedReport <- cbind(newFrame,p_val_summary[,3:6])
  colnames(collapsedReport)[1:4] <- c("Motif","#Motifs", "Shared_Alleles", "#Alleles")
  
  collapsedReport
  
}

# modification of findMotif from SSHAARP without requiring file name
# also accepts alignment instead of creating in the function itself
findMotif<-function(motif, HLAalignment){
  
  motifs<-SSHAARP::getVariantInfo(motif)[[2]]
  
  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)
  
  for(x in 1:length(motifs)) {
    HLAalignment <- HLAalignment[HLAalignment[substr(motifs[x],1,nchar(motifs[x])-1)]==substr(motifs[x],nchar(motifs[x]),nchar(motifs[x])),]
  }
  
  if(nrow(HLAalignment)==0){
    return(c(FALSE, paste(motif, "No alleles possess this motif", sep=": ")))
  }
  #if motifs are found, HLAalignment[[loci[[i]]]] is returned
  return(HLAalignment)
}

# generates a summary of each unique p value, all motifs that share that p value,
# and all dataset alleles that contain the motif
groupUniquePVals <- function(df){
  
  unique_p_vals<-df %>% 
    group_by(p.value) %>% 
    filter(n()>1) %>% 
    ungroup()
  
  p_val_summary<-unique_p_vals %>%
    group_by(p.value) %>%
    summarise(motif = list(motif), DSalleles = unique(alleles), OR = list(unique(OR)), CI.lower = list(unique(CI.lower)), CI.upper = list(unique(CI.upper))) %>%
    relocate(p.value, .after = CI.upper)
  
  return(p_val_summary)
}

# wrapper function for UMLO
filterUnassociatedPositions<-function(positions, combinationNames){
  
  # return a boolean for combination names that start with the unassociated
  # motif or that contain the unassociated motif
  # (^|:) --> starts with the unassociated position(s) or has a colon before the
  # unassociated positions. 
  # EX: -17:4 would match -17:4:30 and -24:-17:4
  # (:|$) --> matches any combination names with a colon after the motif, or
  # indicates the end of the string
  # EX: -17:4 would match -24:-17:4:30 and -24:-17:4, but NOT -24:17:47
  
  res <- foreach(i = 1:length(positions)) %dopar% {
    grepl(sprintf("(^|:)%s(:|$)", positions[[i]]), combinationNames)
  }
  
  containMatch <- Reduce('|', res)
  
  return(containMatch)
}

# wrapper function for MAORI 
getORdiff<-function(BOLOdf, KDLOdf, count, KDLO_list="") {
  #set up BOLO to join with KDLO for the previous iteration's OR value, as well as 
  #the last residue's OR value from the very first KDLO
  #get previous motif positions and motif (n-1)
  #get last residue position and residue 
  BOLOdf<-BOLOdf %>%
    mutate(temp1= strsplit(Locus, ':'),
           prevMotifPositions = sapply(temp1, function(x) paste(x[1:count], collapse = ":")),
           temp2=strsplit(Allele, '~'),
           prevMotif = sapply(temp2, function(x) paste(x[1:count], collapse = "~")),
           lastPosition = sapply(temp1, function(x) x[count+1]),
           lastResidue = sapply(temp2, function(x) x[count+1])) %>%
    select(-c(temp1, temp2))
  
  #OR.diff.A calculation - get difference between current iteration and last iteration OR
  #join with last iteration KDLO by previous motif position and allele name to get OR -- find the difference between the two
  #to populate OR.diff.A (difference in ORs between current and previous iteration)
  BOLOdf<-BOLOdf %>%
    left_join(KDLOdf[,c('Locus', 'Allele', 'OR')], by = join_by(x$prevMotifPositions == y$Locus, x$prevMotif == y$Allele)) %>%
    mutate(OR.x=as.numeric(OR.x),
           OR.y=as.numeric(OR.y)) %>%
    mutate(OR.diff.A = OR.x - OR.y) %>% 
    select(-c(OR.y, prevMotifPositions, prevMotif)) %>%
    rename('OR' = OR.x)
  
  #remove negatives from any OR diffs for OR.diff.A and OR.diff.B - OR diff is
  #should be based on absolute value
  #there are some instances where NA may be present in OR.diff.A or OR.diff.B
  #for OR.diff.A, an nmer may be significant, but the specific residue may not have
  #been significant or moved on. these residues won't be found in the previous
  #iteration for comparison, so the join will return NA. replcae these with 0.
  #EX -- MS_EUR_med.txt dataset for DRB1
  #in evaluation for 3mers, -24:11:26 L~V~L, L~V is not found in the previous KDLO
  #this was because it was filtered out in the previous BOLO -- the motif was binned
  BOLOdf<-BOLOdf %>%
    mutate(OR.diff.A = case_when(grepl('-', .$OR.diff.A) ~ as.numeric(gsub('-', '', OR.diff.A)),
                                 is.na(OR.diff.A) ~ 0,
                                 .default = as.numeric(OR.diff.A)))
  
  #OR.diff.B calculation is not applicable for 0th or first iteration
  if(!count %in% c(0,1)){
    #OR.diff.B calculation - get difference between current iteration and the first iteration for the last residue
    #EX: motif being evaluated -- -24:-7:-1 - L~A~S. find OR value for -1:S
    #join with the first iteration KDLO by last position and residue
    BOLOdf<-BOLOdf %>%
      left_join(KDLO_list[[1]][,c('Locus', 'Allele', 'OR')], by = join_by(x$lastPosition == y$Locus, x$lastResidue == y$Allele))%>%
      mutate(OR.y=as.numeric(OR.y)) %>%
      mutate(OR.diff.B = OR.x - OR.y) %>%
      select(-c(OR.y, lastPosition, lastResidue)) %>%
      rename('OR' = OR.x)
    
    
    BOLOdf<-BOLOdf %>%
      mutate(OR.diff.B = case_when(grepl('-', .$OR.diff.B) ~ as.numeric(gsub('-', '', OR.diff.B)),
                                   is.na(OR.diff.B) ~ 0,
                                   .default = as.numeric(OR.diff.B)))
  }
  
  return(BOLOdf)
  
}
