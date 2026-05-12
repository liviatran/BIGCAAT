##### BIGCAAT #####
# v 0.9
BIGCAAT <- function(loci, GenotypeFile) {
  start<-proc.time()
  if (missing(loci)) {
    return(cat("Please specify a locus, or vector of loci to analyze.")) 
  }
  
  if (missing(GenotypeFile)) {
    GenotypeFile <- fileChoose("Please select a BIGDAWG-formatted genotype datset for analysis.")
  }
  
  cat("-------------------------------------------------------------------\n BIGCAAT: BIGDAWG Integrated Genotype Converted Amino Acid Testing\n-------------------------------------------------------------------\n") ### SJM Banner
  Genotype_Data <- read.table(GenotypeFile, header = TRUE, sep = "\t", quote = "", na.strings = c("", NA), colClasses = "character", check.names = FALSE)
  
  colnames(Genotype_Data) <- gsub("\\_1|\\_2|\\.1|\\.2","",colnames(Genotype_Data)) ### SJM ADDED 8.27.2021 -- strip out locus suffix numbers
  
  AAData <- variantAAextractor(loci, Genotype_Data, NULL)
  
  cat("\nGenotypes Converted to Amino Acid Sequences.\n\n") ### SJM ADDED 8.27.2021
  if(is.list(AAData)==FALSE){
    return(AAData)
  }
  
  CombiData <- vector("list",length(loci))
  names(CombiData) <- loci
  
  for(z in 1:length(CombiData)){
    
    #specifications for predisposing and protective OR analysis added by LT
    CombiData[[loci[[z]]]] <- sapply(c("Predisposing", "Protective"), function(x) NULL)
  
  }

  for(loop in 1:length(CombiData[[loci[[z]]]])){
    
    
    if(loop==1){
      analysis_type <- 'Predisposing'
    } else{
      analysis_type <- 'Protective'
    }
    
    cat(sprintf('%s OR analysis \n', analysis_type))
    
    for (p in 1:length(loci)) {
      
      continue <- TRUE
      
      cat("Analyzing the HLA-",loci[p]," locus \n",sep="") ### SJM added notification
      
      # if an error occurs in a loop for a given locus, populate the list 
      # element with the error message and move onto the next loop or locus
      # analysis 
      tryCatch(

        CombiData[[loci[p]]][[loop]] <- runCombiAnalyzer(loci[p], AAData, loop), 
        
        error = function(e){
           cat(sprintf('%s analysis for HLA-%s failed due to the following error: %s', 
                                                 analysis_type, loci[[p]], e))
          continue<<- FALSE
        }
      )
        
        if(!continue){
          next
        }
      
    }
  }
  
  print(sprintf('TOTAL elapsed time for combiAnalyzer %s seconds', round((proc.time()-start)['elapsed'])))
  
  #}
  CombiData
}


