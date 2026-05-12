#### BIGDAWG DATASET SUMMARY FOR PREDISPOSING AND PROTECTIVE ANALYSES#####
# v 0.9
BIDS <- function(loci, dataset) {
  
  #cluster initialization for parallelization processes
  cl <- initCluster()
  
  # source cluster helper functions to all workers
  clusterEvalQ(cl, {
    source('scripts/clusterHelpers.R')
  })
  
  input_DS <- read.table(
    dataset,
    header = TRUE,
    sep = "\t",
    quote = "",
    as.is = TRUE,
    colClasses = "character",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(input_DS) <- gsub("\\_1|\\_2|\\.1|\\.2", "", colnames(input_DS)) ### SJM ADDED 8.27.2021 -- strip out locus suffix numbers
  
  BIGCAAT_results <- BIGCAAT(loci, dataset)
  
  if (is.list(BIGCAAT_results) == FALSE) {
    return(BIGCAAT_results)
  }
  
  BIDS_start <- proc.time()
  
  for (j in 1:length(loci)) {
    
      locus_start<-proc.time()
      
      locus <- loci[[j]]
      
      cat(sprintf('Summarizing BIGCAAT results for HLA-%s... \n', locus))
      
      # find locus specific alleles in input data
      loci_columns <- input_DS[colnames(input_DS) == locus]
      all_alleles <- unique(unlist(c(loci_columns), use.names = F))
      ds_alleles <- sort(paste(locus, all_alleles[all_alleles != ''], sep = '*'))
      
      # build locus alignment to pass to findMotif
      alignments <- buildAlignments(locus, 'AA')[[1]]$AA
      
      summaries <- c('Predisposing', 'Protective')
      
      locus_results <- BIGCAAT_results[[locus]]

      results <- foreach(k = summaries,
                         .packages = c('dplyr', 'tibble', 'gtools')) %dopar% {
                           
                           tryCatch({
                             summary_results <- locus_results[[k]]

                             #get max nmers from KDLO
                             max_nmer <- length(summary_results[[1]])
                             
                             if (length(summary_results$KDLO) == 0) {
                               return(list(k = k, result = 'No analytical results available'))
                             }
                             
                             else{
                               summary_rows <- summary_results$KDLO[[max_nmer]] %>%
                                 arrange(desc(OR)) %>%
                                 filter(sig == '*') %>%
                                 filter(OR < 1)
                               
                               result_list <- list()
                               
                               summary_rows <- summary_rows %>%
                                 add_column(motif = '') %>%
                                 add_column(alleles = '')
                               
                               # find which dataset alleles contain the significant motifs and add to summary
                               for (i in 1:nrow(summary_rows)) {
                                 motif <- paste(locus, paste(paste(
                                   strsplit(summary_rows[i, ]$Locus, ':')[[1]],
                                   strsplit(summary_rows[i, ]$Allele, '~')[[1]],
                                   sep = ""
                                 ), collapse = '~'), sep = '*')
                                 alleles <- findMotif(motif, alignments)$trimmed_allele
                                 ds_allele <- list(c(unique(alleles[alleles %in% ds_alleles])))
                                 
                                 summary_rows[i, ]$motif <- motif
                                 summary_rows[i, ]$alleles <- ds_allele
                               }
                               
                               result_list <- finalReport(summary_rows)
                               
                               return(list(k = k, result = result_list))

                           }},
                             error = function(e) {
                              return(list(k=k, result = sprintf("Summary generation failed due to the following error: %s", e)))
                             })
                         }
      
      
      names(results) <- sapply(results, function(x)
        x$k)
      
      BIGCAAT_results[[locus]][['Predisposing']]$`Predisposing Summary` <- results[['Predisposing']]$result
      BIGCAAT_results[[locus]][['Protective']]$`Protective Summary` <- results[['Protective']]$result
      
      cat(sprintf('elapsed time for BIDS HLA-%s %s seconds \n', locus, round((
        proc.time() - BIDS_start
      )['elapsed'])))
    }
    
    cat(sprintf('TOTAL elapsed time for BIDS %s seconds', round((
      proc.time() - BIDS_start
    )['elapsed'])))
  
  endCluster(cl)
  
  return(BIGCAAT_results)
}
