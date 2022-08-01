#' QSeq function
#' 
#' A wrapper for QScale that takes phyloseq object as input. Metadata needs to have sample column with scaling factor
#' @param ps = phyloseq object with community data, and sample data with at least one column for scaling
#' @param abundance = name of column in metadata to scale counts by
#' @keywords QSeq, Quantitative Sequencing, Quantitative Scaling
#' @export
#' @examples QSeq(ps, abundance="QPCR_16S")
#' QSeq()

QSeq<-function(ps, abundance){
  
  ps<-physloseq::prune_samples(!is.na(sample_data(ps)[[abundance]]), ps)
  
  print(paste(sum(is.na(physloseq::sample_data(ps)[[abundance]])), " sample(s) missing abundance abundance data. Removing ", sum(is.na(physloseq::sample_data(ps)[[abundance]])), " samples.")) # test to see if works []
  
  scale<-as.numeric(as.character(physloseq::sample_data(ps)[[abundance]])) # make sure scaling data is not a factor
  
  out<-Qscale(ps, scale)
  
  out
}


#' Function that scales otu table by column in metadata
#' 
#' Takes a phyloseq object, and vector for scaling samples in otu table. After scaling, rounds taxon count values to the nearest integer, then returns a phyloseq object with scaled community data.
#' @param ps = phyloseq object with community data, and sample data with at least one column for scaling
#' @param scale = vector that is used for scaling samples
#' @keywords QSeq, Quantitative Sequencing, Quantitative Scaling
#' @export
#' @examples 
#' QScale()
Qscale<-function(ps, scale){
  otu <- data.frame(as.matrix(physloseq::otu_table(physloseq::transform_sample_counts(ps, function(x) x/sum(x))))) # transform count data into relative abundance, output to data table
  if(!taxa_are_rows(ps)){
    otu<-data.frame(t(as.matrix(physloseq::otu_table(physloseq::transform_sample_counts(ps, function(x) x/sum(x))))))
    } # ensure taxa are rows
  
  #scaled<-mapply(`*`, otu, scale)
  #colnames(scaled)<-sample_names(ps) 
  #rownames(scaled)<-taxa_names(ps)
  #scaled<-round(scaled)
  
  otu[] <- mapply(`*`, otu, scale) # Multiply relative abundance by total abundance (gene)
  scaled<-round(otu) # round taxon abundance to nearest integer
  p2<-ps
  physloseq::otu_table(p2) <- physloseq::otu_table(scaled, taxa_are_rows=T)
  
  return(p2)
  
}