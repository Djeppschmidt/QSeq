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
  # check to make sure column exists for abundance data
  if(!any(names(phyloseq::sample_data(ps))==abundance)){
    stop("Error: abundance column does not exist in sample data. Metadata has the following columns: ", names(phyloseq::sample_data(ps)), " ; You asked for: ", abundance)
  }
  
  ps<-phyloseq::prune_samples(!is.na(sample_data(ps)[[abundance]]), ps)
  
  print(paste(sum(is.na(phyloseq::sample_data(ps)[[abundance]])), " sample(s) missing abundance abundance data. Removing ", sum(is.na(phyloseq::sample_data(ps)[[abundance]])), " samples."))
  # test to see if works []
  
  scale<-as.numeric(as.character(phyloseq::sample_data(ps)[[abundance]])) # make sure scaling data is not a factor
  
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
  otu <- as.data.frame(as.matrix(phyloseq::otu_table(phyloseq::transform_sample_counts(ps, function(x) x/sum(x))))) # transform count data into relative abundance, output to data table
  if(!taxa_are_rows(ps)){
    otu<-as.data.frame(t(as.matrix(phyloseq::otu_table(phyloseq::transform_sample_counts(ps, function(x) x/sum(x))))))
    } # ensure taxa are rows
  
  #scaled<-mapply(`*`, otu, scale)
  #colnames(scaled)<-sample_names(ps) 
  #rownames(scaled)<-taxa_names(ps)
  #scaled<-round(scaled)
  
  otu[] <- mapply(`*`, otu, scale) # Multiply relative abundance by total abundance (gene)
  print(rownames(otu))
  scaled<-round(otu) # round taxon abundance to nearest integer
  
  #print(colnames(scaled))
  p2<-ps
  phyloseq::otu_table(p2) <- phyloseq::otu_table(scaled, taxa_are_rows=T)
  
  return(p2)
  
}

#' Import from QIIME2 function
#' 
#' A function that takes output from QIIME2 pipeline and outputs a phyloseq object. Assumes that DADA2 has been run in QIIME2, but you will assign taxonomy in R.
#' @param features = (required) file output from QIIME2 implementation of DADA2, as in "table.qza"
#' @param metadata = (required) text file (and path to file) with metadata to add to phyloseq object. Column with sample identifiers must be labeled "Sample_ID";
#' @param ASV = (required) QIIME2 file with dada2 ASV output. As in: "representative_sequences.qza"
#' @param taxref = (required) path to the training file for taxonomy assignment. We recommend using SILVA.
#' @keywords QIIME2, phyloseq
#' @export
#' @examples QSeq(ps, abundance="QPCR_16S")
#' getQIIME()

getQIIME<-function(features, metadata, ASV, taxref){
  
  # convert QIIME 2 data to phyloseq object.
  # this object does not have sequences, or taxonomy table
  ps<-qiime2R::qza_to_phyloseq(features=features, metadata=metadata)
  
  # get sequences for ps object
  tax<-qiime2R::read_qza(file=ASV)
  
  
  # replace taxa names with sequences
  tnames<-phyloseq::taxa_names(ps)
  for(i in tnames){
    tnames[tnames==i]<-as.character(tax$data[[i]])
  }
  
  phyloseq::taxa_names(ps)<-tnames
  
  # assign taxonomy to sequences:
  taxa <- phyloseq::assignTaxonomy(tnames, taxref, multithread=TRUE) #SILVA V138 / August 15 2020
  
  # add taxonomy table to phyoseq object
  phyloseq::tax_table(ps)<-phyloseq::tax_table(taxa)
  
  # output complete ps object with sequences and metadata for taxonomy and sample data
  return(ps)
}


# Future work:
# abundance-occupancy modeling (abundance vs occupancy)
# correlation (covariance and co-occurance)
# identification of core microbiome using quantitative methods


