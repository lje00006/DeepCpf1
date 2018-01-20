#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DeepCpf1
# 2017
# Jiesi Luo, Wake Forest School of Medicine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' Encode Cpf1 targetable DNA sequence from the a genomic sequence
#' 
#' Targetable sites for Cpf1 are identified by searching for genomic sequence #' matching TTTN-N23 motif. Then, the targetable sites that contained polyT
#' and extreme GC content (<30% or >70%) were removed. Sequence is encoded as #' sixteen binary features (i.e., 1000000000000000 for TT, 0100000000000000 
#' for TC). The output of this function can be used directly for deep learning #' method.
#'
#'
#'
#'
#' @usage encodeOntargetGeno(fastaFileName)
#'
#' @param fastaFileName A character name of the input a gene fasta format 
#' file,including full path to file if it is located outside the current 
#' working directory.
#'
#'@return test_ontarget.csv A matrix containing encoded features. Sequence 
#'features are represented as binary numbers.
#'
#'@author Jiesi Luo
#¡¯
#'@examples
#'encodeOntargetGeno("gene.bed")
#'@export output


encodeOntargetGeno<- function(fastaFileName) {
  
  library(seqinr)
  library(stringr)
  library(IRanges)

  DNAsequence<-read.fasta(fastaFileName,seqtype=c("DNA"),as.string=TRUE,
                          forceDNAtolower=FALSE,seqonly=TRUE)

  for (num in 1:length(DNAsequence)) {
      segment_27mer<-substring(DNAsequence[[num]],1:(nchar(DNAsequence[[num]])
                                                     -26),27:nchar(DNAsequence[[num]]))
  }

  pos_strand<-na.omit(str_extract(segment_27mer,"^TTT[ATCG]{24}"))
  neg_strand<-na.omit(str_extract(segment_27mer,"[ATCG]{24}AAA$"))
  neg_strand<-reverse(neg_strand)
  neg_strand<-chartr("ATCG","TAGC",neg_strand)
  targetsites_all<-c(pos_strand,neg_strand)

  guideRNA<-substr(targetsites_all,5,27)
  GC<-(str_count(guideRNA,"C")+str_count(guideRNA,"G"))/nchar(guideRNA)*100
  polyT<-str_count(guideRNA,"T{4,}")
  targetsites_filter<-targetsites_all[GC>30 & GC<70 & polyT<1]
  
one_hot<-diag(1,16,16)
  
dinuc_type<-c("TT","TC","TA","TG","AT","AC","AA","AG","CT","CC","CA","CG","GT",
                "GC","GA","GG")
  test<-vector()
  for (num in 1:length(targetsites_filter)) {
     dinucleotides<-substring(targetsites_filter[num],1:26,2:27)
     for (dinuc in 1:16) {
            dinucleotides<-str_replace_all(dinucleotides,dinuc_type[dinuc],
                                           str_c(as.character(one_hot[dinuc,]),collapse =""))
                         }
     test[num]<-str_split(str_c(dinucleotides,collapse =""),"")
  }
  
test_array <-t(do.call(rbind,test))

output<-rbind(targetsites_filter,test_array)

write.csv(t(output),file="test_ontarget.csv",row.names=FALSE)

}



