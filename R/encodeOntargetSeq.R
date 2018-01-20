#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DeepCpf1
# 2017
# Jiesi Luo, Wake Forest School of Medicine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' Encode Cpf1 guide RNA targetable sequence
#' 
#' Sequence is encoded as sixteen binary features (i.e., 1000000000000000  for #' TT, 0100000000000000 for TC). The output of this function can be used
#' directly for deep learning method.
#'
#'
#'
#'
#' @usage encodeOntargetSeq(fastaFileName)
#'
#' @param fastaFileName A character name of the input Cpf1 guide RNAs fasta 
#' format file,including full path to file if it is located outside the 
#' current working directory.
#'
#'@return test_ontarget.csv A matrix containing encoded features. Sequence 
#'features are represented as binary numbers.
#'
#'@author Jiesi Luo
#¡¯
#'@examples
#'encodeOntargetSeq("Sequence.bed")
#'@export output


encodeOntargetSeq<- function(fastaFileName) {
  
  library(seqinr)
  library(stringr)
  
  DNAsequence<-read.fasta(fastaFileName,seqtype=c("DNA"),as.string=TRUE,
                          forceDNAtolower=FALSE,seqonly=TRUE)

one_hot<-diag(1,16,16)
  
dinuc_type<-c("TT","TC","TA","TG","AT","AC","AA","AG","CT","CC","CA","CG","GT",
                "GC","GA","GG")
  test<-vector()
  for (num in 1:length(DNAsequence)) {
     dinucleotides<-substring(DNAsequence[num],1:26,2:27)
     for (dinuc in 1:16) {
            dinucleotides<-str_replace_all(dinucleotides,dinuc_type[dinuc],
                                           str_c(as.character(one_hot[dinuc,]),collapse =""))
                         }
     test[num]<-str_split(str_c(dinucleotides,collapse =""),"")
  }
  
test_array <-t(do.call(rbind,test))

output<-rbind(DNAsequence,test_array)

write.csv(t(output),file="test_ontarget.csv",row.names=FALSE)

}



