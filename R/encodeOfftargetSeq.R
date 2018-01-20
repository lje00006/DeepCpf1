#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DeepCpf1
# 2017
# Jiesi Luo, Wake Forest School of Medicine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Encode Cpf1 guide RNA mismatched sequences
#' 
#' Sequence is encoded as twelve binary features (i.e., 100000000000  for 
#' AT, 010000000000 for AC). The output of this function can be used
#' directly for deep learning method.
#'
#'
#' @usage encodeOfftargetSeq(SeqFileName)
#'
#' @param SeqFileName A character name of the input Cpf1 guide RNA on-target #' and off-target sequences file,including full path to file if it is located #' outside the current working directory.
#'
#'@return test_offtarget.csv A matrix containing encoded features. Sequence 
#'features are represented as binary numbers.
#'
#'@author Jiesi Luo
#¡¯
#'@examples
#'encodeOfftargetSeq("mismatch.bed")
#'@export output


encodeOfftargetSeq <- function (SeqFileName) {

     library(stringr)
    
     Seq<-read.table(SeqFileName,header=TRUE)

     offtarget_sites<-Seq$off

     ontarget_site<-Seq$on

     one_hot1<-diag(1,12,12)

     one_hot2<-matrix(0,nrow=4,ncol=12)

     one_hot<-rbind(one_hot1,one_hot2)

       dinuc_type<-c("AT","AC","AG","TA","TC","TG","CA","CT","CG","GA","GT","GC","AA","TT","CC","GG")


       test<-vector()

        matchpairs<-vector()

        on<-substring(ontarget_site,1:27,1:27)


         for (num in 1:length(offtarget_sites)) {

             off<-substring(offtarget_sites[num],1:27,1:27)

             for (posi in 1:27) {
     
            matchpairs[posi]<-str_c(on[posi],off[posi],collapse="")
                            }

        for (dinucle in 1:16) {

            matchpairs<-str_replace_all(matchpairs, dinuc_type[dinucle],str_c(as.character(one_hot[dinucle,]),collapse =""))

                              }

   test[num]<-str_split(str_c(matchpairs,collapse =""),"")

}

test_array<-t(do.call(rbind,test))

output<-rbind(as.character(ontarget_site),as.character(offtarget_sites),test_array)

write.csv(t(output),file="test_offtarget.csv",row.names=FALSE)

}