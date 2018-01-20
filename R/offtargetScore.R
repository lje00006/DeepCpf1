#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DeepCpf1
# 2017
# Jiesi Luo, Wake Forest School of Medicine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' Predict Cpf1 guide RNA off-target specificity using convolutional deep 
#' learning neural networks
#'
#' Training input of CNN classifier is "one-hot" matrix of size 12 x 27.The 
#' first convolutional layer performs 35 convolutions with 7 x 7 filter on 
#' the input layer, producing 35 feature maps of size 6 x 21. The second 
#' pooling layer performs 2 x 2 spatial pooling of each feature map using the #' sum value, producing 35 feature maps of size 3 x 10. A fully connected 
#' layer contains 30 neurons. Finally, the output of the fully connected layer #' is fed to a linear regression layer the assigns a score for the 
#' specificity.
#'
#'
#' @usage offtargetScore(trainFileName,testFileName)
#'
#' @param trainFileName A character name of the training "one-hot" matrix file
#' @param testFileName  A character name of the test "one-hot" matrix file
#'
#'@return result_offtarget.csv A file containing specificity scores
#'@author Jiesi Luo
#¡¯
#'@examples
#'offtargetScore("train_offtarget.csv","test_offtarget.csv")
#'@export result_offtarget
   


offtargetScore<- function (trainFileName,testFileName) {

   library("mxnet")

   train<-read.csv(trainFileName)

   test<-read.csv(testFileName)

   train<-data.matrix(train)

   train_x<-t(train[,-1])

   train_y<-train[,1]

   dim(train_x)<-c(12,27,1,ncol(train_x))

   test_x<-t(test[,3:326])

   on_target<-as.character(test[,1])

   off_targets<-as.character(test[,2])

   dim(test_x)<-c(12,27,1,ncol(test_x))

   data <- mx.symbol.Variable('data')

   conv_1 <- mx.symbol.Convolution(data = data, kernel = c(7,7), num_filter = 35)

   relu_1 <- mx.symbol.Activation(data = conv_1, act_type = "relu")

   pool_1 <- mx.symbol.Pooling(data = relu_1, pool_type = "sum", kernel = c(2, 2), stride = c(2, 2))

   flatten <- mx.symbol.Flatten(data = pool_1)
 
   fc_1 <- mx.symbol.FullyConnected(data = flatten, num_hidden = 300)
 
   relu_2 <- mx.symbol.Activation(data = fc_1, act_type = "relu")
 
   fc_2 <- mx.symbol.FullyConnected(data = relu_2, num_hidden = 1)

   NN_model <- mx.symbol.LinearRegressionOutput(data = fc_2)

   mx.set.seed(100)
 
   devices <- mx.cpu()

   model <- mx.model.FeedForward.create(NN_model,X = train_x, y = train_y, ctx = devices, num.round = 200, array.batch.size = 40,learning.rate = 0.005,momentum = 0.9, eval.metric = mx.metric.rmse, epoch.end.callback = mx.callback.log.train.metric(100))

  predicted<-predict(model,test_x)

  result_offtarget<-rbind(on_target,off_targets,predicted)

   write.csv(t(result_offtarget),file="result_offtarget.csv")

}

