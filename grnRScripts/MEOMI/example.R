source("manipulation_function.R")
source("Conditional_interaction_information_calculation.R")
source("MEOMI.R")
##Set the path to read the file and the path to output the file
path_data <- "../PREMER/insilico_size100_1_knockouts.csv"
#  path_standard <- "D:/DREAM3GoldStandard_InSilicoSize10_Ecoli1.txt"
Path_out <- "result.csv" 
  ##Read standard network and data
mydata <- read.csv(file = path_data)
#print(mydata)
#net <- read.table(file = path_standard,header = F, sep = "\t")
##Delete the extra rows
#if(colnames(mydata)[1] != "G1")
#{
#  mydata<-mydata[-1,-1]
#} 
print(mydata)
  ##The standard network is processed to change direction to undirection
#  for(m in 1: nrow(net)){      
#    if(net[m,]$V3==1){
#      p <- net[m,]$V1
#      q <- net[m,]$V2
#      if(net[(net$V1==q&net$V2==p),]$V3 == 0){
#        net[(net$V1==q&net$V2==p),]$V3 <- 1
#      }
#    }
#  }
  ##Set up parameters
lisan <- 5
lambda <- 0.1
order <- 4
Gval<-MEOMI(mydata = mydata, bins = lisan,lamda = lambda, order = order)
#r <- canshu(mydata, net, Gval[,c(1:3)])
put_data(Path_out, Gval[,c(1:3)])#, data[i])

  
  
  
  
  
  
 
