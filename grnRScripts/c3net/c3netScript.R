#sudo add-apt-repository ppa:jonathonf/gcc-7.1
#sudo apt-get update

#sudo apt-get install gcc-7 g++-7
#sudo apt-get install gfortran-7

#install.packages("tcltk")
#install.packages("igraph")
#install.packages("c3net")

library('c3net')

data<- read.csv("../PREMER/insilico_size100_1_knockouts.csv") # apenas mudar o arquivo de entrada (testado com dream 4)

net <- c3net(data)#, network=TRUE #caso queira visualizar a rede
write.csv(net, file="rede_matrix_c3net.csv")
