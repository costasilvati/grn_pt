library('ennet')

data<- read.csv("../PREMER/insilico_size100_1_knockouts.csv") # apenas mudar o arquivo de entrada (testado com dream 4)

net <-  (data)#, network=TRUE #caso queira visualizar a rede
write.csv(net, file="rede_matrix_ennet.csv")
