require(nets)

path="C:/Users/Swammy/Documents/GitHub/R-Package-nets/demo/data/sp100.csv"

SP=read.csv(path)

#SP[,1] is the column of dates, others columns are stock values
T <- length(SP[,1])
N <- length(SP[1,])-1

#Put the data in a matrix so it can be used by the function nets
y <- matrix(0,T,N)

networks <-list()

for (i in 1:N)
	{y[,i]=SP[,i+1]}

chrono <- matrix (0,45,1)
for (i in 41:45)
{
yy <- y[,1:i]
chrono[i] <- proc.time()[3]
network <- nets( yy, type='pc', lambda=1000)
chrono[i] <- proc.time()[3]-chrono[i]
networks = c(networks,network)
}

write.csv(chrono, file = "time1045.csv")
