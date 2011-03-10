library(MCMCpack)
library(ggplot2)

# There are 47 terms(1953-1999), 3450 cases and 29 Justices.
T <- 47;
K <- 3450;
J <- 29;

# big.data is the matrix whose rows correspond to each case and columns correspond to each justice, and value of each element represents: 1=reverse, 0=affirm, NA=not on that case 
big.data <- matrix(NA, K, J)
colnames(big.data) <- c('Harlan' , 'Black', 'Douglas', 'Stwart',
     'Marshall', 'Brennan', 'White', 'Warren', 'Clark', 'Frankfurter', 'Whittaker', 'Burton', 'Reed', 'Fortas', 'Goldberg', 'Minton', 'Jackson', 'Burger', 'Blackmun', 'Powell', 'Rehnquist', 'Stevens', 'Oconnor', 'Scalia', 'Kennedy', 'Souter', 'Thomas', 'Ginsberg', 'Breyer');
big.data <- data.frame(big.data);


# read data from each year's txt file
vote <- list();
case.term <- rep(NA, T);
for(t in 1:T){
  path = paste('read.table(\'data/vote/', t+1952,'.txt\')', sep='');
  vote.term <- eval(parse(text=path))
  vote[[t]] <- cbind(vote.term, term=t+1952, time=t);
  case.term[t] <- dim(vote.term)[1];
}

# read adjacency matrix that contains information about which justic sits on which term
adj <- read.table("data/adj.txt")

# construct big.data
term <- rep(NA, K);
time <- rep(NA, K);
case.count <- 0;
for(i in 1:T){
  for(j in 1:case.term[i]){
    case.count <- case.count + 1;
    num.justice <- 9;
    if(i==4 | i==9) num.justice <- 10;
    if(i==17) num.justice <- 8;
    for(k in 1:num.justice){
      big.data[case.count,which(adj[i,]==(k-1))] <- vote[[i]][j,k];
    }
    term[case.count] <- 1952+i;
    time[case.count] <- i;
  }
}

# unify the coding
big.data[big.data==9] <- NA;
big.data <- 1 - big.data;  # coding of data and model are inconsistent
big.data <- cbind(big.data,  term=term, time=time);

#d56 <- read.table('data/vote/1956.txt') # 10 Justices in 1956
#d61 <- read.table('data/vote/1961.txt') # 10 Justices in 1961
#d68 <- read.table('data/vote/1969.txt') # 8 Justices in 1969


# setting initial values and priors
tau2.start <- rep(.1, 29);
tau2.start[3] <- .001; #.001 for Douglas, .1 for anyone else
set <- c(1,3,5,6,10,14,21,24,27);
theta.start <- rep(0, 29);
theta.start[c(1,10)] <- 1;
theta.start[3] <- -3;
theta.start[5:6] <- -2;
theta.start[14] <- -1;
theta.start[21] <- 2.0;
theta.start[c(24,27)] <- 2.5;  # 1.0 Har, -3 for Doug, -2 for Mar, -2 for Bren, 1 for Frankf, -1 for Fort, 2 for Rehn, 2.5 for Scalia, 2.5 for Thom

e0 <- rep(0, 29);
E0 <- rep(1, 29);
E0[set] <- .1;  # .1 for the 9 Justices above, and 1 for everyone else

# posterior inference using MCMCdynamicIRT1d function in MCMCpack package
# THIS FUNCTION WILL RUN FOR ABOUT 8 HOURS, BE PATIENT.
out <- MCMCdynamicIRT1d(t(big.data[,1:29]), item.time.map=big.data$time,theta.start=theta.start,mcmc=50000,burnin=20000,thin=5,verbose=500,tau2.start=tau2.start,eo=0,Eo=1,a0=0,A0=1,b0=0,B0=1,c0=-1,d0=-1, store.item=FALSE)

# out is a mcmc object, all the work below is to manipulate the data for drawing
suprecourt <- summary(out);
stat <- suprecourt$statistics;
sum(case.term);

quant <- suprecourt$quantiles;
dim(quant)

term.justice <- apply(adj!=-9, 2, sum);
sum(term.justice)
justice.name <- colnames(big.data)[-c(30,31)];
justice <- rep(justice.name, term.justice)


starting <- c()
count <- rep(1952, J);
term <- rep(NA, sum(term.justice));
for(i in 1:sum(term.justice)){
 
      count[which(justice[i]==justice.name)] <- count[which(justice[i]==justice.name)] + 1;
  
  term[i] <- count[which(justice[i]==justice.name)];
}

years <- function(x){
  which(x!=-9)+1952
}
aa <- apply(adj, 2, years);
time <- NULL;
for(i in 1:length(aa))
  time <- c(time, aa[[i]]);

justice <- factor(justice);
                   
stat.justice.only <- data.frame(stat[1:sum(term.justice),], time=time, justice=justice) # stupid of mistaking data.matrix as data.frame

# throw away justices that are not in the figure of original paper
stat.justice.only <- stat.justice.only[stat.justice.only$justice!='Breyer'&stat.justice.only$justice!='Burton'&stat.justice.only$justice!='Fortas'&stat.justice.only$justice!='Goldberg'&stat.justice.only$justice!='Jackson'&stat.justice.only$justice!='Minton'&stat.justice.only$justice!='Reed'&stat.justice.only$justice!='Whittaker'&stat.justice.only$justice!='Frankfurter',]

# draw the graph using ggplot2 package
a <- qplot(time, Mean, data=stat.justice.only, geom=c('path'))
fig.1 <- a+facet_wrap(~justice,ncol=5)+geom_path(aes(x=time, y=Mean+2*SD), color="grey")+geom_path(aes(x=time, y=Mean-2*SD), color="grey")+ylab('Ideology')
