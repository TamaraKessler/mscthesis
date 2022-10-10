# Translation of Chris Sperber's code

# Clear environment
rm(list=ls())

library('R.matlab')
library('tidyr')
library('dplyr')
library('ggplot2')

# set random seed for reproducibility
set.seed(1)

# jzs function
jzs_corbf <- function(r,n){
  int <- function(r,n,g){
    (1+g)^((n-2)/2)*(1+(1-r^2)*g)^(-(n-1)/2)*g^(-3/2)*exp(-n/(2*g))
    }
  bf10 <- sqrt((n/2))/gamma(1/2)*integrate(int,lower=0,upper=Inf,r=r,n=n)$value
  return(bf10)	
}

test <- function(r,n){
  
  int2 <- function(r,n,g) {
    
    a <- ((n-2)/2) * log(1+g)
    b <- (-(n-1)/2) * log (1+(1-r^2)*g)
    c <- -3/2 * log(g)
    d <- -n/(2*g)
    
    return(ifelse(g == 0, 0, exp(a + b + c + d)))
    
  } 
  bf10 <- sqrt((n/2))/gamma(1/2)*integrate(int2,lower=0,upper=Inf,r=r,n=n)$value
  return(bf10)	
  
}



# define number of permutations in max stat permutation. Importantly, choose a large number 
#(at least 1.000, better 10.000 or even
# more), and make sure that the final p-threshold can be assessed with the
# resolution of p-values provided by the permutation (e.g. mapping results
# at p<0.0001 while having only 1000 permutations is non-sense, as the
# lowest possible p-values are 0.001 and 0)
# perms <- 50000

# read behav. data
behaviour <- read.csv('D:/Tamara/Diskonnektionen/Parcel_Disconnections/female/all_female_behavioural_neg.csv')
behaviour <- subset(behaviour[,2],!is.na(behaviour[,2]))
# read imaging data
#folder <- 'D:/Tamara/Diskonnektionen/Parcel_Disconnections/female/parcels'
folder <- 'D:/Tamara/Diskonnektionen/SSPLs/female/indirect/SSPLs'
fileList <- list.files(folder)


temp <- readMat(paste(folder,'/', fileList[1],sep=""))$idc.matrix
images_2d <- array(NA,c(nrow(temp),nrow(temp),length(fileList)))

for (i in 1:length(fileList)) {
  images_2d[,,i] <- readMat(paste(folder,'/',fileList[i],sep=""))$idc.matrix #idc.matrix
}

# create 2d mask for imaging data
mask <- array(0,c(nrow(temp),nrow(temp)))

# remove all elements in diagonal and blow
for (i in 1:nrow(temp)) {
  for (j in 1:nrow(temp)) {
    if (i<j) {
      mask[i,j] <- 1
    }
  }
}

# then we need to remove rarely/never affected connections

#  binarise the image, every disconnection >0 is set to 1
images_2d_binary <- replace(images_2d,images_2d>0,1)

image_2d_sum <- array(0,c(nrow(temp),nrow(temp)))
for (i in 1:nrow(temp)) {
  for (j in 1:nrow(temp)) {
    image_2d_sum[i,j] <- sum(images_2d_binary[i,j,])
  }
}
image_2d_sum <- image_2d_sum*mask

# set connection N threshold (e.g., minimum 15 patients) 
N <- 20
mask <- replace(mask,image_2d_sum<N,0)
mask_vect <- as.vector(mask)


# vectorise the 2d images
images_vect <- array(0,c(length(fileList),sum(mask_vect)))
for (i in 1:length(fileList)) {
  temp <- images_2d[,,i]
  temp <- as.vector(temp)
images_vect[i,] <- temp[mask_vect==1];
}

# do GLM with original data
stat_vect_orig <- array(0,sum(mask_vect))

# for (i in 1:length(stat_vect_orig)) { # runs GLM with each connection
#    dat <- c(behaviour,images_vect[,i])
#    my.glm <- glm(behaviour ~ images_vect[,i])
#    b <- my.glm[1]
#    dev <- summary(my.glm)[4]
#    # extracting t-values would be too complicated, I can see that it appears to have worked so far
#    stat_vect_orig[i] <- dev; # the minus sign changes the direction, see comment above on polarity
# }


# do bayesian analysis
stat_vect_bf <- stat_vect_orig
stat_vect_bf_orig <- stat_vect_bf



stat_vect_bf <- stat_vect_bf_orig

cors <- c(NA)
for (i in 1:length(stat_vect_orig)) {
  my.bayescor <- test(cor(behaviour,images_vect[,i]),length(behaviour))
  
  if(my.bayescor > 1000) {
    print(my.bayescor)
    print(images_vect[,i])
    print(ggplot(data = cbind(behaviour,as.data.frame(images_vect[,i])), 
                 aes(x=images_vect[,i], y = behaviour)) +
            geom_point() +
            geom_smooth(method="lm")
    )
    readline()
  }

  cors[i] <- cor(behaviour,images_vect[,i])
  #my.bayescor <- jzs_corbf(cor(behaviour,images_vect[,i]),length(behaviour)) #inputs are sample correlation r and number of observations n
  #cors[i] <- cor(behaviour,images_vect[,i])
  stat_vect_bf[i] <- my.bayescor
}

stat_vect_bf <- unlist(stat_vect_bf)

# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 6000000)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 50000)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 5000)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 1000)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 200)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 50)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 10)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf < 3)
# stat_vect_bf <- subset(stat_vect_bf,stat_vect_bf >3)
#stat_vect_bf <- replace(stat_vect_bf,stat_vect_bf<3,0)
# NOW I NEED TO GET THIS INTO A SHAPE WHERE I CAN VISUALISE IT
hist(stat_vect_bf)

# stat_vect_bf <- replace(stat_vect_bf,stat_vect_bf>(1/3),0)
# stat_vect_bf <- replace(stat_vect_bf,stat_vect_bf<(1/3)&stat_vect_bf>3,0)

# log results for them to look nicer
#stat_vect_bf <- log(stat_vect_bf)
#stat_vect_bf <- replace(stat_vect_bf,stat_vect_bf<3,0)

temp <- mask_vect
temp <- replace(temp,temp==1,stat_vect_bf)
results_2d <- matrix(temp,nrow(mask),nrow(mask))

#write.csv(results_2d,paste('D:/Tamara/Diskonnektionen/SSPLs/all/indirect/BF_results/',N,'_BF_nolog.csv',sep=''))
#write.csv(results_2d,paste('D:/Tamara/Diskonnektionen/Parcel_Disconnections/male/BF_results/',N,'_BF_nolog.csv',sep=''))
