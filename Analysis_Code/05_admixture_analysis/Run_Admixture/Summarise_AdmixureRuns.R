#

cv_errors <- data.frame()
logLs <- data.frame()
k_vals <- vector()

for (K_VAL in 1:10){
  #Extract CV error from admixture log files
  system(paste0('grep "CV error" K',K_VAL,'/Run_*/projAdmix.*.log | cut -d " " -f 4 > cv_errors.temp.txt'))
  
  #Extract Loglikelihoods from admixture log files
  system(paste0('grep ^"Loglikelihood:" K',K_VAL,'/Run_*/projAdmix.*.log | cut -d " " -f 2 > Loglikelihoods.temp.txt'))
  
  #Load those files
  cv_errors <- rbind(cv_errors, read.delim("cv_errors.temp.txt", header = FALSE))
  logLs <- rbind(logLs, read.delim("Loglikelihoods.temp.txt", header = FALSE))
  
  #Create vector of K values
  k_vals <- c(k_vals, rep(K_VAL,100))
  
  #Remove temporary files
  system("rm cv_errors.temp.txt")
  system("rm Loglikelihoods.temp.txt")
}

#Add k_value column to data.frame
cv_errors$K <- k_vals
logLs$K <- k_vals

#Update column names
colnames(cv_errors) <- c("cv.error","K")
colnames(logLs) <- c("log.likelihood","K")

#Calculate mean CV error for each K value
library(dplyr)
cv_errors %>% 
  group_by(K) %>%
  summarise_at(vars(cv.error),list(mean.cv = mean)) -> mean.cv.errors

logLs %>%
  group_by(K) %>%
  summarise_at(vars(log.likelihood),list(mean.logLikelihood = mean)) -> mean.logLs

library(ggplot2)
P1 <- ggplot(data = mean.cv.errors, aes(x = K, y = mean.cv)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1,10, by = 1)) +
  labs(y = "Mean CV error")

P2 <- ggplot(data = mean.logLs, aes(x = K, y = mean.logLikelihood)) + 
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1,10, by = 1)) +
  labs(y = "Mean logLikelihood")

library(gridExtra)
pdf("CVerror_logLikelihood_byK.pdf", width = 8, height = 4)
grid.arrange(P1, P2, ncol = 2)
dev.off()


#############################################################

#library(remotes)
#remotes::install_github('royfrancis/pophelper')
library(pophelper)

#Create input for CLUMPP - this will output a new directory with the name "pop_KX"
#where X is the K value being summarised
Q.files <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE, pattern = ".Q")
Q.files <- Q.files[grepl("projAdmix", Q.files)]
clumppExport(readQ(Q.files), exportpath = getwd())

#Run CLUMPP for each K value, excl. K=1
PATH = getwd()
CLUMPP = paste0(PATH, "/CLUMPP_Linux64.1.1.2/CLUMPP")
for (K_VAL in 6:10){
  setwd(paste0(PATH, "/pop_K",K_VAL))
  system(paste0(CLUMPP, " paramfile"))
}
