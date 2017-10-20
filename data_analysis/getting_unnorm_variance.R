#this script will take in unnormalized reads and create distributions for the distance away each pixel
#is given a specific signal.

library(plyr)
library(MASS)
#### Subroutines ####
print_out_data <- function(name, dataframe, o) {
  write.table(dataframe, file = gsub(" ","",paste(o, "/", name, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE, quote = FALSE)
}

#### MAIN ####
file_1 <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/orig_data/count_matrix_chr20_uA1.txt")
file_2 <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/orig_data/count_matrix_chr20_uA2.txt")

file_1 <- na.omit(file_1[,c("Bin1", "Bin2", "Signal")])
file_2 <- na.omit(file_2[,c("Bin1", "Bin2", "Signal")])
file_1$bin <- paste(file_1$Bin1, file_1$Bin2, sep = "-")
file_2$bin <- paste(file_2$Bin1, file_2$Bin2, sep = "-")

#get mean and variance
#comb <- cbind(file_1, file_2)
comb <- merge(file_1, file_2, by = "bin")

#add the original signal values to each pixel
comb$mean_sig <- rowMeans(comb[,c("Signal.x", "Signal.y")])

#order the signals
comb <- arrange(comb, mean_sig)
comb <- comb[,c("Bin1.x", "Bin2.x", "mean_sig", "Signal.x", "Signal.y")]
#check dists
prev=1
med_sig <- c()
var_sig <- c()
for (i in seq(1000,nrow(comb), by=450)) {
  med_sig <- append(med_sig,median(c(comb$Signal.x[prev:i], comb$Signal.y[prev:i])))
  var_sig <- append(var_sig, var(c(comb$Signal.x[prev:i], comb$Signal.y[prev:i])))
  prev = i
}
med_sig <- append(med_sig,mean(c(comb$Signal.x[prev:nrow(comb)], comb$Signal.y[prev:nrow(comb)])))
var_sig <- append(var_sig, var(c(comb$Signal.x[prev:nrow(comb)], comb$Signal.y[prev:nrow(comb)])))
hist(c(comb$Signal.x[(prev-900):(prev-1350)], comb$Signal.y[(prev-900):(prev-1350)]), breaks = 50)
sfxn <- splinefun(med_sig, var_sig)
plot(med_sig,var_sig)
plot(med_sig, sfxn(med_sig))
noise <- data.frame(med_sig, var_sig)
colnames(noise) <- c("Median", "Variance")
noise <- ddply(noise, .(Median), summarize, Variance=median(Variance))
plot(noise)
ggplot(noise, aes(Median, Variance))+
  geom_line()+
  ggtitle("Variance Between Replicates at Median Signal Values")+
  theme( plot.title = element_text(hjust = .5) )+
  labs( x="Median Signal of Pixel", y= "Variance Within Median Pixel Values")

outer <- "/Users/ncolaian/Documents/phanstiel_lab/data/"
print_out_data("noise", noise, outer)
