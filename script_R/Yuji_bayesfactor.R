

library(BayesFactor)


# fig.3a
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ fig3a.txt', stringsAsFactors = FALSE)
input_vector <- result_table$V1
bf = ttestBF(input_vector, rscale=1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.3b
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ fig3b.txt', stringsAsFactors = FALSE)
input_vector <- result_table$V1
bf = ttestBF(input_vector, rscale=1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.4a
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ fig4a.txt', stringsAsFactors = FALSE)
input_vector <- result_table$V1
bf = ttestBF(input_vector, rscale=1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.4b
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ fig4b.txt', stringsAsFactors = FALSE)
input_vector <- result_table$V1
bf = ttestBF(input_vector, rscale=1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.5b
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig5a.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='b'), 1]
vector_lr <- result_table[which(result_table$V2=='lr'), 1]
bf = ttestBF(vector_lr, vector_b, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.5d left
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig5d\ left.txt', stringsAsFactors = FALSE)
result_table <- result_table$V1
bf = ttestBF(result_table, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.5d right
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig5d\ right.txt', stringsAsFactors = FALSE)
result_table <- result_table$V1
bf = ttestBF(result_table, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.6b
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig6b.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='b'), 1]
vector_lr <- result_table[which(result_table$V2=='lr'), 1]
bf = ttestBF(vector_b, vector_lr, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.6c hpc
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig6c\ hpc.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='b'), 1]
vector_lr <- result_table[which(result_table$V2=='lr'), 1]
bf = ttestBF(vector_b, vector_lr, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.6c phc
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig6c\ phc.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='b'), 1]
vector_lr <- result_table[which(result_table$V2=='lr'), 1]
bf = ttestBF(vector_b, vector_lr, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.6c prc
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig6c\ prc.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='b'), 1]
vector_lr <- result_table[which(result_table$V2=='lr'), 1]
bf = ttestBF(vector_b, vector_lr, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.6c erc
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig6c\ erc.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='b'), 1]
vector_lr <- result_table[which(result_table$V2=='lr'), 1]
bf = ttestBF(vector_b, vector_lr, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.7a 
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig7a.txt', stringsAsFactors = FALSE)
vector_b <- result_table[which(result_table$V2=='t1'), 1]
vector_lr <- result_table[which(result_table$V2=='t2'), 1]
bf = ttestBF(vector_b, vector_lr, paired = TRUE, rscale = 1)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])

# fig.7b
result_table <- read.table('/Users/bo/Documents/data_yuji_lab/PPT/manuscript\ egocentric/egocentric\ 2021-10-20\ Neuroimage/rawdata\ Fig7b.txt', stringsAsFactors = TRUE)

for (ith_row in 1:length(result_table[,1])) {
  if(result_table[ith_row, 2]=='mtl'){
    result_table[ith_row, 2]<-1
  }
  if(result_table[ith_row, 2]=='frontal'){
    result_table[ith_row, 2]<-2
  }
  if(result_table[ith_row, 3]=='t1'){
    result_table[ith_row, 3]<-1
  }
  if(result_table[ith_row, 3]=='t2'){
    result_table[ith_row, 3]<-2
  }
}


bf <- anovaBF(V1 ~ V2 + V3, data = result_table, whichModels='bottom', progress=FALSE)
samples = posterior(bf, iterations = 10000)
summary(samples)
median(samples[,3])


classical <- aov(result_table$V1 ~ result_table$V2*result_table$V3, data=result_table)
summary(classical)




classical <- aov(result_table$V1 ~ result_table$V2*result_table$V3, data=result_table)
summary(classical)




