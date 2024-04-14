library(singlet)

# Load the data
data <- read.csv("birds.csv", header = FALSE)
A <- as.matrix(data)
output <- run_nmf(A, rank = 10, verbose = 1, method = "scd")
