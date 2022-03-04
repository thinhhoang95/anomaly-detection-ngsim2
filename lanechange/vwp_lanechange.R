# Read the CSV file
dta <- read.csv('lc.csv')
# Convert to matrix
mat <- data.matrix(dta, rownames.force = NA)
# Remove all the column names and row names from the matrix
colnames(mat) <- NULL
rownames(mat) <- NULL
