
library(reprex)


m <- matrix(1:12, ncol = 4)

row_ind_num <- c(1,3)
# row_ind_chr <- c("1", "3")


m_new <- m[row_ind_num,]
# m_new <- m[row_ind_chr,]

fs::dir_tree("data")
