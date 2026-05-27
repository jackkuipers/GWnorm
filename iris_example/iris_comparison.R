
library(Matrix)
library(BDgraph)
library(WWA)

# Fisher data
data <- iris[101:150, 1:4]

for (seednumber in 102:150) {#<- 101
print(seednumber)
top_i <- 7

resall <- NULL

# BD version bd - doesn't work at the moment
#method_val = "BD_bd"
#source("./iris_BD.R")
#resall <- rbind(resall, res)

# BD version rj
method_val = "BD_rj"
source("./iris_BD.R")
resall <- rbind(resall, res)

# BD MPL version bd
method_val = "BDmpl_bd"
source("./iris_BDmpl.R")
resall <- rbind(resall, res)

# BD MPL version rj
method_val = "BDmpl_rj"
source("./iris_BDmpl.R")
resall <- rbind(resall, res)

# WWA version
source("./iris_WWA.R")
resall <- rbind(resall, res)

print(resall)

save(resall, file = paste0("./iris_out/res", top_i, "seed", seednumber, ".Rdata"))

}