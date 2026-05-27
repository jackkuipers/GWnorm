
library(WWA)

# Fisher data
data <- iris[101:150, 1:4]

G_av <- WWA(data, 1e4)
G_av
