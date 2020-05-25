library(SUtools)
result  <- process_mortran("contrast.m", "contree")
writeLines(result$fortran, "contrast.f")

