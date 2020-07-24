library(SUtools)
result  <- process_mortran("contrast.m", "conTree")
writeLines(result$fortran, "contrast.f")

