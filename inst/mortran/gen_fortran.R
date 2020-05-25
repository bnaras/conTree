library(SUtools)
result  <- process_mortran("contrast_v2.m", "contree")
writeLines(result$fortran, "contrast.f")

