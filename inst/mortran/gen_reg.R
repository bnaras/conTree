## Script for generating registration (naras)
library(SUtools)
lines  <- readLines("../../src/contrast.f")
fs  <- unlist(fortran_subroutines(lines))
## Below is output of grep(".Fortran", *.R), edited for just names of called routines
fns  <- unique(
    c(
        "set_miss",
        "set_trm",
        "set_ntn",
        "set_qint",
        "set_pwr",
        "set_cri",
        "set_samp",
        "set_samp",
        "set_kri",
        "set_qqtrm",
        "set_quant",
        "set_vrb",
        "classin",
        "fcontrast",
        "get_stor",
        "set_miss",
        "prune1",
        "crinode",
        "getnodes1",
        "getlims",
        "set_samp",
        "set_samp",
        "set_kri",
        "set_qqtrm",
        "set_quant",
        "set_vrb",
        "classin",
        "andarm",
        "set_vrb",
        "fintcdf1",
        "set_vrb",
        "cdfpoints1",
        "trans",
        "untie")
)
l  <- sapply(fns, grep, x = fs, value = TRUE)
l$set_samp <- l$set_samp[[1]]
l$andarm  <- l$andarm[[1]]
reg_fns  <- unlist(l)
##registration  <- gen_registration(pkg_name="ConTree", fun_list = stringr::str_trim(reg_fns))
##writeLines(registration, con = "ConTree_init.c")
registration  <- gen_registration(pkg_name="contree", fun_list = stringr::str_trim(reg_fns))
writeLines(registration, con = "contree_init.c")
