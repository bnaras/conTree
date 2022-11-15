#' Save the function f for calling from fortran
#' @param f the R function to be called using `.Fortran`
#' @export save_rfun
save_rfun  <- function(f) {
    .Call("store_rfun", f, PACKAGE = "conTree")
    invisible(TRUE)
}
