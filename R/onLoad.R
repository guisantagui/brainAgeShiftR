.brainAgeShiftR_env <- new.env()

.onLoad <- function(libname, pkgname) {
        sysfile <- system.file("extdata", "sysdata.rda", package = pkgname)

        if (sysfile == "") {
                stop("sysdata.rda file not found in package.")
        }
        load(sysfile, envir = .brainAgeShiftR_env)
        if (interactive()) {
                message("brainAgeShiftR package loaded successfully.")
        }
}
