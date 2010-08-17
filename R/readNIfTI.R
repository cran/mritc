readNIfTI <- function (filename, level = 1) {
    fileparts <- strsplit(filename, "\\.")[[1]]
    nparts <- length(fileparts)
    ext <- tolower(fileparts[nparts])
    if (ext == "gz") {
        filename.nii.gz <- filename
        filename.nii <- paste(c(fileparts[-c(nparts-1 ,nparts)], 
            "nii"), collapse = ".")
        filename.hdr <- paste(c(fileparts[-c(nparts-1 ,nparts)], 
            "hdr"), collapse = ".")
        filename.img <- paste(c(fileparts[-c(nparts-1, nparts)], 
            "img"), collapse = ".")
    }
    else if (ext == "nii") {
        filename.nii <- filename
        filename.nii.gz <- paste(c(fileparts[-nparts], 
            "nii.gz"), collapse = ".")
        filename.hdr <- paste(c(fileparts[-nparts], 
            "hdr"), collapse = ".")
        filename.img <- paste(c(fileparts[-nparts], 
            "img"), collapse = ".")
    }
    else if (ext == "hdr") {
        filename.hdr <- filename
        filename.img <- paste(c(fileparts[-length(fileparts)], 
            "img"), collapse = ".")
    }
    else if (ext == "img") {
        filename.hdr <- paste(c(fileparts[-length(fileparts)], 
            "hdr"), collapse = ".")
        filename.img <- filename
    }
    else {
        filename.nii <- paste(filename, ".nii", sep = "")
        filename.nii.gz <- paste(filename, ".nii.gz", sep = "")
        filename.hdr <- paste(filename, ".hdr", sep = "")
        filename.img <- paste(filename, ".img", sep = "")
    }
    if ((ext != "hdr") && (ext != "img") && (!is.na(file.info(filename.nii)$size) || !is.na(file.info(filename.nii.gz)$size))) {
        if(!is.na(file.info(filename.nii.gz)$size))
            con <- gzfile(filename.nii.gz, "rb")
        else 
            con <- file(filename.nii, "rb")
        header <- fmri:::read.NIFTI.header(con)
        if (!(header$magic == "n+1") && !(header$magic == "ni1")) 
            warning("Hmmm! Dont see the magic NIFTI string! Try to proceed, but maybe some weird results will occur!")
        bytes <- header$voxoffset - 348
        header$extension <- readBin(con, "raw", bytes)
    }
    else {
        if (is.na(file.info(filename.hdr)$size) | (file.info(filename.hdr)$size < 
            348)) 
            stop("Hmmm! This does not seem to be a NIFTI header (hdr/img-pair)! Wrong size or does not exist!")
        con <- file(filename.hdr, "rb")
        header <- fmri:::read.NIFTI.header(con)
        header$extension <- NULL
        close(con)
        if (is.na(file.info(filename.img)$size)) 
            stop("Hmmm! This does not seem to be a NIFTI header (hdr/img-pair)! img-file not found!")
        con <- file(filename.img, "rb")
    }
    dx <- header$dimension[2]
    dy <- header$dimension[3]
    dz <- header$dimension[4]
    dt <- header$dimension[5]
    dd <- if (header$dimension[1] == 5) 
        header$dimension[6]
    else 1
    endian <- header$endian
    if (header$datatype == 1) {
        what <- "raw"
        signed <- TRUE
        size <- 1
    }
    else if (header$datatype == 2) {
        what <- "int"
        signed <- FALSE
        size <- if (header$bitpix) 
            header$bitpix/8/dd
        else 2
    }
    else if (header$datatype == 4) {
        what <- "int"
        signed <- TRUE
        size <- if (header$bitpix) 
            header$bitpix/8/dd
        else 2
    }
    else if (header$datatype == 8) {
        what <- "int"
        signed <- TRUE
        size <- if (header$bitpix) 
            header$bitpix/8/dd
        else 4
    }
    else if (header$datatype == 16) {
        what <- "double"
        signed <- TRUE
        size <- if (header$bitpix) 
            header$bitpix/8/dd
        else 4
    }
    else if (header$datatype == 32) {
        what <- "complex"
        signed <- TRUE
        size <- if (header$bitpix) 
            header$bitpix/8/dd
        else 8
    }
    else if (header$datatype == 64) {
        what <- "double"
        signed <- TRUE
        size <- if (header$bitpix) 
            header$bitpix/8/dd
        else 8
    }
    else {
        what <- "raw"
        signed <- TRUE
        size <- 1
    }
    ttt <- readBin(con, what, n = dx * dy * dz * dt * dd, size = size, 
        signed = signed, endian = endian)
    close(con)
    if (min(abs(header$pixdim[2:4])) != 0) {
        weights <- abs(header$pixdim[2:4]/min(abs(header$pixdim[2:4])))
    }
    else {
        weights <- NULL
    }
    dim(ttt) <- if (dd == 1) 
        c(dx, dy, dz, dt)
    else if (dt == 1) 
        c(dx, dy, dz, dd)
    else c(dx, dy, dz, dt, dd)
    if (dd == 1) {
        mask <- array(TRUE, c(dx, dy, dz))
        mask[ttt[, , , 1] < quantile(ttt[, , , 1], level, na.rm = TRUE)] <- FALSE
        dim(ttt) <- c(prod(dim(ttt)[1:3]), dim(ttt)[4])
        na <- ttt %*% rep(1, dim(ttt)[2])
        mask[is.na(na)] <- FALSE
        ttt[is.na(na), ] <- 0
        dim(mask) <- c(dx, dy, dz)
        z <- list(ttt = writeBin(as.numeric(ttt), raw(), 4), 
            format = "NIFTI", delta = header$pixdim[2:4], origin = NULL, 
            orient = NULL, dim = header$dimension[2:5], dim0 = header$dimension[2:5], 
            roixa = 1, roixe = dx, roiya = 1, roiye = dy, roiza = 1, 
            roize = dz, roit = 1:dt, weights = weights, header = header, 
            mask = mask)
        class(z) <- "fmridata"
    }
    else {
        z <- list(ttt = writeBin(as.numeric(ttt), raw(), 4), 
            format = "NIFTI", delta = header$pixdim[2:4], origin = NULL, 
            orient = NULL, dim = c(dx, dy, dz, dd), dim0 = c(dx, 
                dy, dz, dd), roixa = 1, roixe = dx, roiya = 1, 
            roiye = dy, roiza = 1, roize = dz, roit = 1:dd, weights = weights, 
            header = header)
    }
    invisible(z)
}
