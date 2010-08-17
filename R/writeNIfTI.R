writeNIfTI <- function (ttt, header = NULL, filename){
    if (is.null(header)) 
        header <- list()
    if (!("datatype1" %in% names(header))) 
        header$datatype1 <- paste(rep(" ", 10), collapse = "")
    if (!("dbname" %in% names(header))) 
        header$dbname <- paste(rep(" ", 18), collapse = "")
    if (!("extents" %in% names(header))) 
        header$extents <- c(0)
    if (!("sessionerror" %in% names(header))) 
        header$sessionerror <- c(0)
    if (!("regular" %in% names(header))) 
        header$regular <- "r"
    if (!("diminfo" %in% names(header))) 
        header$diminfo <- " "
    if (!("dimension" %in% names(header))) {
        header$dimension <- rep(0, 8)
        header$dimension <- c(length(dim(ttt)), dim(ttt))
    }
    if (length(header$dimension) < 8) 
        header$dimension[(length(header$dimension) + 1):8] <- 0
    if (!("intentp1" %in% names(header))) 
        header$intentp1 <- 0
    if (!("intentp2" %in% names(header))) 
        header$intentp2 <- 0
    if (!("intentp3" %in% names(header))) 
        header$intentp3 <- 0
    if (!("intentcode" %in% names(header))) 
        header$intentcode <- 0
    if (!("datatype" %in% names(header))) 
        header$datatype <- c(4)
    if (!("bitpix" %in% names(header))) 
        header$bitpix <- c(0)
    if (!("slicestart" %in% names(header))) 
        header$slicestart <- c(0)
    if (!("pixdim" %in% names(header))) 
        header$pixdim <- c(0, 4, 4, 4, rep(0, 4))
    if (length(header$pixdim) < 8) 
        header$pixdim[(length(header$pixdim) + 1):8] <- 0
    if (!("voxoffset" %in% names(header))) 
        header$voxoffset <- c(352)
    if (!("sclslope" %in% names(header))) 
        header$sclslope <- 0
    if (!("sclinter" %in% names(header))) 
        header$sclinter <- 0
    if (!("sliceend" %in% names(header))) 
        header$sliceend <- 0
    if (!("slicecode" %in% names(header))) 
        header$slicecode <- " "
    if (!("xyztunits" %in% names(header))) 
        header$xyztunits <- " "
    if (!("calmax" %in% names(header))) 
        header$calmax <- c(0)
    if (!("calmin" %in% names(header))) 
        header$calmin <- c(0)
    if (!("sliceduration" %in% names(header))) 
        header$sliceduration <- c(0)
    if (!("toffset" %in% names(header))) 
        header$toffset <- c(0)
    if (!("glmax" %in% names(header))) 
        header$glmax <- c(0)
    if (!("glmin" %in% names(header))) 
        header$glmin <- c(0)
    if (!("describ" %in% names(header))) 
        header$describ <- paste(rep(" ", 80), collapse = "")
    if (!("auxfile" %in% names(header))) 
        header$auxfile <- paste(rep(" ", 24), collapse = "")
    if (!("qform" %in% names(header))) 
        header$qform <- 0
    if (!("sform" %in% names(header))) 
        header$sform <- 0
    if (!("quaternb" %in% names(header))) 
        header$quaternb <- 0
    if (!("quaternc" %in% names(header))) 
        header$quaternc <- 0
    if (!("quaternd" %in% names(header))) 
        header$quaternd <- 0
    if (!("qoffsetx" %in% names(header))) 
        header$qoffsetx <- 0
    if (!("qoffsety" %in% names(header))) 
        header$qoffsety <- 0
    if (!("qoffsetz" %in% names(header))) 
        header$qoffsetz <- 0
    if (!("srowx" %in% names(header))) 
        header$srowx <- c(0, 0, 0, 0)
    if (!("srowy" %in% names(header))) 
        header$srowy <- c(0, 0, 0, 0)
    if (!("srowz" %in% names(header))) 
        header$srowz <- c(0, 0, 0, 0)
    if (!("intentname" %in% names(header))) 
        header$intentname <- paste(rep(" ", 16), collapse = "")
    if (!("magic" %in% names(header))) 
        header$magic <- "n+1"
    if (!("extension" %in% names(header))) 
        header$extension <- as.raw(rep(0, header$voxoffset - 
            348))
    dd <- if (header$dimension[1] == 5) 
        header$dimension[6]
    else 1
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
    con <- gzfile(paste(filename, ".nii.gz", sep = ""), "wb")
    writeBin(as.integer(348), con, 4)
    writeChar(header$datatype1, con, 10, NULL)
    writeChar(header$dbname, con, 18, NULL)
    writeBin(as.integer(header$extents), con, 4)
    writeBin(as.integer(header$sessionerror), con, 2)
    writeChar(header$regular, con, 1, NULL)
    writeChar(header$diminfo, con, 1, NULL)
    writeBin(as.integer(header$dimension), con, 2)
    writeBin(header$intentp1, con, 4)
    writeBin(header$intentp2, con, 4)
    writeBin(header$intentp3, con, 4)
    writeBin(as.integer(header$intentcode), con, 2)
    writeBin(as.integer(header$datatype), con, 2)
    writeBin(as.integer(header$bitpix), con, 2)
    writeBin(as.integer(header$slicestart), con, 2)
    writeBin(header$pixdim, con, 4)
    writeBin(header$voxoffset, con, 4)
    writeBin(header$sclslope, con, 4)
    writeBin(header$sclinter, con, 4)
    writeBin(as.integer(header$sliceend), con, 2)
    writeChar(header$slicecode, con, 1, NULL)
    writeChar(header$xyztunits, con, 1, NULL)
    writeBin(header$calmax, con, 4)
    writeBin(header$calmin, con, 4)
    writeBin(header$sliceduration, con, 4)
    writeBin(header$toffset, con, 4)
    writeBin(as.integer(header$glmax), con, 4)
    writeBin(as.integer(header$glmin), con, 4)
    writeChar(header$describ, con, 80, NULL)
    writeChar(header$auxfile, con, 24, NULL)
    writeBin(as.integer(header$qform), con, 2)
    writeBin(as.integer(header$sform), con, 2)
    writeBin(header$quaternb, con, 4)
    writeBin(header$quaternc, con, 4)
    writeBin(header$quaternd, con, 4)
    writeBin(header$qoffsetx, con, 4)
    writeBin(header$qoffsety, con, 4)
    writeBin(header$qoffsetz, con, 4)
    writeBin(header$srowx, con, 4)
    writeBin(header$srowy, con, 4)
    writeBin(header$srowz, con, 4)
    writeChar(header$intentname, con, 16, NULL)
    writeChar(header$magic, con, 4, NULL)
    bytes <- header$voxoffset - 348
    if (bytes != length(header$extension)) 
        warning("header$extension is not of expected size ", 
            bytes, " as given by header$voxoffset. cutting!")
    writeBin(header$extension[1:bytes], con)
    dim(ttt) <- NULL
    writeBin(ttt, con, size)
    close(con)
}
