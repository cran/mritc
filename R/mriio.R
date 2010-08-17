readMRI <- function(file, dim=NULL, format) {
    if(format == "rawb.gz"){
        if(! is.null(dim)){
            if(!is.vector(dim) || length(dim) != 3)
                stop("The dimensions of the file are wrong.")
        }
        f <- gzfile(file, open="rb")
        on.exit(close(f))
        data <- readBin(f,"integer", prod(dim), size = 1, signed = FALSE)
        array(data, dim)
    }
    else if(format == "analyze"){
        f.read.analyze.volume(file)
    }
    else if(format == "nifti"){
        data <- readNIfTI(file)
        array(readBin(data$ttt, what="numeric", n=prod(data$dim), size=4), dim=data$dim)
    }
    else
        stop("Cannot recognize the format of the file.")
}

writeMRI <- function(data, file, header=NULL, format) {
    if(!(is.array(data))){
        stop("Data has to be an array.")
    }
    else if(!(length(dim(data)) ==3 ||
              length(dim(data)) == 4 && dim(data)[4] == 1)){
        stop("Data has to be an array of dimension 3 or dimension 4 with the
             forth dimension equal to 1.")
    }
    if(format == "rawb.gz"){
        f <- gzfile(file, open="wb")
        writeBin(as.vector(data), f, size=1)
        on.exit(close(f))
    }
    else if(format == "analyze"){
        write.ANALYZE(data, header, file)
    }
    else if(format == "nifti"){
        writeNIfTI(data, header, file)
    }
    else
        stop("Cannot recognize the format of the file.")
}
