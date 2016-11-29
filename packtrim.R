## Should handle unpacking a .tar or .tar.gz file and then
## running trim galore

## For use with Rscript

## Inputs: first is the file name for the tar or tar.gz file
##    second is optional and is the working directory to use
##    to write files

## this will write some temporary files to the current working directory
## and then remove them

args <- commandArgs(trailingOnly = TRUE)

## It's assumed that the file.name will be the full path,
##  or in the same directory as the working directory,
##  whether that be the same directory as the script or
##  the specified working directory
file.name = args[1]
## for testing:
## file.name = "testfile.tar.gz"
## file.name = "testfile.tar"

if (length(args)>1){
    thiswd <- args[2]
}
else {
    thiswd <- getwd()
}

setwd(thiswd)


if (length(grep(pattern = "(\\.tar|\\.tar\\.gz)$", x=file.name))==0){
    stop(paste0("Filename is ", file.name, " but it must end in .tar or in .tar.gz"))
}

## Get names of files in archive
tempfile.tarlist <- "tempfile.tarlist.txt"
list.files.command <- paste0("tar -tf ", file.name, " > ", tempfile.tarlist)
system(command=list.files.command)
unpacked.file.base.names <- readLines(con=tempfile.tarlist)
unpacked.file.full.names <- file.path(thiswd, unpacked.file.base.names)

## Unpack archive; note: unpacks them to directory thiswd
unpack.command <- paste0("tar -xvzf ", file.name)
system(unpack.command)


## Now check that file names are suitable names for use as inputs to
## knife.







                                 
