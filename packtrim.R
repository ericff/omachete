## Should handle unpacking a .tar or .tar.gz file and then
## running trim galore

## For use with Rscript
## Example syntax
## Rscript packtrim.R example.tar.gz
## or, to choose directory for the input files and the output files
##   other than the directory containing packtrim.R:
## Rscript packtrim.R example.tar.gz /output/dir/name

## Outputs files from trim galore, and also a .txt file
##  with the names of these two files
##  called "trimgalore.output.file.names.txt"

## ASSUMES THAT .FASTQ FILES DON'T START WITH ./._

## NOTE THAT this will stop if there are errors, so one may want to
## change next setting. But actually subsequent steps will probably
## fail, so one should probably keep this as is, and just do
## try() around calls to this script.
stop.if.errors = TRUE

## NOTE THAT this creates a temporary directory, so one may want
## to program this so as to delete that.

## Setting for trim galore if one wants to change them; could
## add this to args if one wanted it not to be hard-coded

error.rate <- 0.1
minimum.read.length <- 20
stringency<- 1
quality <- 20


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
} else {
    thiswd <- getwd()
}

setwd(thiswd)

current.warning.option <- getOption("warn")

if (stop.if.errors){
    options(warn=2)
}

if (length(grep(pattern = "(\\.tar|\\.tar\\.gz)$", x=file.name))==0){
    warning(paste0("Filename is ", file.name, " but it must end in .tar or in .tar.gz"), immediate. = TRUE)
}

## Get names of files in archive
tempfile.tarlist <- "tempfile.tarlist.txt"
list.files.command <- paste0("tar -tf ", file.name, " > ", tempfile.tarlist)
system(command=list.files.command)
unpacked.file.base.names.1 <- readLines(con=tempfile.tarlist)

## Exclude weird ._ files on mac, assumes that .fastq files don't
## start with ./._ 
unpacked.file.base.names <- unpacked.file.base.names.1[-grep(pattern = "^\\.\\/\\._", x=unpacked.file.base.names.1)]



if (length(unpacked.file.base.names)!=2){
    warning(paste0("ERROR: there are not exactly 2 files in the tar archive, but there should be."), immediate. = TRUE)
}

## Unpack archive; note: unpacks them to directory thiswd
unpack.command <- paste0("tar -xvzf ", file.name)
system(unpack.command)


#######################################################################
#######################################################################
## Make nice names for fastq files for input to trim galore,
## mainly to make sure the reverse and forward files are in
## the right order.
## Define function to do this work.
#######################################################################
#######################################################################

##  Pass this a vector of two filenames
## This will always succeed, i.e. without an error.
## gives an output though of error.in.fastq.file.names = TRUE if there is a problem
determine.reverse.and.forward.fastq.files.with.error.warning <- function(filenames){
    #   first do checks that all take form _1.fastq
    #   and _2.fastq or R1, R2,
    #   .fq, etc.
    # do this by looking for one ending in _R1.fq or _1.fastq.gz etc.
    #   changing it into a 2 and checking that the filenames are the
    #   same
    #
    #
    ## initialize:
    error.in.fastq.file.names <- FALSE
    file.one <- NA
    file.two <- NA
    # first look for which file has the pattern one would expect
    if (length(filenames)!=2){
        error.in.fastq.file.names <- TRUE
        print(paste("ERROR ERROR in determine.reverse.and.forward.fastq.files.with.error.warning: input vector should have length 2 but it does not.\nThe filenames given as input are:\n", paste(filenames, collapse = "\n"), sep=""), sep="")
    }
    # 
    pattern1 = "(_1|_R1)(?=(\\.fastq|\\.fq|\\.fastq.gz|\\.fq.gz)$)"
    outgrep = grep(pattern=pattern1, x =filenames, perl = TRUE)
    {
    if (length(outgrep)==0){
        error.in.fastq.file.names <- TRUE
        print(paste("ERROR in determine.reverse.and.forward.fastq.files.with.error.warning: finding 0 files that\nmatch the pattern for what the forward file name (file number 1) should look like\nThe list of all of the filenames is:\n", paste(filenames, collapse = "\n"), sep=""), sep="")  
    }
    else if (length(outgrep)==2){
        error.in.fastq.file.names <- TRUE
        print(paste("ERROR in determine.reverse.and.forward.fastq.files.with.error.warning: finding 2 files that\nmatch the pattern for what the forward file name (file number 1) should look like\nThe filenames that match are:\n", paste(filenames[outgrep], collapse = "\n"), sep=""), sep="")
    }
    else {
        file.one <- filenames[outgrep]
        # file.two must be the other one
        file.two <- filenames[(3-outgrep[1])]
        ##
        pattern1a = "_1(?=(\\.fastq|\\.fq|\\.fastq.gz|\\.fq.gz)$)"
        pattern1b = "_R1(?=(\\.fastq|\\.fq|\\.fastq.gz|\\.fq.gz)$)"
        filename.that.should.match.with.second.step1 <- gsub(pattern = pattern1a, replacement="_2", x = file.one, perl = TRUE)
        filename.that.should.match.with.second <- gsub(pattern = pattern1b, replacement="_R2", x = filename.that.should.match.with.second.step1, perl = TRUE)
        ## check that it matches with second filename; if not, give error
        ##    if so, then proceed
        if (!(filename.that.should.match.with.second== file.two)){
            error.in.fastq.file.names <- TRUE
            print(paste("ERROR in determine.reverse.and.forward.fastq.files.with.error.warning: after changing _1 or _R1 in first file to _2 or _R2, respectively,\n it does not match the other file.First file is\n", file.one, "\nModified first file is:\n", filename.that.should.match.with.second, "\nSecond file is:\n", file.two, sep=""))  
        }
    }
    } # end if/else
    #
    # 
    list(error.in.fastq.file.names=error.in.fastq.file.names, forward.file=file.one, reverse.file=file.two)
}


# remove _1 or _2 if not right before .fastq* or .fq*
# See https://github.com/lindaszabo/KNIFE
# filename ="asaew_1_ad_1fa_R1we_R1.fastq"
# filename ="asa+ew_1_ad_1fa_R1we_1.fq"
## make.new.file.names.so.nice.for.knife(filename = "700D6AAXX_1_1.fastq")
# make.new.file.names.so.nice.for.knife("unalignedasa+ew_1_ad_1fa_R1we_1.fq")
## From KNIFE readme:
## "The file names for a given sample must start with an identical string that identifies the sample and then have either _1, _R1, _2, or _R2 identifying them as read1 or read2. Any part of the file name after this will be ignored. Reads may be in gzipped files or plain text fastq files, but the file extension must be .fq, .fq.gz, .fastq, or .fastq.gz. Read file names must not contain the string 'unaligned' as this will interfere with logic around identifying de novo junctions."
## Mistake I made, that I noticed in sep 2016
## Conider this: "700D6AAXX_1_1.fastq"
## it would have replaced _1_ first to give 11.fastq
## and that gets rid of the underscore before 1, which is expected
## if have _1_1.fq, want it to keep _ before 1.fq
## so first replace 
make.new.file.names.so.nice.for.knife <- function(filename){
    # replace every instance of _1_ NOT followed by either .fastq or .fq
    ##    or by 1.fastq, 1.fq, 2.fastq, 2.fq, 
    #    with 1; similarly for _1, _R1, _2_, _2, _R2
    ## newname1 <- gsub(pattern="_1_(?!(\\.fastq|\\.fq).*)", replacement = "1", x=filename, perl=TRUE)
    ## changing sep 2016, so don't remove last underscore, e.g.
    ##   from _1_1.fastq
    newname1 <- gsub(pattern="_1_(?!((1|2)\\.fastq|(1|2)\\.fq).*)", replacement = "1", x=filename, perl=TRUE)
    newname2 <- gsub(pattern="_1(?!(\\.fastq|\\.fq).*)", replacement = "1", x=newname1, perl=TRUE)
    newname3 <- gsub(pattern="_R1(?!(\\.fastq|\\.fq).*)", replacement = "1", x=newname2, perl=TRUE)
    ## changing sep 2016, so don't remove last underscore
    newname4 <- gsub(pattern="_2_(?!((1|2)\\.fastq|(1|2)\\.fq).*)", replacement = "2", x=newname3, perl=TRUE)
    newname5 <- gsub(pattern="_2(?!(\\.fastq|\\.fq).*)", replacement = "2", x=newname4, perl=TRUE)
    newname6 <- gsub(pattern="_R2(?!(\\.fastq|\\.fq).*)", replacement = "2", x=newname5, perl=TRUE)
    # remove + signs, bad for SB
    newname7 <- gsub(pattern="\\+", replacement="", x=newname6)
    # check that unaligned is not in the name
    gsub(pattern = "unaligned", replacement = "unaalligned", x=newname7)
}


output.from.determine <- determine.reverse.and.forward.fastq.files.with.error.warning(filenames=unpacked.file.base.names)

if (output.from.determine$error.in.fastq.file.names){
    warning(paste0("ERROR of some type in determine.reverse.and.forward.fastq.files.with.error.warning for file names\n", paste0(unpacked.file.base.names, collapse = "\n"), "\n"), immediate. = TRUE)
}

## put forward file first, reverse second, as one would expect
ordered.fastq.base.names.1 <- c(output.from.determine$forward.file, output.from.determine$reverse.file) 


## Now rename in case the file names have _1 or _2 or _R1 or _R2
##   in the middle of the name
ordered.fastq.base.names <- vector("character", length =2)
for (tti in 1:2){
    ordered.fastq.base.names[tti] <- make.new.file.names.so.nice.for.knife(filename= ordered.fastq.base.names.1[tti])
}


ordered.fastq.full.names <- file.path(thiswd, ordered.fastq.base.names)

#######################################################################
#######################################################################
## Now run trim galore
## Puts output in a new directory, so as to minimize
## any chance of mixing up with files from other trim
## galore runs that may be running in parallel
## 
#######################################################################
#######################################################################

## Include a time stamp on end of directory to minimize chance of
## duplication
nicetime <- format(Sys.time(),"%H%M%S")
tempoutputdir <- file.path(thiswd,paste0(gsub(pattern ="(_|\\.|fq|fastq)", replacement="", x=ordered.fastq.base.names[1]),nicetime))

dir.create(tempoutputdir)

## Saves output files to output directory; this
## differs from command call on seven
## bridges

trim.galore.command <- paste0("/Users/awk/gl/trimgalore/trim_galore --gzip --fastqc --paired ", ordered.fastq.full.names[1], " ", ordered.fastq.full.names[2], " --quality ", quality, " --stringency ", stringency, " -e ", error.rate, " --length ", minimum.read.length, " --output_dir ", tempoutputdir)

## example
## trim_galore --gzip --fastqc --paired  140416_UNC11-SN627_0353_BC42B0ACXX_CGTACG_L005_1.fastq 140416_UNC11-SN627_0353_BC42B0ACXX_CGTACG_L005_2.fastq --quality 20 --stringency 1 -e 0.1 --length 20

system(trim.galore.command, wait= TRUE)

#######################################################################
#######################################################################
## Now check that file names are suitable names for use as inputs to
## knife, and change if necessary.
#######################################################################
#######################################################################

## Get file names produced by trim galore
## Assumes they have a certain format
## E.g.
## 1111129_UNC14-SN744_0193_BD09J3ACXX_ATCACG_L004_1_val_1.fq.gz
## 1111129_UNC14-SN744_0193_BD09J3ACXX_ATCACG_L004_2_val_2.fq.gz

## forward file first, reverse second, as one would expect
trim.output.reads <- c(dir(path=tempoutputdir, pattern= "_val_1\\.fq\\.gz"), dir(path=tempoutputdir, pattern= "_val_2\\.fq\\.gz"))


##  Rename val files, because of e.g. _1_val before _1.fq.gz
##    or _R1_val before _R1.fq.gz, which would give knife problems
## NOTE THAT the inputs must be ordered, 1 is forward, 2 reverse:
## Also moves the resulting files, whether the name is changed or
##  not, up from tempoutputdir to wdir
make.nice.names.of.validated.fastq.files <- function(trimmed.filenames.raw, wdir, tempoutputdir){
    error.in.fastq.file.names <- FALSE
    trimmed.paths.new <- c(NA,NA)
    ##
    ## Rename val files, because of e.g. _1_val before _1.fq.gz
    ##    or _R1_val before _R1.fq.gz
    trimmed.filenames.new <- vector("character", length=2) 
    trimmed.filenames.new[1] <- gsub(pattern="(_1_val|_R1_val)", replacement="val", x=trimmed.filenames.raw[1])
    trimmed.filenames.new[2] <- gsub(pattern="(_2_val|_R2_val)", replacement="val", x=trimmed.filenames.raw[2])
    ## Check that these new file names match
    ## This fails if they do not match; actually does more, but 
    ##  we do not care about the output. We just want
    ##   to know if there is a problem.
    out.determine <- determine.reverse.and.forward.fastq.files.with.error.warning(filenames = trimmed.filenames.new)
    {
        if (out.determine$error.in.fastq.file.names == TRUE) {
            error.in.fastq.file.names <- TRUE
        }
        else {
            {
                if (!identical(trimmed.filenames.new, trimmed.filenames.raw)){
                    ## Rename the files
                    print(paste0("Renaming files, specifically to get rid of strings like\n_1_val or _R1_val inserted into the file name by trim galore:\n Original names are:\n", paste(trimmed.filenames.raw, collapse="\n"), "\nNew names are:\n", paste(trimmed.filenames.new, collapse="\n")))

                }
                ## move the files to a new name possibly,
                ##  and always to wdir from temp directory
                for (ii in 1:2){
                    cmd.mv <- paste0("mv ", file.path(tempoutputdir,trimmed.filenames.raw[ii]),  " ", file.path(wdir, trimmed.filenames.new[ii]))
                    system(cmd.mv)
                }
            }
        } ## end else
    } ## end if/else
    list(error.in.fastq.file.names=error.in.fastq.file.names, trimmed.filenames.new=trimmed.filenames.new)
}


nice.names.output <- make.nice.names.of.validated.fastq.files(trimmed.filenames.raw=trim.output.reads, wdir=thiswd, tempoutputdir=tempoutputdir)

full.paths.of.trim.output <- file.path(thiswd, nice.names.output$trimmed.filenames.new)

writeLines(full.paths.of.trim.output, con =file.path(thiswd, "trimgalore.output.file.names.txt"))



## set warning option back to what it was before running this:
options(warn = current.warning.option)
