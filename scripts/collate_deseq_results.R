# PiGx RNAseq Pipeline.
#
# Copyright © 2019 Bora Uyar <bora.uyar@mdc-berlin.de>
#
# This file is part of the PiGx RNAseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.




# Collate DESeq results into one big matrix

args <- commandArgs(trailingOnly = TRUE)
mapper <- args[1]
inpDir <- args[2]
outDir <- args[3]

if (mapper=="hisat2"||mapper=="star") {
  mapper <- ""
}

strpat <- paste0(mapper, ".deseq_results.tsv$")
filenames <- dir(inpDir, pattern = strpat, full.names = TRUE)

#filenames <- list.files(pattern=".deseq_results.tsv$")

all_files <- lapply(filenames, function(x) {
  file <- read.table(x)
  sub_file <- subset(file, select = c(0,2,6))
  # Get the start of filename prefix 
  prefix = sub(".deseq_results.tsv", "", x)
  inpstr = paste(inpDir, "/", sep = "")
  prefix = sub(inpstr, "", prefix)
  colnames(sub_file) <- paste(prefix, colnames(sub_file), sep='_')
  return(sub_file)
})


multimerge <- function (mylist) {
  ## mimics a recursive merge or full outer join
  
  unames <- unique(unlist(lapply(mylist, rownames)))
  
  n <- length(unames)
  
  out <- lapply(mylist, function(df) {
    
    tmp <- matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] <- as.matrix(df)
    rm(df); gc()
    
    return(tmp)
  })
  
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
  
  bigout <- do.call(cbind, out)
  colnames(bigout) <- paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "")
  return(bigout)
}


collated_dataframe <- multimerge(all_files)


# save results to out file
if (mapper=="genes"||mapper=="transcripts") {
  mapper <- paste0(".", mapper)
}
finalname <- paste("collated", mapper, ".deseq_results.tsv", sep="")

collatedFile <- file.path(outDir, finalname)

write.table(x = as.data.frame(collated_dataframe), file = collatedFile, 
            quote = FALSE, sep = '\t')

#write.table(collated_dataframe, file ="collated.deseq_results.tsv", quote = FALSE, sep = '\t')
