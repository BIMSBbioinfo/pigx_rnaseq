# PiGx RNAseq Pipeline.
#
# Copyright © 2017, 2018 Jona Ronen <yablee@gmail.com>
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

# R script takes outputs from salmon and uses tximport package to create a matrix
# that is genes x samples

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
kallisto_output_folder <- args[1]
colDataFile <- args[2]

writeCounts <- function(colDataFile, kallisto_output_folder) {
  colData <- read.table(colDataFile)
  files <- file.path(kallisto_output_folder, rownames(colData), "abundance.h5")
  names(files) <- rownames(colData)
  txi <- tximport::tximport(files, type = "kallisto", txOut=TRUE)

  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData, ~group)

  write.table(x = DESeq2::counts(dds),
              file = file.path(kallisto_output_folder, paste0("counts_from_KALLISTO.tsv")),
              quote = F, sep = '\t')
}

writeCounts(colDataFile, kallisto_output_folder)
