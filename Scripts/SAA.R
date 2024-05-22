#! /usr/bin/env Rscript

########################################################################
# Integrating Functional Annotation information by the genomic BLUP    #
# with Multiple random effects (FIAM)                                  #
# Copyright (C) 2024  Zhenshuang Tang, Lilin Yin and Xiaolei Liu       #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
########################################################################

cat(paste0("Start time: ", Sys.time(), "\n"))
options(stringsAsFactors=F)
library("optparse")

## Parameter setting
args_list <- list(
  make_option("--map", type = "character", default = NULL,
              help = "INPUT: the filename of map file", metavar = "character"),
  make_option("--map_format", type = "character", default = "bim",
              help = "INPUT: there two formats for map, one is the 'bim' like *.bim file in PLINK binary format, 
              another is the 'map3' (Chr SNP Pos), there no header for both types", metavar = "character"),
  make_option("--anno", type = "character", default = NULL,
              help = "INPUT: the filename of annotation files", metavar = "character"),
  make_option("--outPath", type="character", default=NULL,
              help="INPUT: the output path", metavar="character"),
  make_option("--output_prefix", type = "character", default = "IFAM",
              help = "INPUT: the prefix of output (default:IFAM)", 
              metavar = "character") 
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## check the options
if (!file.exists(opt$map)){
  cat(paste0("ERROR: ", opt$map, " does not exist! Please check!\n"))
  q()
}
anno_str <- unlist(strsplit(opt$anno, ","))
for (i in 1:length(anno_str)){
  if (!file.exists(anno_str[i])){
    cat(paste0("ERROR: ", anno_str[i], " does not exist! Please check!\n"))
    q()
  }
}
if (!file.exists(opt$outPath)){
  cat(paste0("ERROR: ", opt$outPath, " does not exist! Please check!\n"))
  q()
}

## map files
map <- read.delim(opt$map, head=FALSE)
if (opt$map_format == "bim"){
    map_SNP <- map[,2]
    cat(paste0("The map file has ", nrow(map), " SNPs!\n"))
} else if(opt$map_format == "map3"){
    map_SNP <- map[,2]
    cat(paste0("The map file has ", nrow(map), " SNPs!\n"))
} else{
  cat(paste0("ERROR: the format of map file doesn't meet the requirements! Please check!\n"))
}

## annotation files
cat(paste0("Analysis ", length(anno_str), " annotations: \n"))
for (i in 1:length(anno_str)){
  cat(paste0(anno_str[i], "\n"))
}

cat(paste0("------------ Annotation information ------------\n"))
anno_list <- list(Annotations=NULL, Number_annotaions=NULL, Number_overlap_anno_map=NULL, Anno_cover_percent=NULL)
SNP_list <- list()
for (i in 1:length(anno_str)){
  anno <- read.delim(anno_str[i], head=FALSE)
  anno_str_str <- unlist(strsplit(anno_str[i], "/"))
  anno_list$Annotations[i] <- gsub(".txt", "", anno_str_str[length(anno_str_str)])
  anno_list$Number_annotaions[i] <- nrow(anno)
  anno_list$Number_overlap_anno_map[i] <- length(intersect(anno[,1], map_SNP))
  anno_list$Anno_cover_percent[i] <- anno_list$Number_overlap_anno_map[i]/length(map_SNP)
  SNP_list[[i]] <- anno
}
anno_list$Anno_cover_percent <- paste(round(100*anno_list$Anno_cover_percent, 4), "%", sep="")
print(do.call(cbind, anno_list))
write.table(anno_list, paste0(opt$outPath, opt$output_prefix, ".annotaion.information.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat(paste0("------------ Overlap between annotations ------------\n"))
overlap_matrix <- matrix(NA, nrow = length(anno_str), ncol = length(anno_str), dimnames=list(anno_list$Annotations, anno_list$Annotations))
for (i in 1:length(anno_str)){
  for (j in 1:length(anno_str)){
    if (i == j){
      overlap_matrix[i,j] <- "--"
    }else if (i < j){
      overlap_matrix[i,j] <- length(intersect(unlist(SNP_list[[i]]), unlist(SNP_list[[j]])))
    }
  }
}
print(overlap_matrix)
overlap_matrix[is.na(overlap_matrix)] <- " "
write.table(overlap_matrix, paste0(opt$outPath, opt$output_prefix, ".overlaped.SNPs.between.anntations.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

cat(paste0("End time: ", Sys.time(), "\n"))
