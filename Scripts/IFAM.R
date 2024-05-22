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

## making GRMs for some sets of SNPs lists
## args
## set_snplists: sets of SNPs list for each annotation region
## bfile: the prefix of genotype file (plink binary format)
## weight: optional, the filename of SNP weight file (default: NULL)
## outPath: the output path
## thread: the number of threads (default: 1)
make_GRM <- function(set_snplists=NULL, bfile=NULL, weight=NULL, outPath=NULL, thread=1){
  set_snplists_str <- unlist(strsplit(set_snplists, ","))
  GRM_list <- c()
  for(i in 1:length(set_snplists_str)){    
    snplist <- set_snplists_str[i]
    GRM_name <- unlist(strsplit(snplist, "/"))
    GRM_name <- gsub(".txt", "", GRM_name[length(GRM_name)])
    GRM_name <- paste0(outPath, GRM_name) 
    if (!is.null(opt$weight)){
      makeGRM_cmd <- paste0("hiblup --bfile ", bfile, " --extract ", snplist, " --make-xrm --snp-weight ", weight, " --threads ", thread, " --out ", GRM_name, ".w")
      if (!file.exists(paste0(GRM_name, ".w.GA.bin"))){
        system(makeGRM_cmd)
      }
      GRM_list <- c(GRM_list, paste0(GRM_name, ".w.GA"))
    } else{
      makeGRM_cmd <- paste0("hiblup --bfile ", bfile, " --extract ", snplist, " --make-xrm --threads ", thread, " --out ", GRM_name)
      if (!file.exists(paste0(GRM_name, ".GA.bin"))){
        system(makeGRM_cmd)
      }
      GRM_list <- c(GRM_list, paste0(GRM_name, ".GA"))
    } 
  }
  GRM_list <- paste(GRM_list, collapse=",") 
  return(GRM_list)
}

## Optimizing the number of random effects based on the values of variance components
## args
## vc_list: the values of variance components for each annotations
## randomMax: the maximium number of random effects in the model (default:5)
opt_random <- function(vc_list=NULL, randomMax=5){
  var_max <- max(vc_list[,2])
  random_list <- list()
  for (j in 1:randomMax){
    index <- which(vc_list[,2] <= var_max/10^(j-1) & vc_list[,2] > var_max/10^j)
    if (j == randomMax) {
      index <- which(vc_list[,2] <= var_max/10^(j-1) & vc_list[,2] >= 0)
    }
    if (length(index)==0){
        random_list[[j]] <- NULL
    } else {
        random_list[[j]] <- vc_list[index, 1]
    }
  }
  random_list <- random_list[vapply(random_list, Negate(is.null), NA)]
  return(random_list)
}

## Parameter setting
args_list <- list(
  make_option("--bfile", type = "character", default = NULL,
              help = "INPUT: the prefix of genotype file (plink binary format)", metavar = "character"),
  make_option("--pheno", type = "character", default = NULL,
              help = "INPUT: the filename of phenotype file", metavar = "character"),
  make_option("--anno", type = "character", default = NULL,
              help = "INPUT: the filename of annotation files", metavar = "character"),
  make_option("--anno_spec", type = "character", default = NULL,
              help = "INPUT: the filename of special annotation file, this file doesn't participate 
              in the optimization of random effects in the model", metavar = "character"),
  make_option("--anno_GRM", type = "character", default = NULL,
              help = "INPUT: the filename of GRM for each annotation", metavar = "character"),             
  make_option("--weight", type = "character", default = NULL,
              help = "INPUT: the filename of SNP weight file", metavar = "character"),
  make_option("--VCfile", type = "character", default = NULL,
              help = "INPUT: the filename of the results of variance component estimation 
              for all anotations", metavar = "character"),
  make_option("--outPath", type="character", default=NULL,
              help="INPUT: the output path", metavar="character"),
  make_option("--output_prefix", type = "character", default = "IFAM",
              help = "INPUT: the prefix of output (default:IFAM)", 
              metavar = "character"), 
  make_option("--pheno_pos", type = "integer", default = "2",
              help = "INPUT: the position of the analyzed phenotype in columns of phenotype file (default:2)", 
              metavar = "character"),
  make_option("--randomMax", type = "integer", default = "5",
              help = "INPUT: the maximium number of random effects in the model (default:5)", 
              metavar = "character"),
  make_option("--VCmethod", type = "character", default = "AI",
              help = "INPUT: the method of variance component estimation (default: AI method)", 
              metavar = "character"),
  make_option("--thread", type = "integer", default = "1",
              help = "INPUT: the number of threads (default: 1)", 
              metavar = "character"),
  make_option("--tmp_files", type = "logical", default = TRUE,
              help = "INPUT: Whether temporary files are stored (default: TRUE)")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## check the options
bfile_str <- paste0(opt$bfile, c(".bed", ".bim", ".fam"))
if (!file.exists(bfile_str[1])){
  cat(paste0("ERROR: ", opt$bfile, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$pheno)){
  cat(paste0("ERROR: ", opt$pheno, " does not exist! Please check!\n"))
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
if (!is.null(opt$anno_spec)){
  anno_spec_str <- unlist(strsplit(opt$anno_spec, ","))
  if (length(anno_spec_str) > 1){
    for (i in 1:length(anno_spec_str)){
      if (!file.exists(anno_spec_str[i])){
        cat(paste0("ERROR: ", anno_spec_str[i], " does not exist! Please check!\n"))
        q()
      }
    }
  } else {
    if (!file.exists(opt$anno_spec)){
      cat(paste0("ERROR: ", opt$anno_spec, " does not exist! Please check!\n"))
      q()
    }
  }
}
if (!is.null(opt$anno_GRM)){
  anno_GRM_str <- unlist(strsplit(opt$anno_GRM, ","))
  if (is.null(opt$anno_spec)){
    if (length(anno_str) != length(anno_GRM_str)) {
      cat(paste0("ERROR: the number of GRMs and annotations must be the same! Please check!\n"))
      q()
    }
  } else {
    if ((length(anno_str) + length(anno_spec_str)) != length(anno_GRM_str)) {
      cat(paste0("ERROR: the number of GRMs and all annotations must be the same! Please check!\n"))
      q()
    }
  }
  for (i in 1:length(anno_GRM_str)){
    anno_GRM_str_str <- paste0(anno_GRM_str[i], c(".bin", ".id"))
    for(j in 1:length(anno_GRM_str_str)){
      if (!file.exists(anno_GRM_str_str[j])){
        cat(paste0("ERROR: ", anno_GRM_str_str[j], " does not exist! Please check!\n"))
        q()
      }
    }
  }
}
if (!is.null(opt$weight)){
  if (!file.exists(opt$weight)){
    cat(paste0("ERROR: ", opt$weight, " does not exist! Please check!\n"))
    q()
  }
}
if (!is.null(opt$VCfile)){
  if (!file.exists(opt$VCfile)){
    cat(paste0("ERROR: ", opt$VCfile, " does not exist! Please check!\n"))
    q()
  }
}

## phenotype file 
phe_header <- unlist(strsplit(readLines(opt$pheno, n=1), "\t"))
trait_name <- phe_header[opt$pheno_pos]
cat(paste0("Analysis Trait: ", trait_name, "\n"))

## annotation files
cat(paste0("Analysis ", length(anno_str), " annotations: \n"))
anno_str_names <- c()
for (i in 1:length(anno_str)){
  anno_str_str <- unlist(strsplit(anno_str[i], "/"))
  anno_str_names <- c(anno_str_names, gsub(".txt", "", anno_str_str[length(anno_str_str)]))
  cat(paste0(anno_str[i], "\n"))
}
anno_str2 <- matrix(anno_str, ncol=1, dimnames = list(anno_str_names, "Annotations"))

if (!is.null(opt$anno_spec)){
  anno_spec_names <- c()
  if (length(anno_spec_str) > 1){
    cat(paste0(length(anno_spec_str), " annotations don't participate in random effects optimization: \n"))
    for (i in 1:length(anno_spec_str)){
      anno_spec_str_str <- unlist(strsplit(anno_spec_str[i], "/"))
      anno_spec_names <- c(anno_spec_names, gsub(".txt", "", anno_spec_str_str[length(anno_spec_str_str)]))
      cat(paste0(anno_spec_str[i], "\n"))
    }
  } else {
    cat(paste0(length(length(anno_spec_str)), " annotation doesn't participate in random effects optimization: \n"))
    anno_spec_str_str <- unlist(strsplit(anno_spec_str, "/"))
    anno_spec_names <- c(anno_spec_names, gsub(".txt", "", anno_spec_str_str[length(anno_spec_str_str)]))
    cat(paste0(anno_spec_str, "\n"))
  }
  anno_spec_str2 <- matrix(anno_spec_str, ncol=1, dimnames = list(anno_spec_names, "Annotations"))
}

if (is.null(opt$anno_spec)){
  anno_all_str2 <- anno_str2
} else {
  anno_all_str2 <- rbind(anno_str2, anno_spec_str2)
}

## GRMs files
if (!is.null(opt$anno_GRM)){
  anno_GRM_names <- c()
  for (i in 1:length(anno_GRM_str)){
    anno_GRM_str_str <- unlist(strsplit(anno_GRM_str[i], "/"))
    anno_GRM_names <- c(anno_GRM_names, gsub(".GA", "", anno_GRM_str_str[length(anno_GRM_str_str)]))
  }
  anno_GRM_str2 <- matrix(anno_GRM_str, ncol=1, dimnames = list(anno_GRM_names, "GRMs"))
}

start <- proc.time()
## estimating variance components for optimizing the number of random effects 
if (!is.null(opt$VCfile)){
  vars <- read.delim(opt$VCfile, head=TRUE)
} else {
  ## making GRMs for each annotation for variance component estimation
  if (is.null(opt$anno_GRM)){
    if (is.null(opt$anno_spec)){
      GRMs <- make_GRM(set_snplists=opt$anno, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, thread=opt$thread)
    } else{
      GRMs_all <- make_GRM(set_snplists=c(opt$anno, opt$anno_spec), weight=NULL, bfile=opt$bfile, outPath=opt$outPath, thread=opt$thread)
      GRMs_all_str <- unlist(strsplit(GRMs_all, ","))
      GRMs_all_names <- c()
      for (i in 1:length(GRMs_all_str)){
        GRMs_all_str_str <- unlist(strsplit(GRMs_all_str[i], "/"))
        GRMs_all_names <- c(GRMs_all_names, gsub(".GA", "", GRMs_all_str_str[length(GRMs_all_str_str)]))
      }
      GRMs_all_str2 <- matrix(GRMs_all_str, ncol=1, dimnames = list(GRMs_all_names, "GRMs"))
      index <- match(row.names(anno_spec_str2), row.names(GRMs_all_str2))
      GRMs <- paste(GRMs_all_str[-index], collapse=",")   
    } 
  } else {
    if (is.null(opt$anno_spec)){
      GRMs <- opt$anno_GRM
    } else{
      index <- match(row.names(anno_spec_str2), row.names(anno_GRM_str2))
      GRMs <- paste(anno_GRM_str2[-index], collapse=",") 
    }
  }
  VC_cmd <- paste0("hiblup --single-trait --threads ", opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
  " --xrm ", GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc")
  system(VC_cmd)
  rm_cmd <- paste0("rm ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc.beta; rm ", 
                    opt$outPath, opt$output_prefix, "_", trait_name, "_vc.rand")
  system(rm_cmd)
  vars <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", trait_name, "_vc.vars"), head=TRUE)
  vars$Item <- gsub(".GA", "", vars$Item)
  vars <- vars[vars$Item!="e", 1:2]
}
if (length(anno_str) != nrow(vars)){
  cat(paste0("ERROR: the variance components doesn't match the number of annotations, Please check!\n"))
  q()
}

## Optimizing random effects in two scenarios
if (is.null(opt$anno_spec)){
  # No special annotations
  # Optimizing random effects for all annotations
  random_list <- opt_random(vc_list=vars, randomMax=opt$randomMax)
}else{
  # having special annotations, which doesn't participate in the optimization of random effects in the model
  # Optimizing random effects for other annotations
  random_list <- opt_random(vc_list=vars, randomMax=opt$randomMax - length(anno_spec_str))
  for(i in 1:length(anno_spec_str)){
    anno_spec_str_s <- anno_spec_str[i]
    anno_spec_str_s <- unlist(strsplit(anno_spec_str_s, "/"))
    anno_spec_str_s <- gsub(".txt", "", anno_spec_str_s[length(anno_spec_str_s)])
    random_list[[length(random_list) + i]] <- anno_spec_str_s
  }
}
end1 <- proc.time()
cat("Optimization time: ", end1[3]-start[3], "s \n")
cat(paste0("The number of random effects after optimization: ", length(random_list), "\n"))

## predicting additive genetic values using multiple random effects model

random_set_snplist_files <- list()
for (i in 1:length(random_list)){
  random_set <- random_list[[i]]
  random_set_snplist <- list()
  random_set_snplist_files[[i]] <- opt$outPath
  for(j in 1:length(random_set)){
    random_set_snplist_files[[i]] <- paste0(random_set_snplist_files[[i]], random_set[j])
    random_set_file <- anno_all_str2[rownames(anno_all_str2)==random_set[j]]
    random_set_snplist[[j]] <- read.delim(random_set_file, head=FALSE)
  }
  random_set_snplist <- unique(unlist(random_set_snplist))
  cat(paste0("Random effects set ", i, " including ", length(random_set_snplist), " SNPs!\n"))
  random_set_snplist_files[[i]] <- paste0(random_set_snplist_files[[i]], ".txt")
  if (!file.exists(random_set_snplist_files[[i]])){
    write.table(random_set_snplist, random_set_snplist_files[[i]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
}
random_set_snplist_files <- unlist(random_set_snplist_files)
if (is.null(opt$weight)){
  random_GRMs <- make_GRM(set_snplists=random_set_snplist_files, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, thread=opt$thread)
} else {
  random_GRMs <- make_GRM(set_snplists=random_set_snplist_files, bfile=opt$bfile, weight=opt$weight, outPath=opt$outPath, thread=opt$thread)
}
multipleBLUP_cmd <- paste0("hiblup --single-trait --threads ",  opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
" --xrm ", random_GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, opt$output_prefix, "_", trait_name)
system(multipleBLUP_cmd)
rand <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", trait_name, ".rand"), head=TRUE)
genetic_values <- apply(rand[2:(length(random_set_snplist_files)+1)], 1, sum)
output <- cbind(rand[, 1:(length(random_set_snplist_files)+1)], genetic_values, rand[,ncol(rand)])
colnames(output)[ncol(output)] <- "residuals"

## handling temporary files
if (!opt$tmp_files){
  rm_cmd <- paste0("rm ", opt$outPath, opt$output_prefix, "**")
  system(rm_cmd)
  rm_cmd <- paste0("rm ", opt$outPath, "**.bin")
  system(rm_cmd)
  rm_cmd <- paste0("rm ", opt$outPath, "**.id")
  system(rm_cmd)
  rm_cmd <- paste0("rm ", opt$outPath, "**.log")
  system(rm_cmd)
  rm_cmd <- paste0("rm ", opt$outPath, "**.txt")
  system(rm_cmd)
}

write.table(output, paste0(opt$outPath, opt$output_prefix, "_", trait_name, ".rand"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

end2 <- proc.time()
cat("Estimating genetic values time: ", end2[3]-end1[3], "s \n")
cat(paste0("End time: ", Sys.time(), "\n"))
