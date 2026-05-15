library(LAVA)
library(foreach)
library(doParallel)
registerDoParallel(cores = 100)


input = process.input(input.info.file="Documents/OP_SAR/lava/input.info.txt",           # input info file
                      sample.overlap.file="Documents/OP_SAR/lava/sample.overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="Documents/OP_SAR/lava/g1000_eur",                    # reference genotype data prefix
                      phenos=c("op","sar"))       # subset of phenotypes listed in the input info file that we want to process

loci = read.loci("/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

result = foreach(i = 1:nrow(loci), .packages='foreach') %dopar%{
  locus = process.locus(loci[i,], input, drop.failed = F)
  
  result = list()
  if(!is.null(locus)){
    res_univ = cbind(loci=i, run.univ(locus))
    result[[1]] = res_univ
    
    if(sum(locus$failed)==0){
      res_bivar = cbind(loci=i, run.bivar(locus))
      result[[2]] = res_bivar
    }
  }
 
  result
}





