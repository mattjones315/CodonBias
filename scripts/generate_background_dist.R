require(biomaRt)
require(Biostrings)
require("BSgenome.Hsapiens.UCSC.hg19")
setwd("~/Box Sync/ZaitlenLab/CodonBias/")

args = commandArgs(trailing=T)
fp = args[[1]]
output = args[[2]]

data = as.data.frame(data.table::fread(fp, skip=1))
data = data[,-1]
cn = as.matrix(read.table(fp, nrow=1))[1,]
colnames(data) = cn

# get hg19 / grch37 human genome build
ensembl=useMart("ensembl")
hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                      path="/biomart/martservice",dataset="hsapiens_gene_ensembl")

getflank <- function(position, alleles="[N/N]", chr="chr12", offset=10) {
  leftflank  <- getSeq(Hsapiens,chr,position-offset,position-1)
  rightflank <- getSeq(Hsapiens,chr,position+1,position+offset)
  snp = getSeq(Hsapiens, chr, position, position)
  paste(leftflank,snp,rightflank,sep="")
}

SNP <- function(chr="", pos="", ref="", alt="", alt_is_anc="", AF="", context="") {
  
  me = list(chr=chr, pos=pos, ref=ref, alt=alt, alt_is_anc=alt_is_anc, AF=AF, context=context)
  class(me) = append(class(me), "SNP")
  
  return(me)
  
}

generatePossibleContexts = function(K=3, prefix="", contexts=c()) {
  bp = c("A", "G", "T", "C")
  N = length(bp)
  if (K == 0) {
    contexts = c(contexts, prefix)
    return(contexts)
  }
  for (i in 1:N) {
    nprefix = paste0(prefix, bp[[i]])
    contexts = c(contexts, generatePossibleContexts(K=K-1, prefix=nprefix))
  }
  return(contexts)
}

contexts = generatePossibleContexts(K=3)
count_matrix = matrix(0L, nrow=4, ncol=length(contexts))
colnames(count_matrix) = contexts 
rownames(count_matrix) = c("A", "C", "G", "T")


for (i in 1:nrow(data)) {
  chr = data[i,"CHROM"]
  pos = data[i,"POS"]
  ref = data[i,'REF']
  alt = data[i,"ALT"]
  
  if (nchar(alt) > 1 || nchar(ref) > 1) {
    next 
  }
  if (is.na(alt) || is.na(ref)) {
    next
  }
  
  context = getflank(as.numeric(pos), chr=paste0("chr", chr), offset=1)
  der_context = context
  if (data[i,"AF"] < .05 || is.na(data[i,"AF"])) {
    # if allele frequency is less than 5%, then alt is considered to be the derived allele
    substr(der_context, 3, 3) = alt
    count_matrix[ref, der_context] = count_matrix[ref, der_context] + 1
  } else {
    substr(der_context, 3, 3) = ref
    count_matrix[alt, der_context] = count_matrix[alt, der_context] + 1
  }
}

write.table(count_matrix, file=output)