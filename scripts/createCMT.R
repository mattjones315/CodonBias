#### Script for creating a Codon Mutational Table (CMT) from the ExAC dataset. 
require(biomaRt)
require(Biostrings)

ARGS = commandArgs(trailingOnly = T)
exac.fp = ARGS[[1]]

exac = data.table::fread(exac.fp, header=T)

## Only get mutations that are "Passing"
exac.passing = exac[which(exac[,"filter"] == 'PASS'),]

if (length(ARGS) == 3) {
  saveFile = ARGS[[3]]
} else {
  saveFile = "cmt_out.rds"
}

# map codon to amino acid
CODON_DICT = list("TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
                  "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
                  "TAT" = "Y", "TAC" = "Y", "TAA" = "STOP", "TAG" = "STOP",
                  "TGT" = "C", "TGC" = "C", "TGA" = "STOP", "TGG" = "W",
                  "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
                  "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
                  "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
                  "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
                  "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
                  "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
                  "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
                  "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
                  "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
                  "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
                  "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
                  "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G")

CodonMutTable <- function(name="",
                      frequencies=list(),
                      counts = list(),
                      data = matrix(0L, nrow=1, ncol=1)){

  if (length(counts) == 0) {
    count.mat = matrix(0L, nrow=64, ncol=64)
    rownames(count.mat) = names(CODON_DICT)
    colnames(count.mat) = names(CODON_DICT)
  } else {
    count.mat = counts
  }

  me = list(name = name,
          frequences=frequencies,
          counts = count.mat,
          data=data)
  class(me) <- append(class(me), "CodonFreqTable")
  return(me)
}

CodonMutTable.countMutations <- function(x, indices = list()) {

  message(paste0("Counting Mutations for CodonMutTable = ", x$name, "..."))

  if (length(indices) != 0) {
    toCount = x$data[indices,]
  } else {
    toCount = x$data
  }

  for (i in 1:nrow(toCount)) { 
      from = toupper(toCount[i,"codon_ref"])
      to = toupper(toCount[i,"codon_alt"])
      print(i)

      ## Skip NA's
      if (is.na(from) || is.na(to)) {
        next
      }

      ## Skip "N's" 
      if (grepl("N", from) || grepl("N", to)) {
        next 
      }

      ## Skip non-codons
      if (nchar(from) != 3 || nchar(to) != 3) {
        next
      }

      x$counts[from, to] = x$counts[from, to] + 1
    
  }

  return(x)
}

CodonMutTable.computeFreqs <- function(x) {
    
    x$freqs = apply(x$counts, 1, function(x) x / sum(x))
    return(x)
}

CodonMutTable.normalizeByGroup <- function(x, norm.group = list()) {
  message(paste0("Normalizing Codon Counts by Group, CodonFreqTable = ", x$name, "..."))

  if (length(norm.group) == 0) {
    x$frequences = x$counts / sum(x$counts)
    return(x)
  }

  for (g in norm.group) {
    x$frequencies[g] = x$counts[g] / sum(x$counts[g])
  }

  return(x)
}

amino_acids = unlist(unique(unname((CODON_DICT))))

# map amino acids -> list of codons
rev_codon_list = list()
for (aa in amino_acids) {

  laa = lapply(CODON_DICT, function(x) ifelse(x == aa, T, F))
  codons = names(which(laa == T))
  rev_codon_list[[aa]] = unlist(codons)

}


### Create Codon Mutation Table 
cmt = CodonMutTable(name="ExAC Codon Mutation Table", data=exac.passing)

cmt = CodonMutTable.countMutations(cmt)
cmt = CodonMutTable.computeFreqs(cmt)

saveRDS(list(counts = cmt$counts,
              freqs = cmt$frequencies), saveFile)

