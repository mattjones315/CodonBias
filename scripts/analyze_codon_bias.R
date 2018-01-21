#setwd("~/Box Sync/ZaitlenLab/codon_bias/Codons/")

require(biomaRt)
require(Biostrings)

ARGS = commandArgs(trailingOnly = T)
canonicalTranscripts = ARGS[[1]]
transcriptSequences = ARGS[[2]]
if (length(ARGS) == 3) {
  saveFile = ARGS[[3]]
} else {
  saveFile = "cft_out.rds"
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

CodonFreqTable <- function(name="",
                      frequencies=list(),
                      counts = list(),
                      ct = "",
                      seqfile = "",
                      sequences = list()) {
  if (length(counts) == 0) {
    count.list = rep(0, 64)
    names(count.list) = names(CODON_DICT)
  } else {
    count.list = counts
  }

  me = list(name = name,
          frequences=frequencies,
          counts = count.list,
          ct = ct,
          seqfile = seqfile,
          sequences = sequences)
  class(me) <- append(class(me), "CodonFreqTable")
  return(me)
}

CodonFreqTable.readSequences <- function(x, indices=list()) {

  message(paste0("Reading Sequences for CodonFreqTable = ", x$name, "..."))

  if (length(indices) != 0) {
    print('only read in subset of sequences')
  }
  kct = read.table(x$ct, sep="\t", header=T)
  seqset = readDNAStringSet(x$seqfile)

  transcriptNames = as.vector(kct[,"hg38.knownToEnsembl.value"])
  total.names = names(seqset)

  total.TranscriptNames = unlist(lapply(total.names, function(x) {
    spl = strsplit(x, "|", fixed=T)
    return(unlist(spl)[[1]])
  }))

  inter = intersect(total.TranscriptNames, transcriptNames)
  m = match(inter, total.TranscriptNames)

  sequences = seqset[m]
  x$sequences <- sequences

  return(x)
}

CodonFreqTable.countCodons <- function(x, indices = list()) {

  message(paste0("Counting Codons for CodonFreqTable = ", x$name, "..."))

  if (length(indices) != 0) {
    toCount = x$sequences[indices]
  } else {
    toCount = x$sequences
  }

  for (i in 1:length(toCount)) { 
    tf = trinucleotideFrequency(toCount[[i]], step=3)
    x$counts = x$counts + tf
  }

  return(x)
}

CodonFreqTable.normalizeByGroup <- function(x, norm.group = list()) {
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

## Create an object that will function in a graph structure, connecting codons (nodes) to other codons with 
## edges that indicate that synonymous mutation events can take place between them. Edge weights will represent the
## frequency of these transitions. 
CodonNode <- function(name="", 
                      mutations=list(),
                      gene="", 
                      pos="", 
                      siginif=list(),
                      frequencies=list()) {
  ## codon is the name of the codon
  ## mutations is a list mapping mutation -> counts of mutation
  ## frequencies is a list mapping mutation -> frequency of mutation
  me = list(name=name, 
            mutations=mutations,
            frequencies=frequencies,
            gene=gene, 
            pos=pos, 
            signif=signif)
  class(me) <- append(class(me), "CodonNode")
  return(me)
}

CodonNode.incrementMutation <- function(x, der_codon) {
  x$mutations[[der_codon]]$counts = x$mutations[[der_codon]]$counts + 1
  return(x)
}

CodonNode.addMutation <- function(x, der_codon, syn) {
  muts = list()
  muts[[der_codon]]$counts = 1
  muts[[der_codon]]$syn = syn
  x$mutations = c(x$mutations, muts)
  return(x)
}

CodonNode.runHGT <- function(x) {
  #' Runs a Hypergeometric test (HGT) on the given transition probabilities to assess 
  #' statistical significance.
  
  return(1.0)
}

amino_acids = unlist(unique(unname((CODON_DICT))))

# map amino acids -> list of codons
rev_codon_list = list()
for (aa in amino_acids) {

  laa = lapply(CODON_DICT, function(x) ifelse(x == aa, T, F))
  codons = names(which(laa == T))
  rev_codon_list[[aa]] = unlist(codons)

}

#Read in known Canonical Transcript table

cft = CodonFreqTable(name="Canonical Transcripts", ct=canonicalTranscripts, seqfile=transcriptSequences)
cft = CodonFreqTable.readSequences(cft)

cft = CodonFreqTable.countCodons(cft)

cft = CodonFreqTable.normalizeByGroup(cft, norm.group=rev_codon_list)

saveRDS(list(counts = cft$counts,
              freqs = cft$frequencies), saveFile)




