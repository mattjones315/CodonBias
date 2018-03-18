
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

GeneCodonFreqTable <- function(name="",
                      frequencies=matrix,
                      counts = matrix(),
                      ct = "",
                      seqFile = "",
                      sequences = list()) {
  if (nrow(counts) == 0) {
    count.mat = matrix(0L, ncol=length(names(CODON_DICT)), nrow=length(sequences))
    colnames(count.mat) = names(CODON_DICT)
  } else {
    count.mat = counts
  }

  me = list(name = name,
          freq.mat=frequencies,
          counts = count.mat,
          ct = ct,
          seqfile = seqFile,
          sequences = sequences)
  class(me) <- append(class(me), "GeneCodonFreqTable")
  return(me)
}

GeneCodonFreqTable.readSequences <- function(x, indices=list()) {

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

  total.TranscriptGenes = unlist(lapply(total.names, function(x) {
    spl = strsplit(x, "|", fixed=T)
    return(unlist(spl)[[2]])
  }))

  total.TranscriptGenes = unlist(lapply(total.TranscriptGenes, function(x) {
    spl = strsplit(x, ".", fixed=T)
    return(unlist(spl)[[1]])
  }))

  total.geneNames = as.vector(kct[,"hg38.knownCanonical.protein"])
  total.geneNames = unlist(lapply(total.geneNames, function(x) {
    spl = strsplit(x, ".", fixed=T)
    return(unlist(spl)[[1]])
  }))

  inter = intersect(total.TranscriptNames, transcriptNames)
  m = match(inter, total.TranscriptNames)
  m2 = match(inter, transcriptNames);

  sequences = seqset[m]

  total.geneNames = total.geneNames[m2];
  x$count.mat = matrix(0L, nrow=length(total.geneNames), ncol=length(names(CODON_DICT)));
  rownames(x$count.mat) = total.geneNames
  colnames(x$count.mat) = names(CODON_DICT)

  x$sequences <- sequences

  return(x)
}

GeneCodonFreqTable.countCodons <- function(x, indices = list()) {
  
  message(paste0("Counting Codons for CodonFreqTable = ", x$name, "..."))
  
  if (length(indices) != 0) {
    toCount = x$sequences[indices]
  } else {
    toCount = x$sequences
  }
  
  for (i in 1:length(toCount)) {
    tf = trinucleotideFrequency(toCount[[i]], step=3)
    n = unlist(strsplit(names(toCount)[[i]], "|", fixed=T))[[2]]
    n = unlist(strsplit(n, ".", fixed=T))[[1]]
    
    x$count.mat[n, names(tf)] = x$count.mat[n, names(tf)] + tf
  }
  
  return(x)
}


GeneCodonFreqTable.normalizeByGroup <- function(x, norm.group = list()) {
  message(paste0("Normalizing Codon Counts by Group, CodonFreqTable = ", x$name, "..."))

  if (length(norm.group) == 0) {
    x$freq.mat = x$count.mat / sum(x$count.mat)
    return(x)
  }
  
  x$freq.mat = matrix(0L, nrow=nrow(x$count.mat), ncol=ncol(x$count.mat));
  rownames(x$freq.mat) = rownames(x$count.mat)
  colnames(x$freq.mat) = colnames(x$count.mat)
  for (g in norm.group) {
    
    for (n in 1:nrow(x$count.mat)) {

      x$freq.mat[n, g] = x$count.mat[n, g] / sum(x$count.mat[n, g])
    }
  }

  return(x)
}
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

amino_acids = unlist(unique(unname((CODON_DICT))))

# map amino acids -> list of codons
rev_codon_list = list()
for (aa in amino_acids) {

  laa = lapply(CODON_DICT, function(x) ifelse(x == aa, T, F))
  codons = names(which(laa == T))
  rev_codon_list[[aa]] = unlist(codons)

}

#Read in known Canonical Transcript table

cft = GeneCodonFreqTable(name="Canonical Transcripts", ct=canonicalTranscripts, seqFile=transcriptSequences)
cft = GeneCodonFreqTable.readSequences(cft)

cft = GeneCodonFreqTable.countCodons(cft)

cft = GeneCodonFreqTable.normalizeByGroup(cft, norm.group=rev_codon_list)

saveRDS(list(counts = cft$counts,
              freqs = cft$frequencies), saveFile)
