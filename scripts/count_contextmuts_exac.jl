using DataFrames;
using FastaIO;

exac_fp = ARGS[2];

# Read in all chromosome fasta files
chr1 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr1.fa.gz");
chr2 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr2.fa.gz");
chr3 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr3.fa.gz");
chr4 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr4.fa.gz");
chr5 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr5.fa.gz");
chr6 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr6.fa.gz");
chr7 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr7.fa.gz");
chr8 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr8.fa.gz");
chr9 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr9.fa.gz");
chr10 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr10.fa.gz");
chr11 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr11.fa.gz");
chr12 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr12.fa.gz");
chr13 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr13.fa.gz");
chr14 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr14.fa.gz");
chr15 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr15.fa.gz");
chr16 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr16.fa.gz");
chr17 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr17.fa.gz");
chr18 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr18.fa.gz");
chr19 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr19.fa.gz");
chr20 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr20.fa.gz");
chr21 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr21.fa.gz");
chr22 = readfasta("/ye/netapp/jimmie.ye/ref/genomes/hg19/chr22.fa.gz");

const chrdict = Dict(1 => chr1, 2 => chr2, 3 => chr3, 4 => chr4,
                5 => chr5, 6 => chr6, 7 => chr7, 8 => chr8,
                9 => chr9, 10 => chr10, 11 => chr11, 12 => chr12,
                13 => chr13, 14 => chr14, 15 => chr15, 16 => chr16,
                17 => chr17, 18 => chr18, 19 => chr19, 20 => chr20,
                21 => chr21, 22 => chr22);

const COMP = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');
exac_fp = "Codons/codon_table_v1.txt";
exac = readtable(exac_fp, separator='\t', header=true, nrows=3000);
exac_passing = exac[exac[:,:filter] .== "PASS", :];
exac_synonymous = exac_passing[exac_passing[:,:consequence] .== "synonymous_variant", :];

function generatePossibleContexts(K::Number, prefix::String, contexts=[])

    bp = ["A", "G", "T", "C"];
    N = length(bp);

    if K == 0
        push!(contexts, prefix);
        return(contexts);
    end

    for i in 1:N
        nprefix = string(prefix, bp[i]);
        push!(contexts, generatePossibleContexts(K-1, nprefix));
    end

    return(contexts)
end

function unfold(A)
    V = [];
    for x in A
        if length(x) == 1
            push!(V, x[1]);
        else
            append!(V, unfold(x));
        end
    end
    return(V);
end

function get_context(pos::Number,chrnum::Number,  K::Number)
    offset = Int((K-1) / 2);
    #chr = chrdict[chrnum];
    chr = chr1;
    return uppercase(chr[1][2][pos-offset:pos+offset])
end

function reverse_complement(seq)
    rseq = reverse(seq);
    join([COMP[c] for c in rseq])
end



contexts = unfold(generatePossibleContexts(3, ""));
nt = DataFrame(NT=["A", "C", "G", "T"]);
count_matrix = repmat([0,0,0,0], 1, length(contexts));
count_matrix = DataFrame(count_matrix);
names!(count_matrix, [Symbol("$i") for i in contexts]);
count_matrix = hcat(nt, count_matrix);

gly_exac = exac_synonymous[exac_synonymous[:,:aa_ref] .== "G",:];

for i in 1:length(gly_exac[:,:pos_id])
    pos_id = split(gly_exac[i,:pos_id], '_');
    chr = parse(Int, pos_id[1]);
    pos = parse(Int, pos_id[2]);
    if pos == 1 || isna(pos)
        continue;
    end

    context = get_context(pos, chr, 3)

    if isna(gly_exac[i, :alt_is_ancestral])

        if gly_exac[i,:freq_alt] < 0.05
            mut = string(context[1], gly_exac[i,:alt], context[3])
            count_matrix[count_matrix[:,:NT] .== gly_exac[i,:ref], Symbol(mut)] += 1
        else
            mut = string(context[1], gly_exac[i,:ref], context[3])
            count_matrix[count_matrix[:,:NT] .== gly_exac[i,:alt], Symbol(mut)] += 1
        end

    else
        if gly_exac[i, :alt_is_ancestral] == false
            mut = string(context[1], gly_exac[i,:alt], context[3])
            count_matrix[count_matrix[:,:NT] .== gly_exac[i,:ref], Symbol(mut)] += 1
        else
            mut = string(context[1], gly_exac[i,:ref], context[3])
            count_matrix[count_matrix[:,:NT] .== gly_exac[i,:alt], Symbol(mut)] += 1
        end
    end

end

writetable("exac.gly.context.txt", count_matrix, separator='\t')
