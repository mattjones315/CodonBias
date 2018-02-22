using DataFrames;
using FastaIO;

# file pointer for variant file
variant_fp = ARGS[1];

# file pointer for where to write the count matrix
out_fp = ARGS[2];

# read in variants
df = readtable(variant_fp, separator='\t',header=true);

# read in and change to appropriate working directory
wd = ARGS[3];
cd(wd)

chr1 = readfasta("human_ancestor_1.fa");
chr2 = readfasta("human_ancestor_2.fa");
chr3 = readfasta("human_ancestor_3.fa");
chr4 = readfasta("human_ancestor_4.fa");
chr5 = readfasta("human_ancestor_5.fa");
chr6 = readfasta("human_ancestor_6.fa");
chr7 = readfasta("human_ancestor_7.fa");
chr8 = readfasta("human_ancestor_8.fa");
chr9 = readfasta("human_ancestor_9.fa");
chr10 = readfasta("human_ancestor_10.fa");
chr11 = readfasta("human_ancestor_11.fa");
chr12 = readfasta("human_ancestor_12.fa");
chr13 = readfasta("human_ancestor_13.fa");
chr14 = readfasta("human_ancestor_14.fa");
chr15 = readfasta("human_ancestor_15.fa");
chr16 = readfasta("human_ancestor_16.fa");
chr17 = readfasta("human_ancestor_17.fa");
chr18 = readfasta("human_ancestor_18.fa");
chr19 = readfasta("human_ancestor_19.fa");
chr20 = readfasta("human_ancestor_20.fa");
chr21 = readfasta("human_ancestor_21.fa");
chr22 = readfasta("human_ancestor_22.fa");
chrX = readfasta("human_ancestor_X.fa");

const chrdict = Dict(1 => chr1, 2 => chr2, 3 => chr3, 4 => chr4,
                5 => chr5, 6 => chr6, 7 => chr7, 8 => chr8,
                9 => chr9, 10 => chr10, 11 => chr11, 12 => chr12,
                13 => chr13, 14 => chr14, 15 => chr15, 16 => chr16,
                17 => chr17, 18 => chr18, 19 => chr19, 20 => chr20,
                21 => chr21, 22 => chr22, 23 => chrX);

const COMP = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');

function generatePossibleContexts(K::Number, prefix::String, contexts=[])

    bp = ["A", "G", "T", "C"];
    N = length(bp);

    if K == 0
        push!(contexts, prefix)
        return(contexts)
    end

    for i in 1:N
        nprefix = string(prefix, bp[i])
        push!(contexts, generatePossibleContexts(K-1, nprefix))
    end

    return(contexts)
end

function unfold(A)
    V = []
    for x in A
        if length(x) == 1
            push!(V, x[1])
        else
            append!(V, unfold(x))
        end
    end
    return(V)
end

function get_context(chrnum, pos, K)
    offset = Int((K-1) / 2);
    chr = chrdict[chrnum];
    return uppercase(chr[1][2][pos-offset:pos+offset])
end

contexts = unfold(generatePossibleContexts(3, ""));
nt = DataFrame(NT=["A", "C", "G", "T"]);
count_matrix = repmat([0,0,0,0], 1, length(contexts));
count_matrix = DataFrame(count_matrix);
names!(count_matrix, [Symbol("$i") for i in contexts]);
count_matrix = hcat(nt, count_matrix);

# For each position, let's count mutations from the ancestral sequence
for i in 1:length(df[:,:POS])

    if length(df[i,:REF]) > 1 || length(df[i,:ALT]) > 1
        continue
    end
    if df[i,:CHROM] == "X"
        chrnum = 23;
    else
        chrnum = parse(Int, df[i,:CHROM]);
    end
    pos = parse(Int, df[i,:POS]);

    context = get_context(chrnum, pos, 3);
    count_matrix[count_matrix[:,:NT] .== df[:,:alt], Symbol(context)] += 1

end

writetable(out_fp, count_matrix, separator='\t')
