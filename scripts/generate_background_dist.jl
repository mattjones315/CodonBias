using DataFrames;
using FastaIO;
using JuliaDB;

println(ARGS[1])
chr = readfasta("chr1.fa.gz");
df = readtable("../chr1.intergenic.small.v2.txt", separator='\t')

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

function get_context(pos::Number, K::Number)
    offset = Int((K-1) / 2)
    return uppercase(chr[1][2][pos-offset:pos+offset])
end

contexts = unfold(generatePossibleContexts(3, ""));
nt_lookup = Dict("A"=>1, "C"=>2, "G"=>3, "T"=>4);
count_matrix = repmat([0,0,0,0], 1, length(contexts));
count_matrix = DataFrame(count_matrix);
names!(count_matrix, [Symbol("$i") for i in contexts]);

for i in 1:length(df[:,:POS])

    if length(df[i,:REF]) > 1 || length(df[i,:ALT]) > 1
        continue
    end

    alt_is_der = (df[i,:AF] < 0.05)
    context = get_context(df[i,:POS], 3)
    if isna(alt_is_der) || alt_is_der
        mut = string(context[1], df[i,:ALT], context[3])
        count_matrix[nt_lookup[df[i,:REF]], Symbol(mut)] += 1
    else
        mut = string(context[1], df[i,:REF], context[3])
        count_matrix[nt_lookup[df[i,:ALT]], Symbol(mut)] += 1
    end

end

count_matrix[:,:NT] = ["A", "C", "G", "T"]
writetable("chr1.count_matrix.txt", count_matrix, separator='\t')
