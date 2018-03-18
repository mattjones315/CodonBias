using DataFrames;
using Distributions;
using CSV;
using Gadfly;

exac_counts_fp = "context_gene_3nt/gly.context_gene.3nt.txt";
cmap_fp = "contexts/gly.context.txt";
gene_counts = "confounders/gc.gly.context_gene.3nt.txt";
#exac_counts_dir = ARGS[1];
#cmap_fp = ARGS[2];

B = 100;

function parse_context_file(cmap_fp)
    # Takes in a file that maps each codon to all possible contexts for this
    # count matrix. Mapping is helpful for running the Chi-squared test for
    # significance later.
    cmap = open(cmap_fp);
    lines = readlines(cmap);
    map_dict = Dict()

    for l in lines
        lsplit = split(l)
        codon = lsplit[1]
        map_dict[codon] = lsplit[2:end]
    end

    map_dict

end

function compute_confounder_similarity(gs, gc)

    simmat = zeros(Int, length(gs), length(gs))

    for i in 1:length(gs)
        for j in 1:length(gs)
            simmat[i, j] = abs(g_weights[gs[i]] - g_weights[gs[j]])
        end
    end

    simmat

end

context_codon_map = parse_context_file(cmap_fp);

all_data = readtable(exac_counts_fp, separator='\t', header=true);
gene_counts = readtable(gene_counts, separator='\t', header=true);

genes_ = unique(all_data[:Gene]);
genes = [];

for g in genes_
    if length(split(g, ";")) == 1
        push!(genes, g)
    end
end

# Convert gene_counts Dataframe to hash table for quicker lookup
g_weights = Dict(i => 0 for i in genes);
for g in 1:size(gene_counts, 1)
    g_weights[gene_counts[:Gene][g]] = gene_counts[:Count][g]
end

gene_sim = compute_confounder_similarity(genes, gene_counts)


function select_contexts(exac, bg, cmap)

    map_keys = collect(keys(cmap));
    syn = []
    for m in 1:length(map_keys)
        for c in cmap[map_keys[m]]
            push!(syn, c)
        end
    end

    exac_syn = exac[[Symbol("$k") for k in syn]]
    bg_syn = bg[[Symbol("$k") for k in syn]]

    exac_syn = hcat(DataFrame(NT = convert(Array, bg[:,:NT])),
                                exac_syn);
    bg_syn = hcat(DataFrame(NT = convert(Array, bg[:, :NT])),
                                bg_syn);

    # Hack to account for when degrees of freedom < 2
    for i in 1:nrow(exac_syn)
        for j in 2:ncol(exac_syn)
            if exac_syn[i, j] == 0
                bg_syn[i, j] = 0
            end
        end
    end

    exac_syn, bg_syn
end

function select_indices(N, K)
    """
    Select indices from N possible indices for two sets S1 and S2 such that both
    S1 and S2 are of length K, and their weight (as determined by the gene_sim defined
    above) are as equal as possible.
    """

    bg_ii = zeros(Int, K);
    o_ii = zeros(Int, K);

    pii = [i for i in 1:N]

    for i in 1:K
        b = rand(pii, 1)[1]
        filter!(e -> e != b, pii)
        o = findmin(gene_sim[b, pii])[2];
        o = pii[o];
        filter!(e -> e != o, pii)
        bg_ii[i], o_ii[i] = b, o
    end

    bg_ii, o_ii
end

function compare_set_weights(s1, s2)

    w1 = 0.0;
    w2 = 0.0;

    for g in s1
        w1 += g_weights[g]
    end

    for g in s2
        w2 += g_weights[g]
    end

    abs(w1 - w2)
end


function generate_bootstrap_samples(all_counts)
    """
    Let's sample half of the genes with replacement to be our "background" dist and
    the other half to be the "observed" counts.
    """

    contexts = unique(all_counts[:Context]);

    # bg_ii = rand(1:length(genes), fld(length(genes), 2))
    # o_ii = setdiff([i for i in 1:length(genes)], bg_ii)
    bg_ii = rand(1:length(genes), 5000);
    o_ii = setdiff([i for i in 1:length(genes)], bg_ii);
    o_ii = rand(o_ii, 5000);
    #bg_ii, o_ii = select_indices(length(genes), 5000);

    bg_genes = genes[bg_ii];
    o_genes = genes[o_ii];

    println(string("Confounder Diff: ", compare_set_weights(bg_genes, o_genes)));

    bg_counts = all_counts[findin(all_counts[:Gene], bg_genes), :];
    o_counts = all_counts[findin(all_counts[:Gene], o_genes), :];

    bg_df = DataFrame(zeros(4, length(contexts)))
    o_df = DataFrame(zeros(4, length(contexts)))

    names!(bg_df, [Symbol("$k") for k in contexts])
    names!(o_df, [Symbol("$k") for k in contexts])

    for c in contexts

        bg = colwise(sum, bg_counts[bg_counts[:Context] .== c, 3:end])
        o = colwise(sum, o_counts[o_counts[:Context] .== c, 3:end])

        for i in 1:4
            bg_df[i, Symbol(c)] = bg[i][1];
            o_df[i, Symbol(c)] = o[i][1];
        end
    end

    o_df = hcat(DataFrame(NT = ["A", "T", "C", "G"]), o_df);
    bg_df = hcat(DataFrame(NT = ["A", "T", "C", "G"]), bg_df);

    o_df, bg_df
end

function generate_bootstrap_samples_chr(exac_files)
    """
    Let's select half of the chromosomes to be our 'background' dist and the other
    half to be the 'observed' counts.
    """

    bg_ii = rand(1:length(exac_files), fld(length(exac_files), 2))
    o_ii = setdiff([i for i in 1:length(exac_files)], bg_ii)

    background_counts = CSV.read(joinpath(exac_counts_dir, exac_files[bg_ii[1]]), delim='\t', header=true, nullable=false);
    for b in 2:length(bg_ii)
        b2 = CSV.read(joinpath(exac_counts_dir, exac_files[bg_ii[b]]), delim='\t', header=true, nullable=false);

        for i in 1:nrow(background_counts)
            for j in 2:ncol(background_counts)
                background_counts[i, j] += b2[i,j]
            end
        end
    end

    observed_counts = CSV.read(joinpath(exac_counts_dir, exac_files[o_ii[1]]), delim='\t', header=true, nullable=false);
    for o in 2:length(o_ii)
        o2 = CSV.read(joinpath(exac_counts_dir, exac_files[o_ii[o]]), delim='\t', header=true, nullable=false);

        for i in 1:nrow(observed_counts)
            for j in 2:ncol(observed_counts)
                observed_counts[i, j] += o2[i,j]
            end
        end
    end


    observed_counts, background_counts
end

function apply_colwise_chisq(o_counts, bg_counts, df)

    #=
    For each column, we'll run a chi squared test:
        1. Multiply the freq (empirical prob) of the e_freqs column / nt pair
        by the number of times you observed mutations in o_counts
        2. Apply a chi squared test
    E_FREQS and O_COUNTS should be the same dimension - NT x N_CONTEXT
    =#

    e_freqs = zeros(Float64, size(bg_counts,1), size(bg_counts, 2)-1)
    for i in 1:nrow(bg_counts)
        for j in 2:ncol(bg_counts)
            e_freqs[i, j-1] = bg_counts[i, j] / sum(bg_counts[:, j])
        end
    end
    e_freqs = DataFrame(e_freqs)
    rename!(e_freqs, f => t for (f,t) = zip(names(e_freqs),
                                names(bg_counts)[2:end]));
    e_freqs = hcat(DataFrame(NT = convert(Array, bg_counts[:,:NT])),
                                e_freqs);

    o_total_counts = count_colwise(o_counts[:,2:end]);
    e_counts = zeros(Float64, nrow(e_freqs), ncol(e_freqs)-1);
    for i in 1:nrow(e_freqs)
        for j in 2:ncol(e_freqs)
            e_counts[i, j-1] = e_freqs[i, j] * o_total_counts[j-1][1];
        end
    end

    # Keep track of our chisquared test statistics and pvalues
    chisq = zeros(Float64, ncol(e_freqs)-1, 1);
    pvalues = zeros(Float64, ncol(e_freqs)-1, 1);
    for j in 1:size(e_counts, 2)

        tstat = 0.0
        for i in 1:size(e_counts, 1)
            if o_counts[i, j+1] == 0
                continue
            end

            tstat += ((e_counts[i, j] - o_counts[i, j+1])^2) / e_counts[i, j][1]
        end
        chisq[j] = tstat
        pvalues[j] = 1.0 - cdf(Chisq(df), tstat)
    end

    chisq, pvalues, e_counts, o_counts
end


function count_colwise(df)

    #=
    Return column-wise sums for the dataframe DF.

    Useful because the column wise operations can be a little confusing in Julia,
    and it might be nice to abstract out this function for changing later.
    =#

    counts = colwise(sum, df);
    return counts;

end

function count_rowwise(df)

    counts = df[:,2:end];
    n_counts = [sum(Array(counts)[i,:]) for i in 1:size(counts, 1)];

    n_counts

end



function apply_contextwise_chisq(o_counts, bg_counts, mapping)

    C = length(mapping)
    chisq = zeros(C)
    pvalues = zeros(C)
    observed = DataFrame(Float64, 0, C);
    expected = DataFrame(Float64, 0, C);

    names!(observed, [Symbol(n) for n in o_counts[:,:NT]])
    names!(expected, [Symbol(n) for n in o_counts[:,:NT]])

    codons = DataFrame(Codon = collect(keys(mapping)));

    map_keys = collect(keys(mapping));

    for m in 1:length(map_keys)
        cs = mapping[map_keys[m]];
        bg = bg_counts[:, [Symbol(k) for k in cs]]
        o = o_counts[:, [Symbol(k) for k in cs]];

        bg = [sum(Array(bg)[i,:]) for i in 1:size(bg, 1)];
        o = [sum(Array(o)[i,:]) for i in 1:size(o, 1)];

        e_freqs = bg / (sum(bg))
        e_counts = sum(o) * e_freqs

        push!(observed, o);

        push!(expected, e_counts);

        tstat = 0.0
        for i in 1:length(e_counts)
            if e_counts[i] == 0
                continue
            end
            tstat += (e_counts[i] - o[i])^2 / (e_counts[i])
        end
        chisq[m] = tstat
        pvalues[m] = 1.0 - cdf(Chisq(2), tstat)
    end

    observed = hcat(codons, observed);
    expected = hcat(codons, expected);

    expected_freqs = DataFrame(zeros(size(expected,1),
                                    size(expected,2) - 1));
    rename!(expected_freqs, f => t for (f,t) = zip(names(expected_freqs),
                                names(expected)[2:end]));
    expected_freqs = hcat(DataFrame(Codon = convert(Array, expected[:,:Codon])),
                                expected_freqs);

    # fill in transition freqs from background_counts
    for i in 1:size(expected, 1)
        for j in 2:size(expected, 2)
            expected_freqs[i, j] = expected[i, j] /
                                            sum(convert(Array, expected[i,2:end]));
        end
    end

    observed_freqs = DataFrame(zeros(size(expected,1),
                                    size(expected,2) - 1));
    rename!(observed_freqs, f => t for (f,t) = zip(names(observed_freqs),
                                names(expected)[2:end]));
    observed_freqs = hcat(DataFrame(Codon = convert(Array, expected[:,:Codon])),
                                observed_freqs);

    # fill in transition freqs from background_counts
    for i in 1:size(observed, 1)
        for j in 2:size(observed, 2)
            observed_freqs[i, j] = observed[i, j] /
                                            sum(convert(Array, observed[i,2:end]));
        end
    end

    chisq, pvalues, observed, expected, observed_freqs, expected_freqs
end

# Now let's bootstrap for B terationss
pvalues = Array{Float64, 2}(B, 16);
tvalues = Array{Float64, 2}(B, 16);
for b in 1:B

    # give updates every 100 iterations of bootsrapping
    if b % 1 == 0
        println(string("Bootstrap iteration: ", b))
    end

    exac_counts, background_counts = generate_bootstrap_samples(all_data);
    exac_counts, background_counts = select_contexts(exac_counts, background_counts, context_codon_map)

    results = apply_colwise_chisq(exac_counts, background_counts, 2)
    #results = apply_colwise_chisq(exac_counts, background_counts)
    pvalues[b,:] = results[2];
    tvalues[b,:] = results[1];


end

pvalues = DataFrame(pvalues);
names!(pvalues, [Symbol("$c") for c in collect(keys(context_codon_map))]);

tvalues = DataFrame(tvalues);
names!(tvalues, [Symbol("$c") for c in collect(keys(context_codon_map))]);


# Set plot = true if you'd like to plot the results
plot = false;
if plot
    plot(pvalues, x="GGT", Geom.histogram, Guide.title("P value Histogram, GGT"),
            Guide.xlabel("P values"), Guide.ylabel("Count"))

    plot(pvalues, x="GGA", Geom.histogram, Guide.title("P value Histogram, GGA"),
        Guide.xlabel("P values"), Guide.ylabel("Count"))

    plot(pvalues, x="GGC", Geom.histogram, Guide.title("P value Histogram, GGC"),
        Guide.xlabel("P values"), Guide.ylabel("Count"))

    plot(pvalues, x="GGG", Geom.histogram, Guide.title("P value Histogram, GGG"),
            Guide.xlabel("P values"), Guide.ylabel("Count"))

    all_pv = collect(Iterators.flatten(convert(Array, pvalues)))
    all_pv = DataFrame(PV = all_pv)
    plot(all_pv, x="PV", Geom.histogram, Guide.title("P value Histogram, All"),
            Guide.xlabel("P values"), Guide.ylabel("Count"))

end
