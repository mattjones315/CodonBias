using DataFrames;
using Distributions;
using CSV;
using Gadfly;

exac_counts_dir = ARGS[1];
cmap_fp = ARGS[2];

B = 2000;

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

context_codon_map = parse_context_file(cmap_fp);

exac_files = readdir(exac_counts_dir);

function generate_bootstrap_samples(exac_files)
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

function apply_colwise_chisq(o_counts, bg_counts)

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
    e_counts = zeros(Float64, 4, 64);

    for i in 1:nrow(e_freqs)
        for j in 2:ncol(e_freqs)
            e_counts[i, j-1] = e_freqs[i, j] * o_total_counts[j-1][1];
        end
    end

    # Keep track of our chisquared test statistics and pvalues
    chisq = zeros(Float64, 64, 1);
    pvalues = zeros(Float64, 64, 1);
    for j in 1:size(e_counts, 2)

        tstat = 0.0
        for i in 1:size(e_counts, 1)
            if e_counts[i, j] == 0
                continue
            end

            tstat += ((e_counts[i, j] - o_counts[i, j+1])^2) / e_counts[i, j][1]
        end
        chisq[j] = tstat
        pvalues[j] = 1.0 - cdf(Chisq(2), tstat)
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
pvalues = Array{Float64, 2}(B, 4);
tvalues = Array{Float64, 2}(B, 4);
for b in 1:B

    # give updates every 100 iterations of bootsrapping
    if b % 100 == 0
        println(string("Bootstrap iteration: ", b))
    end

    exac_counts, background_counts = generate_bootstrap_samples(exac_files);

    results = apply_contextwise_chisq(exac_counts, background_counts, context_codon_map)
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
