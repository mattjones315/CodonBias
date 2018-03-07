using DataFrames;
using Distributions;
using CSV;

exac_counts_fp = ARGS[1];
bg_fp = ARGS[2];
cmap_fp = ARGS[3];
out_fp = ARGS[4];

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

exac_counts = CSV.read(exac_counts_fp, delim='\t', header=true);
background_counts = CSV.read(bg_fp, delim='\t', header=true);

background_freqs = DataFrame(zeros(size(background_counts)));
rename!(background_freqs, f => t for (f,t) = zip(names(background_freqs),
                            names(background_counts)));
background_freqs[:,:NT] = background_counts[:,:NT];

# fill in transition freqs from background_counts
for i in 1:size(background_counts, 1)
    for j in 2:size(background_counts, 2)
        background_freqs[i, j] = background_counts[i, j] /
                                        sum(background_counts[:, j]);
    end
end

function apply_colwise_chisq(e_freqs, o_counts)

    #=
    For each column, we'll run a chi squared test:
        1. Multiply the freq (empirical prob) of the e_freqs column / nt pair
        by the number of times you observed mutations in o_counts
        2. Apply a chi squared test
    E_FREQS and O_COUNTS should be the same dimension - NT x N_CONTEXT
    =#

    # Keep track of our chisquared test statistics and pvalues
    chisq = zeros(Float64, 64, 1)
    pvalues = zeros(Float64, 64, 1)

    for j in 2:ncol(e_freqs)
        if sum(o_counts[:,j]) < 5
            chisq[j] = NaN;
            pvalues[j] = 1.0;
            continue;
        end
        observed = zeros(Float64, 4, 1);
        expected = zeros(Float64, 4, 1);
        for i in 1:nrow(e_freqs)
            observed[i] = convert(Float64, o_counts[i, j]);
            expected[i] = convert(Float64, o_counts[i, j]) * e_freqs[i, j];
        end

        # Calculate chisq test statistic
        tstat = 0.0
        for k in 1:length(observed)
            if expected[k] == 0
                continue
            end
            tstat += ((observed[k] - expected[k])^2 / expected[k]);
        end
        chisq[j-1] = tstat
    end

    # Calculate p value from chi squared test stat
    for c in 1:length(chisq)
        # if isnan, then we set the chisq test stat to NaN b/c there
        # were not enough counts in column. This is an approximation,
        # we could rather use Fisher's test to account for this

        if isnan(chisq[c])
            continue;
        end
        pvalues[c] = 1.0 - cdf(Chisq(2), chisq[c]);
    end

    chisq, pvalues
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

function count_rowise(df)

    counts = df[:,2:end];
    n_counts = [sum(Array(counts)[i,:]) for i in 1:size(counts, 1)];

    n_counts

end



function apply_contextwise_chisq(e_freqs, o_counts, e_counts, mapping)

    C = length(mapping)
    chisq = zeros(C)
    pvalues = zeros(C)
    bg = DataFrame(Float64, 0, C);
	observed = DataFrame(Float64, 0, C);
	expected = DataFrame(Float64, 0, C);

	names!(observed, [Symbol(n) for n in e_freqs[:,:NT]])
	names!(expected, [Symbol(n) for n in e_freqs[:,:NT]])

	codons = DataFrame(Codon = collect(keys(mapping)));

    map_keys = collect(keys(mapping));

    for m in 1:length(map_keys)
        cs = mapping[map_keys[m]];
        e = e_freqs[:, [Symbol(k) for k in cs]];
        o = o_counts[:, [Symbol(k) for k in cs]];
        b = o_counts[:, [Symbol(k) for k in cs]];
        o_total_counts = count_colwise(o);
        e_counts = zeros(Float64, 4, length(cs));

        for i in 1:length(cs)
            e_counts[:, i] = e[:, i] * o_total_counts[i];
        end
        #e_counts = DataFrame(e_counts);
        e_counts = [sum(e_counts[i,:]) for i in 1:size(e_counts, 1)];
        o = [sum(Array(o)[i,:]) for i in 1:size(o, 1)];
        b = [sum(Array(b)[i,:]) for i in 1:size(b, 1)];

		push!(observed, o);

		push!(expected, e_counts);

        push!(bg, b);

        tstat = 0.0
        for i in 1:length(e_counts)
            if e_counts[i] == 0
                continue
            end
            tstat += (e_counts[i] - o[i])^2 / (e_counts[i])
        end
        chisq[m] = tstat
        pvalues[m] = 1.0 - cdf(Chisq(2), chisq[m])
    end

	observed = hcat(codons, observed);
	expected = hcat(codons, expected);

    chisq, pvalues, observed, expected
end

results = apply_contextwise_chisq(background_freqs, exac_counts, background_counts, context_codon_map);
tcounts_background = DataFrame(background_counts = count_rowise(background_counts));
tcounts_exac = DataFrame(exac = count_rowise(results[3]));
pvalues_df = DataFrame(pvalues = results[2]);
tstat_df = DataFrame(tstat = results[1]);
codons = DataFrame(codon = collect(keys(context_codon_map)))
resultdf = hcat(codons, pvalues_df, tstat_df, tcounts_background, tcounts_exac);

writetable(out_fp, resultdf, separator='\t');
