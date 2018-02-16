using DataFrames;
using Distributions;

exac_counts_fp = ARGS[1];
bg_fp = ARGS[2];
out_fp = ARGS[3];

exac_counts = readtable(exac_counts_fp, separator='\t', header=true);
background_counts = readtable(bg_fp, separator='\t', header=true);

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
        by the number of times you observed the mutation event in o_counts
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
            expected[i] = o_counts[i, j] * e_freqs[i, j];

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

results = apply_colwise_chisq(background_freqs, exac_counts);
pvalues_df = DataFrame(pvalues = results[2][:,1]);
tstat_df = DataFrame(tstat =results[1][:,1]);
context = DataFrame(context = names(background_freqs)[2:end]);
resultdf = hcat(context, pvalues_df, tstat_df);

writetable(out_fp, resultdf, separator='\t');
