using DataFrames;
using Distributions;
using CSV;

## Pass a file path that contains the counts of synonymous mutations found in ExAC
exac_counts_fp = ARGS[1];

## Pass a file path that contains the counts of "synonymous" mutations that will be treated as the null background
bg_fp = ARGS[2];

## Pass a context-codon mapping file, as seen in the contexts/ directory
cmap_fp = ARGS[3];

## The file path you wish to write out to.
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

exac_counts = CSV.read(exac_counts_fp, delim='\t', header=true, nullable=false);
background_counts = CSV.read(bg_fp, delim='\t', header=true, nullable=false);

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

function generate_possible_contexts(K, prefix, contexts=[])
    """
    This function will generate all possible K-length contexts. This is useful
    for populating the initial count matrix.
    """

    bp = ["A", "G", "T", "C"];
    N = length(bp);

    if K == 0
        push!(contexts, prefix);
        return(contexts);
    end

    for i in 1:N
        nprefix = string(prefix, bp[i]);
        push!(contexts, generate_possible_contexts(K-1, nprefix));
    end

    contexts
end

function unfold(A)
    """
    A basic utility function for "unfolding" a list of lists.
    """

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

function collapse_sequence(counts, k)
    """
    A utility function for collapsing some lmers to kmers count matrix. Very useful
    for when you have counts for 5-mers, for example, but would like counts for 3-mers.
    Takes the initial counts matrix COUNTS and your desired length of context K.
    """

    contexts = unfold(generate_possible_contexts(k, ""))
    nt = DataFrame(NT=["A", "C", "G", "T"]);
    count_matrix = repmat([0,0,0,0], 1, length(contexts));
    count_matrix = DataFrame(count_matrix);
    names!(count_matrix, [Symbol("$i") for i in contexts]);
    count_matrix = hcat(nt, count_matrix);

    l = length(string(names(counts)[2]));
    mid = convert(Int, ceil(l / 2));
    bal = convert(Int, fld(k, 2));

    for c1 in names(count_matrix)[2:end]
        for c2 in names(counts)[2:end]
            if string(c2)[mid-bal:mid+bal] == string(c1)
                for i in 1:nrow(counts)
                    count_matrix[i, Symbol(c1)] += counts[i, Symbol(c2)]
                end
            end
        end
    end

    count_matrix
end

function convert_to_R_dataframe(ocounts, ecounts)
    """
    Utility function for converting a julia dataframe, as populated in this script,
    to an R dataframe format that will be used for advanced plotting.
    """

    odf = DataFrame(context = string.(names(ocounts)[2:end]),
                    counts = convert(Array, ocounts[1, 2:end])[:],
                    Sub = fill(ocounts[1,:NT], ncol(ocounts) - 1))
    for i in 2:nrow(ocounts)
        df2 = DataFrame(context = string.(names(ocounts)[2:end]),
                        counts = convert(Array, ocounts[i, 2:end])[:],
                        Sub = fill(ocounts[i,:NT], ncol(ocounts) - 1))
        odf = vcat(odf, df2)
    end

    edf = DataFrame(context = string.(names(ocounts)[2:end]),
                    counts = convert(Array, ecounts[1, 2:end])[:],
                    Sub = fill(ocounts[1,:NT], ncol(ocounts) - 1))

    for i in 2:nrow(ocounts)
        df2 = DataFrame(context = string.(names(ocounts)[2:end]),
                        counts = convert(Array, ecounts[i, 2:end])[:],
                        Sub = fill(ocounts[i,:NT], ncol(ocounts) - 1))
        edf = vcat(edf, df2)
    end

    odf, edf

end

function apply_colwise_chisq(o_counts, bg_counts, df)

    """
    For each column, we'll run a chi squared test:
        1. Multiply the freq (empirical prob) of the e_freqs column / nt pair
        by the number of times you observed mutations in o_counts
        2. Apply a chi squared test with DF degrees of freedom
    BG_COUNTS and O_COUNTS should be the same dimension - NT x N_CONTEXT
    """

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

    """
    Return column-wise sums for the dataframe DF.

    Useful because the column wise operations can be a little confusing in Julia,
    and it might be nice to abstract out this function for changing later.
    """

    counts = colwise(sum, df);
    return counts;

end

function count_rowise(df)
    """
    Return row-wise sums for the dataframe DF.

    Useful because the row wise operations can be a little confusing in Julia,
    and it might be nice to abstract out this function for changing later.
    """

    counts = df[:,2:end];
    n_counts = [sum(Array(counts)[i,:]) for i in 1:size(counts, 1)];

    n_counts

end

function select_contexts(exac, bg, cmap)
    """
    A utility function for only choosing some contexts to be considered for the future
    statistical tests -- in particular, choosing synonymous contexts. It will choose
    contexts according to the context-codon mapping in CMAP in both the observed count matrix
    EXAC and background count matrix BG.
    """

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


function apply_codonwise_chisq(e_freqs, o_counts, e_counts, mapping)
    """
    We'd like to apply a codon-wise chisq test, which is comprised of collapsing
    all possible contexts for each codon and running the chi-squared test on the
    codon level rather than the context level.

    Takes E_FREQS as the expected frequencies, O_COUNTS as the observed counts, E_COUNTS
    as the total number of "synonymous" mutations in the expected group, and the contect-codon
    mapping MAPPING.
    """

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
            if o[i] == 0
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

function collapse_to_codon(data, mapping)
    """
    Utility funciton for collapsing to the codon level from contexts using MAPPING.
    Will aggregate counts in DATA accordingly.
    """

    C = length(mapping)
    codons = DataFrame(Codon = collect(keys(mapping)));
    ndf = DataFrame(Float64, 4, 0)

    map_keys = collect(keys(mapping));

    for m in 1:length(map_keys)

        cs = mapping[map_keys[m]];
        counts = zeros(Float64, 4, length(cs));

        for i in 1:length(cs)
            counts[:, i] = data[:, Symbol(cs[i])];
        end

        counts = [sum(counts[i,:]) for i in 1:size(counts, 1)];

		ndf = hcat(ndf, counts);

    end
    names!(ndf, [Symbol("$k") for k in map_keys])
    ndf = hcat(DataFrame(NT = data[:NT]), ndf)
    ndf

end

function compare_to_codons(results, mapping)
    """
    Utility function for taking the test-specific RESULTS and collapsing to codons from
    contexts according to MAPPING. Useful if you'd like to tie some biological meaning
    to the differences between group mutation rates at the codon level.
    """

    expected, observed = DataFrame(results[3]), results[4]
    names!(expected, names(observed)[2:end])
    expected = hcat(DataFrame(NT = observed[:NT]), expected)

    o_codons = collapse_to_codon(observed, mapping)
    e_codons = collapse_to_codon(expected, mapping)

    o_codons, e_codons

end

function run_3mer_analysis(exac_counts, bg_counts, k, mapping, df)
    """
    Very simple pipeline for taking 5-mer counts and running a context-wise chi squared
    test at the 3-mer level. Takes observed counts EXAC_COUNTS, background counts BG_COUNTS, which are
    both 5-mer level counts and K the desired length of the k-mer to be studied (presumed to be 3 here).

    Will collapse according to MAPPING and will run the chi-squared test with DF degress of freedom.

    Will return the "observed" and "expected" dataframe after having been converted to
    R-compatible dataframes for future plotting. 
    """

    c3nt = collapse_sequence(exac_counts,k);
    b3nt = collapse_sequence(bg_counts, k);

    c_syn, c_bg = select_contexts(c3nt, b3nt, mapping);

    results = apply_colwise_chisq(c_syn, c_bg, df);
    odf, edf = compare_to_codons(results, mapping);

    odf, edf = convert_to_R_dataframe(odf, edf);

    odf, edf

end
