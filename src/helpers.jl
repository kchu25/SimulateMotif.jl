###### function that reads a fasta file ################
function reading(filepath::String, max_entries=max_num_read_fasta)
    # read the file
    reads = read(filepath, String);
    # process the fasta file to get the DNA part
    dna_reads = [join(split(i, "\n")[2:end]) for i in split(reads, '>') if !isempty(i)]; 
    dna_reads = length(dna_reads) > max_entries ? dna_reads[1:max_entries] : dna_reads;
    return dna_reads
end

function read_fasta(filepath::String, max_entries=max_num_read_fasta)
    #= read a fasta file =#
    dna_reads = reading(filepath, max_entries);   
    # convert all DNA seqeunce to uppercase
    return [uppercase(i) for i in dna_reads]
end

function read_fasta_that_has_ground_truth(filepath::String, max_entries=max_num_read_fasta)
    #= read the fasta file that contains the ground truth 
    (ground truths are in capital letters, as fasta files 
    from public databases typically contains it), return 
    1. dna reads and 2. ground truth 

    note: ok maybe it's not really ground truth since 
    by defn ground truth is unknown but we still keep that 
    labeled info
    =#
    
    #= read a fasta file =#
    dna_reads = reading(filepath, max_entries);   
    # ground truth motifs
    ground_truth_motif_strings = [join([j for j in i if isuppercase(j)]) for i in dna_reads];
    # convert all DNA seqeunce to uppercase
    dna_reads = [uppercase(i) for i in dna_reads];
    return (str = dna_reads, 
            labeled_ground_truth = ground_truth_motif_strings)
end

function get_JASPAR_PFM(jaspar_path::String)
    read_jaspar = read(jaspar_path, String);    
    s = split(read_jaspar, "\n");
    parsed = [[parse(Float64, i.match) for i in eachmatch(nmbr_regex, s[j])] for j = 2:length(s)-1];
    mat = reduce(hcat, parsed);
    return mat ./ sum(mat, dims=2)
end
########################################################

###### reverse complement of a string ##################
function reverse_complement(s::String)    
    join(islowercase(s[si]) ? s[si] : DNA_complement[s[si]] for si = length(s):-1:1)
end
########################################################

###### create data-set for as specified in DNAdataset###
function dna2dummy(dna_string::String; F=Float32)
    v = Array{F,1}(undef, 4*length(dna_string));
    for (index, alphabet) in enumerate(dna_string)
        start = (index-1)*4+1;
        v[start:start+3] = dummy[uppercase(alphabet)];
    end
    return v
end

function sim_data_2_dummy(dna_sim_data_vec::Vector{sim_dna_str_w_motif}; F=Float32)
    how_many_strings = length(dna_sim_data_vec);
    @assert F <: Real "input F must be subtype of Real"
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    len = length(dna_sim_data_vec[1].str); # length of each dna string in data    
    # note that for the simulated data here, each string always has 
    # the same length
    S = Array{F, 2}(undef, (4*len, how_many_strings));
    for i = 1:how_many_strings
        @inbounds S[:, i] = dna2dummy(dna_sim_data_vec[i].str);
    end
    return S
end

"""
get the set of dummy-vectors from a set of dna-strings
the dummy-vectors are all of same length (for now)
"""
function DNA_strings_2_dummy_vectors_same_length(dna_strings::Vector{String}; F=Float32)
    how_many_strings = length(dna_strings);
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    len = length(dna_strings[1]); # length of each dna string in data    
    S = Array{F, 2}(undef, (4*len, how_many_strings));
    for i = 1:how_many_strings
        @inbounds S[:, i] = dna2dummy(dna_strings[i]);
    end
    return S
end
########################################################