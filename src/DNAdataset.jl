struct DNAdataset{S <: Real}
    #=    
    data can be 
        1. simulated data: Vector{sim_dna_str_w_motif}
        2. fasta that has labeled ground truth 
            (those contigugous substring with CAP latters):
            NamedTuple{(:str, :labeled_ground_truth), 
                Tuple{Vector{String}, Vector{String}}}
        3. fasta file without labeled ground tuth:
                Vector{String}
    =#
    raw_data::Union{Vector{sim_dna_str_w_motif},
                Vector{String},
                Vector{public_dna_str_w_motif},
                Nothing
                }    

    # motif: 
    # - if it's a simulated data, then it has a motif type
    # - if data is loaded from fasta as in case 2 or 3 above, 
    #   then it is a Nothing type
    motif::motif_type

    # number of data points
    N::Int

    # number of rows is the length of each transformed data_pt
    # number of columns is the total number of data entries
    # currently work for data of the same length only
    dat::Matrix{S}

    """
    Input: 
        motif: a motif type that generates the data
        num_data_pts: how many data points to generate
        GPU: (true/false) use GPU to store the wyk-data
            (so can use GPU for training)
    """
    function DNAdataset{S}(motif::motif_type, 
                        num_data_pts::Int
        ) where {S <: Real}
        raw_data = sample_backgound_with_motif_multiple(motif, num_data_pts);
            new(
                raw_data,
                motif,
                num_data_pts,
                sim_data_2_dummy(raw_data;F=S)
            )
    end

    """
    Input:
        filepath: fasta filepath 

        has_ground_truth: (true/false) the input fasta 
            file has each string contained capped-letter 
            DNA-alphabets that are labeled as motifs.

        same_len: (true/false) each string in the input 
            fasta file has the same length.

        GPU: (true/false) use GPU to store the wyk-data
            (so can use GPU for training)
    """    
    function DNAdataset{S}(dna_read::Vector{String}, 
                        motif::motif_type,
                        raw_data::Vector{public_dna_str_w_motif}
        ) where {S <: Real}       
        new(raw_data, 
            motif, 
            length(dna_read),
            DNA_strings_2_dummy_vectors_same_length(dna_read;F=S)
        );
    end 
    function DNAdataset{S}(dna_read::Vector{String}) where {S <: Real}
        new(nothing, 
        nothing, 
        length(dna_read),
        DNA_strings_2_dummy_vectors_same_length(dna_read; F=S)
        );
    end
end
    