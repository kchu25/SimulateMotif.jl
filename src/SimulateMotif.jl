module SimulateMotif

######################### Dependencies #########################
using Distributions
################################################################

################# Exported methods and types ###################
export single_part_motif, 
       k_parts_motif,
       gapped_k_parts_motif, 
       mixture_k_parts_motifs, 
       mixture_gapped_k_parts_motifs, 
       read_fasta_that_has_ground_truth,
       DNAdataset
################################################################

###################### Load files ##############################
include("types.jl")
include("helpers.jl")
include("sample.jl")
include("DNAdataset.jl")


end
