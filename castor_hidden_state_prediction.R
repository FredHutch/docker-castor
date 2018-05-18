#!/usr/bin/env Rscript

# Load libraries
library(assertthat)
library(argparse)
library(parallel)
library(castor)
library(ape)

parser <- ArgumentParser(description='Perform hidden state prediction')

# Arguments passed in by the user
parser$add_argument('--newick', help='Newick tree')
parser$add_argument('--annotations', help='Tab-delimited table with known annotations')
parser$add_argument('--output', help='Path for output (tab-delimited table)')
parser$add_argument('--threads', help='Number of processors to use', default=1, type="integer")
parser$add_argument('--method', help='Method to use for HSP', default="pic", type="character")

args <- parser$parse_args()

# Read in the tree
print(paste("Reading in newick tree from", args$newick))
assert_that(is.readable(args$newick))
full_tree <- read.tree(args$newick)

# Read in the known annotations for a subset of the tips of the tree
print(paste("Reading in leaf traits from", args$annotations))
assert_that(is.readable(args$annotations))
trait_values <- read.delim(args$annotations, check.names=FALSE, row.names=1)

# Function to get HSP state probabilities for study (i.e. "unknown" tips only).
# Adds rownames of sequences and colnames of counts. 
# Also remove columns that are all zeros (no probability of that state).
get_sorted_prob <- function(in_likelihood, study_tips_i, tree_tips) {
  
  # Subet to study sequences only and set as rownames.
  tmp_lik <- in_likelihood[study_tips_i, , drop=FALSE]
  rownames(tmp_lik) <- tree_tips[study_tips_i]
  
  # Set column names to be 0 to max num of counts.
  colnames(tmp_lik) <- c(0:(ncol(tmp_lik)-1))
  
  # Remove columns that are 0 across all sequences.
  col2remove <- which(colSums(tmp_lik) == 0)
  if(length(col2remove) > 0) {
    tmp_lik <- tmp_lik[, -col2remove, drop=FALSE]
  }
  
  return(tmp_lik)
  
}


# Order the trait table to match the tree tip labels. Set all tips without a value to be NA.
unknown_tips_index <- which(! full_tree$tip.label %in% rownames(trait_values))
unknown_tips <- full_tree$tip.label[unknown_tips_index]
print(paste("Number of tips with unknown traits to be predicted: ", length(unknown_tips)))
assert_that(length(unknown_tips) > 0)

unknown_df <- as.data.frame(matrix(NA,
                                   nrow=length(unknown_tips),
                                   ncol=ncol(trait_values)))

rownames(unknown_df) = unknown_tips
colnames(unknown_df) = colnames(trait_values)

trait_values_ordered <- rbind(trait_values, unknown_df)

trait_values_ordered <- trait_values_ordered[full_tree$tip.label, , drop=FALSE]

num_tip <- nrow(trait_values_ordered)

if (args$method == "pic" | args$method == "scp" | args$method == "subtree_average") {
  
  if (args$method == "pic") {
    predict_out <- mclapply(trait_values_ordered,
                            hsp_independent_contrasts,
                            tree=full_tree,
                            weighted=TRUE,
                            check_input=FALSE,
                            mc.cores = args$threads)
    
  } else if (args$method == "scp") {
    
    predict_out <- mclapply(trait_values_ordered,
                            hsp_squared_change_parsimony,
                            tree=full_tree,
                            weighted=TRUE,
                            check_input=FALSE,
                            mc.cores = args$threads)
    
  } else if (args$method == "subtree_average") {
    
    predict_out <- mclapply(trait_values_ordered,
                            hsp_subtree_averaging,
                            tree = full_tree,
                            check_input = FALSE,
                            mc.cores = args$threads)
  }
  
  predicted_values <- mclapply(predict_out, function(x) { x$states[unknown_tips_index] }, mc.cores = args$threads)
  
} else if(args$method == "emp_prob" | args$method == "mp") {
  
  # Add 1 to all input counts because because traits states need to start at 1.
  trait_states_mapped <- trait_values_ordered + 1
  
  if (args$method == "emp_prob") {
    
    hsp_out_models <- mclapply(trait_states_mapped,
                               hsp_empirical_probabilities,
                               tree = full_tree,
                               check_input = FALSE,
                               mc.cores = args$threads)
    
  } else if (args$method == "mp") {
    
    hsp_out_models <- mclapply(trait_states_mapped,
                               hsp_max_parsimony,
                               tree = full_tree,
                               check_input = FALSE,
                               transition_costs = "proportional",
                               weight_by_scenarios = TRUE,
                               mc.cores = args$threads)
    
  }
  
  # Get subset of likelihood matrices for previously unknown tips only and output RDS file.
  num_unknown <- length(unknown_tips)
  
  hsp_out_models_unknown_lik <- mclapply(names(hsp_out_models), 
                                         function(x) { get_sorted_prob(hsp_out_models[[x]]$likelihoods,
                                                                       study_tips_i=unknown_tips_index, 
                                                                       tree_tips=full_tree$tip.label)},
                                         mc.cores = args$threads)
  
  names(hsp_out_models_unknown_lik) <- names(hsp_out_models)
  
  # Get state with highest probability in each case.
  predicted_values <- mclapply(hsp_out_models_unknown_lik,
                               function(x) { colnames(x)[max.col(x)] },
                               mc.cores = args$threads)  
}

# Format the predicted values as a data.frame
predicted_values <- data.frame(predicted_values, check.names = FALSE)

print("Checking that the input and output have the same number of traits")
assert_that(ncol(predicted_values) == ncol(trait_values))

# Add "sequence" as first column of predicted_values.
row.names(predicted_values) <- full_tree$tip.label[unknown_tips_index]
predicted_values <- predicted_values[, colnames(trait_values_ordered)]

# Add the predicted values to the previously known values
output_values <- rbind(trait_values[, colnames(predicted_values)], predicted_values)

# Write out predicted values in tabular format
write.table(output_values, file=args$output, row.names=TRUE, quote=FALSE, sep="\t")
