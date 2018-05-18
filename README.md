# docker-castor
Docker image running castor

Castor is a set of utilities for dealing with large phylogenetic
treest, and documentation can be found 
[here](https://cran.r-project.org/web/packages/castor/castor.pdf).

#### Scripts:

`castor_hidden_state_prediction.R`: Utility for performing hidden state
prediction using `castor`. The inputs are (1) a Newick tree (2) a tab-
delimited file with states known for leaves of that tree and (3) a path
for the output table with predicted hidden states. Note that large portions
of the code for this script was found in `castor_hsp.R` within `picrust/picrust2` 
