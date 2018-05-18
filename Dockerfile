FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
    apt-get install -y python3 python3-pip r-base bats

# R dependencies
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('Rcpp'); install.packages('assertthat'); install.packages('argparse'); install.packages('parallel'); install.packages('castor'); install.packages('ape')"

# Use /share as the working directory
RUN mkdir /share && mkdir /scratch
WORKDIR /share

# Add a wrapper to help execution via SciLuigi
RUN pip3 install bucket_command_wrapper==0.2.0 

# Add a wrapper script to do hidden state prediction with castor
ADD castor_hidden_state_prediction.R /usr/local/bin/

# Run tests and then remove the folder
ADD tests /usr/castor/tests
RUN bats /usr/castor/tests/ && rm -r /usr/castor/tests/
