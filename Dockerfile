FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
    apt-get install -y python3 python3-pip r-base bats

# Use /share as the working directory
RUN mkdir /share && mkdir /scratch
WORKDIR /share

# Add a wrapper to help execution via SciLuigi
RUN pip3 install bucket_command_wrapper==0.2.0 

# Run tests and then remove the folder
ADD tests /usr/castor/tests
RUN bats /usr/castor/tests/ && rm -r /usr/castor/tests/
