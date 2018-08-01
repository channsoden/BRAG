# Use an official Python runtime as a parent image
FROM python:3.6.4

# Create BRAG directory
RUN mkdir /home/BRAG

# Set the working directory
WORKDIR /home/BRAG

# Copy the current directory contents into the container working directory
ADD . /home/BRAG

# Install any Python dependencies specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt


# Set the working directory for the sample analysis
WORKDIR /home/BRAG/sample_data

# Create a directory for the results
RUN mkdir /home/BRAG/sample_data/results

# Run sample_data/sample_analysis.py when the container launches
CMD exec bash sample_analysis.sh 1> results/sample_analysis.out 2> results/sample_analysis.err
