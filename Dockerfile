# Use an official Python runtime as a parent image
FROM python:3

# Create BRAG directory
RUN mkdir /home/BRAG

# Copy the current directory contents into the container working directory
ADD . /home/BRAG

# Set the working directory
WORKDIR /home/BRAG/sample_data

# Install any Python dependencies specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Create a directory for the results
RUN mkdir /home/BRAG/sample_data/results

# Run sample_data/sample_analysis.py when the container launches
CMD "bash sample_analysis.sh 1> sample_data/results/sample_analysis.out 2> sample_data/results/sample_analysis.err"
