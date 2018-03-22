# Use an official Python runtime as a parent image
FROM python:3

# Create BRAG directory
RUN mkdir /home/BRAG

# Set the working directory
WORKDIR /home/BRAG

# Copy the current directory contents into the container working directory
ADD . /home/BRAG

# Install any Python dependencies specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Run sample_data/sample_analysis.py when the container launches
#CMD ["bash", "./sample_data/sample_analysis.sh"]
