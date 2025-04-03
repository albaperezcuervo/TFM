#!/usr/bin/python
import os
import csv

# Retrieve all .fastq.gz files in the current directory
files = [f for f in os.listdir() if f.endswith('.fastq.gz')]

# Open or create 'manifest.csv' for writing
with open('manifest.csv', 'w') as file:
    writer = csv.writer(file)
    # Write the header row
    writer.writerow(["sample-id", "absolute-filepath", "direction"])

    # Iterate over each file to determine its direction and write the details to CSV
    for j in files:
        if '_1.fastq.gz' in j:
            name = j.split('_')[0]  # Extract sample ID
            path = os.path.join(os.getcwd(), j)  # Get absolute file path
            direction = 'forward'
            writer.writerow([name, path, direction])
            print(f"Adding: {name}, {path}, {direction}")  # Print details for verification

        elif '_2.fastq.gz' in j:
            name = j.split('_')[0]  # Extract sample ID
            path = os.path.join(os.getcwd(), j)  # Get absolute file path
            direction = 'reverse'
            writer.writerow([name, path, direction])
            print(f"Adding: {name}, {path}, {direction}")  # Print details for verification