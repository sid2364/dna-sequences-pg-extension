"""
For parsing data files downloaded from National Library of Medicine (for e.g., https://www.ncbi.nlm.nih.gov/nuccore/HQ287898.1)
"""
def extract_dna_sequences(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        first_line_written = False  # To ensure all output is on a single line
        for line in infile:
            if line.startswith('N'):
                sequence = line[1:].strip()
                if not first_line_written:
                    outfile.write(sequence)
                    first_line_written = True
                else:
                    outfile.write(sequence)
        outfile.write('\n')

# Example usage
input_file = 'SAMN01780187.fastq'
output_file = 'dna.txt'
#extract_dna_sequences(input_file, output_file)

"""
Generating random DNA sequences of the specified length
"""
import random

def generate_random_nucleotides(output_file, num_nucleotides=10000000):
    nucleotides = ['A', 'T', 'C', 'G']

    with open(output_file, 'w') as outfile:
        # Generate random nucleotides and write them directly to the file
        for _ in range(num_nucleotides):
            outfile.write(random.choice(nucleotides))
        outfile.write('\n') 

output_file = 'random_nucleotides_10M_100Mb.txt'
generate_random_nucleotides(output_file, 10_000_000)

output_file = 'random_nucleotides_200K_2Mb.txt'
generate_random_nucleotides(output_file, 200_000)

output_file = 'random_nucleotides_1M_10Mb.txt'
generate_random_nucleotides(output_file, 1_000_000)

output_file = 'random_nucleotides_100K_1Mb.txt'
generate_random_nucleotides(output_file, 100_000)

output_file = 'random_nucleotides_1K_100Kb.txt'
generate_random_nucleotides(output_file, 1000)

# Side note: place the output files in /tmp/ or some other place (not home dir), otherwise postgres cannot access them!
