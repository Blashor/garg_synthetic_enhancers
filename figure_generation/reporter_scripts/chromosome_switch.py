def replace_fasta_sequence(input_file, output_file, target_header, new_sequence):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        in_target = False

        for line in infile:
            if line.startswith('>'):
                if line.strip() == target_header:
                    in_target = True
                    outfile.write(line)  # write the matching header
                    # Write the new sequence wrapped at 60 characters
                    for i in range(0, len(new_sequence), 60):
                        outfile.write(new_sequence[i:i+60] + '\n')
                else:
                    in_target = False
                    outfile.write(line)
            else:
                if not in_target:
                    outfile.write(line)
                # skip original sequence if in_target is True

#Fill in input and output file names and specify which chromosome is being replaced
input_fasta = "GRCm38.primary_assembly.genome.fa"
output_fasta = "custom_genome.fa"
header_to_replace = ">chr11 11"  #or whichever chromosome is being swapped
with open("chr11_edited.txt") as f: #replace with edited txt file of chromosome
    new_seq = f.read().replace("\n", "")

replace_fasta_sequence(input_fasta, output_fasta, header_to_replace, new_seq)