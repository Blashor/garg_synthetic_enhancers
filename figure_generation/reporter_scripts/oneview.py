#Bigwig files have to be named as sample_name_1.reporterfastafilename.bw <-- swap out "sample" and "reporterfastafilename" for each sample

import pyBigWig
import os

def write_chr1(entries, out_file):
    # Open a new BigWig file for writing
    bw_out = pyBigWig.open(out_file, "w")
    print(chrom_length)
    # Add chromosome definition for 'chr1' with the specified length
    bw_out.addHeader([("chr1", chrom_length)])

    # Check if there are any entries
    print(entries)
    starts = [e[0] for e in entries]
    ends = [e[1] for e in entries]
    values = [e[2] for e in entries]
    print(starts)
    print(ends)
    print(values)
    c = ["chr1" for e in entries]
    # Write all entries under 'chr1'
    bw_out.addEntries(c, starts, ends=ends, values=values)

    # Close the BigWig file
    bw_out.close()


def get_all_entries(input_file, scaffold_id):
    # Open the BigWig file
    bw = pyBigWig.open(input_file)
    print(bw)
    global chrom_length
    # Get all chromosomes in the BigWig file
    chromosomes = bw.chroms()
    print(chromosomes)
    # print(len(chromosomes), chromosomes)
    """
    for chrom, length in chromosomes.items():
        print(chrom)

        if chrom in i5_to_DNA_name:
            good_chrom = chrom
            dna_name = i5_to_DNA_name[chrom]
            print(i5_to_DNA_name[chrom], input_file)
    """
    # Initialize a dictionary to store all original entries
    original_entries = {}
    print(original_entries)
    # Iterate over each chromosome
    for chrom, length in chromosomes.items():
        if str(chrom) == "chr1":
            chrom_length = length
            print(chrom_length)
        if chrom == scaffold_id:
            # Retrieve all entries for the current chromosome
            entries = bw.intervals(chrom, 0, length)
            print(entries)
            # Store original entries in the dictionary
            original_entries[chrom] = entries
            print(original_entries[chrom])

    # Close the BigWig file
    bw.close()

    return scaffold_id, entries

def main():
    directory = "/main/directory"
    os.system(f"mkdir -p {directory}/bws_onescreen")
    # i5_to_DNA_name = {}
    """
    with open(f"{demux_info}/{atac_table}", encoding="latin-1") as file:
        for line in file:
            cell = line.strip().split("\t")
            i5 = str(Seq(cell[4]).reverse_complement())
            scaffold_id = cell[6].split(" ")[0].split(".dna")[0]
            i5_to_DNA_name[i5] = scaffold_id
    """

    bws = os.listdir("/main/directory/bigwig")
    for bw in bws:
        print(bw)
        scaffold_id = bw.split(".")[1]
        print(scaffold_id)
        input_file = f"/main/directory/bigwig/{bw}"
        print(input_file)
        dna_name, entries = get_all_entries(input_file, scaffold_id)
        write_chr1(entries, f"{directory}/bws_onescreen/{bw}")


main()