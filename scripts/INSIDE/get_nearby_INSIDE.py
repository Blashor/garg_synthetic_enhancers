from collections import defaultdict


def parse_region(region_str):
    chrom, coords = region_str.replace(",", "").split(":")
    start, end = map(int, coords.split("-"))
    return chrom, start, end


def region_to_str(region_tuple):
    chrom, start, end = region_tuple
    return f"{chrom}:{start}-{end}"


def within_20kb(region1, region2):
    c1, s1, e1 = region1
    c2, s2, e2 = region2
    if c1 != c2:
        return False
    return abs(s1 - e2) <= 20000 or abs(s2 - e1) <= 20000


# Your gene list (name, region)
gene_regions_raw = [
    ("TBX1", "chr16:18,588,019-18,588,901"),
    ("SHH", "chr5:28,466,644-28,467,072"),
    ("LIS1", "chr11:74,724,100-74,724,371"),
    ("SALL1", "chr8:89,042,707-89,042,847"),
    ("TWIST1", "chr12:33,957,671-33,958,006"),
    ("RAI1", "chr11:60,105,009-60,105,396"),
    ("TCF4", "chr18:69,344,509-69,344,958"),
    ("KMT2D", "chr15:98,871,173-98,871,212"),
    ("CHD7", "chr4:8,690,415-8,690,611"),
    ("ASXL1", "chr2:153,345,834-153,345,872"),
    ("NIPBL", "chr15:8,444,150-8,444,388"),
    ("ARID1B", "chr17:4,995,075-4,996,491"),
    ("ZEB2", "chr2:45,110,346-45,110,384"),
    ("EYA1", "chr1:14,309,671-14,310,107"),
    ("SOX9", "chr11:112,782,209-112,782,249"),
]

gene_regions = [(name, parse_region(region)) for name, region in gene_regions_raw]

# Dictionary to collect matches by gene
matches_by_gene = defaultdict(list)
name_counter = defaultdict(int)

# Input file
input_file = "/Users/blake/Documents/gargLab/Figures_for_Box/Figure_1_Ranking/data/3runs_ce_rank.txt"

with open(input_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        try:
            enh1 = parse_region(parts[0])
            enh2 = parse_region(parts[1])
        except ValueError:
            continue

        for gene_name, gene_region in gene_regions:
            if within_20kb(enh1, gene_region) or within_20kb(enh2, gene_region):
                name_counter[gene_name] += 1
                label = f"{gene_name}-{name_counter[gene_name]}" if name_counter[gene_name] > 1 else gene_name
                matches_by_gene[gene_name].append((label, region_to_str(enh1), region_to_str(enh2)))
                break

# Sort and print
print("enhs = [")
for gene in sorted(matches_by_gene.keys()):
    for entry in matches_by_gene[gene]:
        print(f"    {entry},")
print("]")
