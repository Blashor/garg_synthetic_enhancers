# this will read dna files
# create fastas for each corresponding demux from dna file
# align

import subprocess
import os
import sys


def IUPAC_codes(s, s2, fa):
    sam_file = pysam.AlignmentFile(s, "r")
    out = pysam.AlignmentFile(s2, "wb", template=sam_file)
    fasta = pysam.FastaFile(fa)

    for read in sam_file.fetch():
        if read.is_unmapped:
            # out.write(read)
            continue
        if read.query_sequence != None:
            y_list = list(read.query_sequence)
            r_list = list(read.query_sequence)
            Y = 0
            R = 0
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                sam_n = read.query_sequence[query_pos].upper()
                fasta_n = fasta.fetch(read.reference_name, ref_pos, ref_pos + 1).upper()
                if fasta_n == "C" and sam_n == "T":
                    y_list[query_pos] = "Y"
                    Y += 1
                if fasta_n == "G" and sam_n == "A":
                    r_list[query_pos] = "R"
                    R += 1
            denom = Y + R
            if denom != 0:
                y_prop = Y / denom
                if y_prop > 0.9:
                    read.query_sequence = "".join(y_list)
                elif y_prop < 0.1:
                    read.query_sequence = "".join(r_list)
            out.write(read)
    # print(Y, R, s)
    out.close()
    sam_file.close()
    fasta.close()


def pipeline(sysarg):
    demux_info = "/home/bmt26/garg/Atac_reporter/demux_files/fastas"
    fastq_to_fasta = {
        # "g0": "mm10_SS_6_1_G1.fa",
        # "g30": "mm10_SS_6_1_G30.fa",
        # "g90": "mm10_SS_6_1_G90.fa",
        "g120": "mm10_SS_6_1_G120.fa",
        "g150": "mm10_SS_6_1_G150.fa",
    }

    fastq_to_fasta = {
        "0nM": "mm10_SS_6_1_G150.fa",
        "200nM": "mm10_SS_6_1_G150.fa",
        "400nM": "mm10_SS_6_1_G150.fa",
        "2uM": "mm10_SS_6_1_G150.fa",
        "4uM": "mm10_SS_6_1_G150.fa",
    }
    fastq_to_fasta = {
        "G1.2": "mm10_SS_6_1_G1.fa",
        "G30.2": "mm10_SS_6_1_G30.fa",
        "G60.2": "mm10_SS_6_1_G60.fa",
        "G90.0": "mm10_SS_6_1_G90.fa",
        "G90.2": "mm10_SS_6_1_G90.fa",
        "G90.4": "mm10_SS_6_1_G90.fa",
        "G120.2": "mm10_SS_6_1_G120.fa",
        "G150.2": "mm10_SS_6_1_G150.fa",
    }
    fastq_to_fasta = {
        "sep160": "G1_G60.fa",
        "sep601": "G60_G1.fa",
        "sep901": "G90_G1.fa",
        "sep190": "G1_G90.fa",
        "sepr6g1_0": "mm10_SS_6_1_G1.fa",
        "sepr6g1": "mm10_SS_6_1_G1.fa",
        "sepr6g60_0": "mm10_SS_6_1_G60.fa",
        "sepr6g60": "mm10_SS_6_1_G60.fa",
        "sepr6g90": "mm10_SS_6_1_G90.fa",
    }
    fastq_to_fasta = {
        "plG1": "mm10_SS_6_1_G1.fa",
        "plG60": "mm10_SS_6_1_G60.fa",
        "plG90": "mm10_SS_6_1_G90.fa",
        "plG60NE": "mm10_SS_6_1_G60.fa",
        "plG90NE": "mm10_SS_6_1_G90.fa",
        "jk_r6g1": "mm10_SS_6_1_G1.fa",
        "jk_r6g30": "mm10_SS_6_1_G30.fa",
        "jk_r6g60": "mm10_SS_6_1_G60.fa",
        "jk_r6g90": "mm10_SS_6_1_G90.fa",
        "jk_r6g120": "mm10_SS_6_1_G120.fa",
        "jk_r6g150": "mm10_SS_6_1_G150.fa",
        "jk_r6g90p58": "mm10_SS_6_1_G90.fa",
        "jk_rg6g60NE": "mm10_SS_6_1_G60.fa",
        "jk_rg6g90NE": "mm10_SS_6_1_G90.fa",
    }
    fastq_to_fasta = {
        "shh_1": "shh_one.fasta",
        "shh_20": "shh_twenty.fasta",
        "shh_50": "shh_fifty.fasta",
        "shh_1ne": "shh_one.fasta",
        "shh_20ne": "shh_one.fasta",
    }
    fastq_to_fasta = {
        "shhg1": "shh_one.fasta",
        "shhg1_NE": "shh_one.fasta",
        "shhg10_5": "shh_ten.fasta",
        "shhg10_3": "shh_ten.fasta",
        "shhg10_NE5": "shh_ten.fasta",
        "shhg10_NE3": "shh_ten.fasta",
        "shhg20": "shh_twenty.fasta",
        "shhg50": "shh_fifty.fasta",
    }
    fastq_to_fasta = {
        "meox1_g0": "meox1_one.fasta",
        "meox1_g10": "meox1_ten.fasta",
        "meox1_g50": "meox1_fifty.fasta",
    }

    fastq_to_fasta = {
        "mar160": "G1_G60.fa",
        "mar601": "G60_G1.fa",
        "mar901": "G90_G1.fa",
        "mar190": "G1_G90.fa",
        "marg1": "mm10_SS_6_1_G1.fa",
        "marg60": "mm10_SS_6_1_G60.fa",
        "marg90": "mm10_SS_6_1_G90.fa",
    }
    fastq_to_fasta = {
        "meox_0_mar15": "meox1_one.fasta",
        "meox_10_mar15": "meox1_ten.fasta",
        "meox_50_mar15": "meox1_fifty.fasta",
    }
    fastq_to_fasta = {
        "shh0_mar21": "shh_one.fasta",
        "shh20_mar21": "shh_twenty.fasta",
        "shh50_mar21": "shh_fifty.fasta",
    }

    fastq = "/home/bmt26/garg/blake/daf_seq/fastq"
    sam_raw = "/home/bmt26/garg/blake/daf_seq/sam_preIUPAC"
    sam_final = "/home/bmt26/garg/blake/daf_seq/sam"
    bam_m6a = "/home/bmt26/garg/blake/daf_seq/bam"
    bam_nuc = "/home/bmt26/garg/blake/daf_seq/bam_nuc"
    #
    os.makedirs(sam_raw, exist_ok=True)
    os.makedirs(sam_final, exist_ok=True)
    os.makedirs(bam_m6a, exist_ok=True)
    os.makedirs(bam_nuc, exist_ok=True)
    config_out = []
    for sample in fastq_to_fasta:
        fa = f"{demux_info}/{fastq_to_fasta[sample]}"
        fq = f"{fastq}/{sample}.fastq"
        #
        sam_temp = f"{sam_raw}/{sample}.sam"
        sam = f"{sam_final}/{sample}.bam"
        m6a = f"{bam_m6a}/{sample}.bam"
        nuc = f"{bam_nuc}/{sample}.bam"
        # -Y
        mm2 = f"minimap2 --MD -Y -a -y -x map-pb {fa} {fq} > {sam_temp}"
        if sysarg == "ft":
            os.system(mm2)
            IUPAC_codes(sam_temp, sam, fa)
            os.system(f"ft ddda-to-m6a {sam} {m6a}")
            os.system(f"ft add-nucleosomes -n 60 -c 70 --min-distance-added 15 -d 10 {m6a} > {nuc}")
            os.system(f"samtools sort {nuc} -o {nuc}.sorted.bam")
            os.system(f"samtools index {nuc}.sorted.bam")

            os.system(f"ft extract {nuc} --m6a {nuc}.m6a.bed")
            os.system(f"ft extract {nuc} --msp {nuc}.msp.bed")
            os.system(f"ft extract {nuc} --nuc {nuc}.nuc.bed")
        elif sysarg == "fire":
            with open("/home/bmt26/garg/FIRE_snakemake/FIRE/config/config.yaml", "w") as file:
                file.writelines([f"ref: {fa}\nref_name: mm10\nmanifest: config/config.tbl\nmax_t: 32\n"])
            with open("/home/bmt26/garg/FIRE_snakemake/FIRE/config/config.tbl", "w") as file:
                file.writelines(["sample\tbam\n", f"{sample}\t{nuc}.sorted.bam\n"])
            os.system("./fire --configfile /home/bmt26/garg/FIRE_snakemake/FIRE/config/config.yaml --unlock")
            os.system("./fire --configfile /home/bmt26/garg/FIRE_snakemake/FIRE/config/config.yaml")


"""
if sys.argv[1] == "ft":
    import pysam
pipeline(sys.argv[1])
"""
import pysam

pipeline("ft")
