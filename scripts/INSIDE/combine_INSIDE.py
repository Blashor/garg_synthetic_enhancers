import json
import os


def read_sbatch_splits():
    gen_obj = {}
    line_obj = {}
    lines = []
    how_many_reps_in_gen = 60  # Run twice, first time to see what to set this to
    gen_max = 151
    s_file_combine = "INSIDE_genome_sep12_75gens.txt"
    s_file_combine = "INSIDE_jun13_160gens.txt"
    # s_files in rep order
    s_files = [
        "/Users/blake/Downloads/INSIDE_split_0.out",
        "/Users/blake/Downloads/INSIDE_split_1.out",
        "/Users/blake/Downloads/INSIDE_split_2.out",
        "/Users/blake/Downloads/INSIDE_split_3.out",  # first splits
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34128895.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34128896.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34128897.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34128898.out",  # up to 250
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211463.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211464.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211465.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211466.out",  # 250 only
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211694.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211695.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211696.out",
        "/Users/blake/Downloads/INSIDE_split/INSIDE_34211697.out",  # current
    ]
    s_file_combine = "/Users/blake/Downloads/INSIDE_margin.txt"
    s_files = [
        "/Users/blake/Downloads/INSIDE_margin1.out",
        "/Users/blake/Downloads/INSIDE_margin2.out",
        "/Users/blake/Downloads/INSIDE_margin3.out",
        "/Users/blake/Downloads/INSIDE_margin4.out",
        "/Users/blake/Downloads/INSIDE_margin5.out",
        "/Users/blake/Downloads/INSIDE_margin6.out",
    ]
    """
    s_files = [
        "/Users/blake/Documents/gargLab/INSIDE/INSIDE_gen80_rep40.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34354680.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34354701.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34354728.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34354745.out",  #
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34533312.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34533315.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34533316.out",
        "/Users/blake/Downloads/INSIDE_splitB/INSIDE_34533320.out",
    ]
    s_files = list(
        map(
            lambda x: f"/Users/blake/Documents/gargLab/INSIDE_shrinking_window/{x}",
            sorted(os.listdir("/Users/blake/Documents/gargLab/INSIDE_shrinking_window")),
        )
    )
    """

    # s_file_combine = "INSIDE_combine.txt"
    # s_file_combine = "INSIDE_shrinking_window.txt"

    for s in s_files:
        last_gen = "N/A"
        print(s)
        if "INSIDE" in s.split("/")[-1]:
            with open(s) as file:
                for line in file:
                    if line.startswith("Gen:"):
                        cell = line.split("\t")
                        if cell[1] not in line_obj:
                            line_obj[cell[1]] = []
                        line_obj[cell[1]].append(line)
                    else:
                        pass
                        # print(line)
    #
    # Fixes wonky generations
    #
    num_survivors = 4
    for line in line_obj:
        current_gen = 0  # will immediately add 1
        for s_i, sample_in_time in enumerate(line_obj[line]):
            if s_i % num_survivors == 0:
                current_gen += 1
            cell = sample_in_time.split("\t")
            cell[0] = f"Gen: {current_gen}"
            new_line = "\t".join(cell)
            if cell[0] not in gen_obj:
                gen_obj[cell[0]] = []
            gen_obj[cell[0]].append(new_line)
    for gen in gen_obj:
        gen_len = len(gen_obj[gen])
        print(gen, gen_len)
        if int(gen.split(": ")[-1]) < gen_max:
            if gen_len == how_many_reps_in_gen:
                lines += gen_obj[gen]
    # print(lines)
    print(len(lines))

    with open(s_file_combine, "w") as file:
        file.writelines(lines)


read_sbatch_splits()
