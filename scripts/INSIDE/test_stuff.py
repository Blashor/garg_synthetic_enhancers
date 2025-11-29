def shrinking_window(seq):
    seq_left, seq_right = seq[: len(seq) // 2], seq[len(seq) // 2 :]
    new_seqs = []
    print(seq_left, seq_right, "\n")
    for small_seq, no_change, reverse, flip_seqs in [
        (seq_left, seq_right, False, False),
        (seq_left, seq_right, True, False),
        (seq_right, seq_left, False, True),
        (seq_right, seq_left, True, True),
    ]:
        shrink_seq = ""
        ran_into_T = True
        if reverse == True:
            small_seq = small_seq[::-1]
        for i, base in enumerate(small_seq):
            if ran_into_T == True:
                if base == "T":
                    # print(i, base)
                    shrink_seq += base
                else:
                    # print(i, "T")
                    shrink_seq += "T"
                    ran_into_T = False
            elif ran_into_T == False:
                # print(i, base)
                shrink_seq += base
        if reverse == True:
            shrink_seq = shrink_seq[::-1]
        if flip_seqs == False:
            print(shrink_seq, no_change)
            new_seqs.append(shrink_seq + no_change)
        else:
            print(no_change, shrink_seq)
            new_seqs.append(no_change + shrink_seq)
    print(new_seqs)
    return new_seqs


shrinking_window("TACGCGTTTACGTTTT")
