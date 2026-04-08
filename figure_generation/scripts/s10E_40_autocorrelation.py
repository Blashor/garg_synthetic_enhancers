import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import INSIDE_1000_enhs
import random


# lines.append(f"{chrom1}\t{begin}\t{ending}\n")
def seq_to_gc_binary(seq):
    return np.array([1 if b in "GC" else 0 for b in seq])


def scramble_string(s, seed=None):
    if seed is not None:
        random.seed(seed)

    # split into 2-bp chunks
    chunk = 1
    pairs = [s[i : i + chunk] for i in range(0, len(s), chunk)]

    # shuffle the chunks
    random.shuffle(pairs)

    # rejoin
    return "".join(pairs)


def upsample(signal, factor=2):
    return np.repeat(signal, factor)


def seq_to_S_binary(b, filters):
    # diffs = []
    diffs = np.zeros(len(b))
    # filters = {"AA", "TT", "AT", "TA"}
    # filters = {"G", "C"}
    chunk = len(next(iter(filters)))
    for i in range(0, len(b)):
        # if 65 < random.randint(0, 100):
        #    diffs.append((i, "G_C"))
        dinuc = "".join(b[i : i + chunk])
        # filters = {"TTT", "AAA", "TAT", "AAT", "ATT", "TTA", "TAA", "ATA"}
        # filters = {"TTTT", "AAAA"}  # {"AA", "TT", "AT", "TA"}
        # filters = {"G", "C"}
        # filters = {"G", "C"}
        if dinuc in filters:
            diffs[i] = 1

    return np.array((diffs))


def seq_to_W_binary(b, filters):
    # diffs = []
    diffs = np.zeros(len(b))
    # filters = {"AA", "TT", "AT", "TA"}
    # filters = {"G", "C"}
    chunk = len(next(iter(filters)))
    for i in range(0, len(b)):
        # if 65 < random.randint(0, 100):
        #    diffs.append((i, "G_C"))

        dinuc = "".join(b[i : i + chunk])
        #
        # filters = {"TTT", "AAA", "TAT", "AAT", "ATT", "TTA", "TAA", "ATA"}

        # filters = {"GG", "CC", "GC", "CG"}
        # filters = {"G", "C"}
        # filters = {"GG", "CC", "GC", "CG"}
        if dinuc in filters:
            diffs[i] = 1

    return np.array((diffs))


# autocorrelation
def pearson_autocorr(seq, max_lag=None, filter_set=None):
    x = np.array(seq_to_S_binary(seq[:max_lag], filter_set))
    y = np.array(seq_to_W_binary(seq[:max_lag], filter_set))  # renamed to avoid overwrite
    N = len(x)
    if max_lag is None:
        max_lag = N - 1
    result = []
    mean_x = np.mean(x)
    std_x = np.std(x)
    mean_y = np.mean(y)
    std_y = np.std(y)

    for lag in range(max_lag + 1):
        if lag == 0:
            result.append(1.0)  # correlation with itself
        else:
            x1 = x[:-lag] - mean_x
            y1 = y[lag:] - mean_y
            corr = np.sum(x1 * y1) / (len(x1) * std_x * std_y)
            result.append(corr)
    return np.array(result)


from scipy.signal import correlate
import numpy as np


def string_diff(a, b):
    diffs = []

    for i in range(1, len(b)):
        # if 65 < random.randint(0, 100):
        #    diffs.append((i, "G_C"))
        dinuc = "".join(b[i - 1 : i + 1])

        if dinuc == "AA" or dinuc == "AT" or dinuc == "TA":
            diffs.append((i, "G_C"))

    return diffs


def set_style():
    fig = plt.figure()
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


enh_strings = np.zeros([40, 300, 2]).tolist()
with open("../data/INSIDE_jun13.txt") as file:
    for line in file:
        cell = line.split("\t")
        # print(cell)
        enh_num = int(cell[1].split("rep")[1]) - 1

        gen_num = int(cell[0].split(": ")[-1]) - 1
        if gen_num < 300:
            e1 = cell[2][:512]
            e2 = cell[2][512:]

            enh_strings[enh_num][gen_num][0] = e1
            enh_strings[enh_num][gen_num][1] = e2
        # print(cell[0], len(e1), len(e2))
        # print(cell)
        # print(cell)
    # print(enhs)


scale_size = 862

change_types = ["G_C"]
current_diff_across = {c: np.zeros(scale_size) for c in change_types}
total_coverage = np.zeros(scale_size)
max_len = 0


def multi_sequence_periodicity(seqs, max_lag=512, label=0, filter_set=None):
    total_corr = np.zeros(max_lag)
    total_power = np.zeros(max_lag)
    n = 0

    for seq in seqs:
        # mid = len(seq) // 2
        # seq = seq[mid - max_lag // 2 : mid + max_lag // 2]

        # print(signal)
        # max_lag = min(max_lag, len(signal))
        # pearson_autocorr(signal)
        # autocorr = np.correlate(signal, signal, mode="full")  # full returns 2N-1 values
        # autocorr = autocorr[autocorr.size // 2 :]
        # print(autocorr[:max_lag])
        corr = pearson_autocorr(seq, max_lag, filter_set)[:max_lag]
        if np.isnan(corr).any() == False:
            total_corr += corr
            n += 1

            # centered_signal = signal - np.mean(signal)
            # fft_vals = np.fft.fft(centered_signal, n=max_lag)
            # power = np.abs(fft_vals) ** 2
            # total_power += power

        # print(signal)
        # max_lag = min(max_lag, len(signal))
        # pearson_autocorr(signal)
        # autocorr = np.correlate(signal, signal, mode="full")  # full returns 2N-1 values
        # autocorr = autocorr[autocorr.size // 2 :]

    avg_corr = total_corr / n
    avg_power = total_power / n

    lags = np.arange(max_lag)
    freqs = np.fft.fftfreq(max_lag, d=1)  # 1 bp per step
    periods = np.zeros_like(freqs)
    nonzero = freqs != 0
    periods[nonzero] = 1 / freqs[nonzero]

    # Keep only positive periods
    positive = periods > 0

    # --- Plot FFT power spectrum ---
    linestyle = "--"
    if "random" in label:
        return avg_corr
        alpha = int(label.split("_")[0])
    else:
        linestyle = "-o"
        alpha = int(label)
    print(n)
    plt.plot(
        np.array(range(len(avg_corr))),
        avg_corr,
        linestyle,
        label=f"{alpha+1}",
        color="black",
        alpha=(alpha + 50) / 350,
        markersize=2,
    )
    # plt.plot(periods[positive], avg_power[positive], label=label, color="black", alpha=(label + 10) / 60)


def motif_pair_positions_degenerate(seq, motif_set1, motif_set2, lag):
    motif_len = len(next(iter(motif_set1)))  # assume all motifs same length
    return [
        i
        for i in range(len(seq) - motif_len - lag - motif_len + 1)
        if seq[i : i + motif_len] in motif_set1 and seq[i + lag : i + lag + motif_len] in motif_set2
    ]


for f_name, filter_set in (("WW", {"AT", "AA", "TA", "TT"}), ("S", {"G", "C"})):
    set_style()
    seqs0 = []
    seqs50 = []

    for i in range(0, 301, 50):
        all_pos = []
        gen = max(0, i - 1)
        seqs = []
        for ei, e in enumerate(enh_strings):
            e1_len, e2_len = (512, 512)
            max_len = max(max_len, e1_len, e2_len)
            # generation 0 vs generation 50

            seq_1 = e[gen][0]
            seq_2 = e[gen][1]

            seqs.append(seq_1)
            seqs.append(seq_2)
        multi_sequence_periodicity(seqs, label=f"{gen}", filter_set=filter_set)
        # for seq in seqs:
        #    positions = motif_pair_positions_degenerate(seq, {"GC", "CG", "CC", "GG"}, {"AA", "TT", "AT", "TA"}, 16)
        #    for p in positions:
        #        all_pos.append(p)
        # sns.histplot(all_pos, label=f"{gen}")
    # plt.legend(loc="center right")
    # plt.show()
    #
    scrambles = []
    for w in range(17):
        for i in range(0, 301, 50):
            gen = max(0, i - 1)
            seqs = []
            for ei, e in enumerate(enh_strings):
                e1_len, e2_len = (512, 512)
                max_len = max(max_len, e1_len, e2_len)
                # generation 0 vs generation 50

                seq_1 = e[gen][0]
                seq_2 = e[gen][1]

                seqs.append(scramble_string(seq_1))
                seqs.append(scramble_string(seq_1))

            scrambles.append(multi_sequence_periodicity(seqs, label=f"{gen}_random", filter_set=filter_set))
    scrambles = np.array(scrambles[:100])
    print(scrambles.shape)
    scrambles = np.percentile(scrambles, [2.5, 97.5], axis=0)
    print(scrambles.shape)
    plt.fill_between(np.array(range(0, len(scrambles[1]))), scrambles[0], scrambles[1], color="black", alpha=0.5)
    plt.plot(
        np.array(range(0, len(scrambles[0]))),
        scrambles[0],
        "--",
        label="Randomized Sequences",
        color="black",
        alpha=1,
    )
    plt.plot(
        np.array(range(0, len(scrambles[1]))),
        scrambles[1],
        "--",
        color="black",
        alpha=1,
    )
    # plt.axhline(0)
    plt.legend(loc="upper right")
    plt.ylim([-0.15, 0.20])
    plt.xlim([0, 45])
    plt.savefig(f"../figures/s10e_40_{f_name}_ac.svg")
# plt.show()
# print(max_len)
