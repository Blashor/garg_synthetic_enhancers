def make_pairs(boundaries, n_samples):
    pairs = []
    end_pairs = []
    last_b = -1
    for b in boundaries:
        for i in range(last_b + 1, b):
            if i + 1 < b:
                pairs.append((i, i + 1))
            if i + 1 == b:
                end_pairs.append((i, i + 1))
        last_b = b
    t_indices, t_plus_1_indices = zip(*pairs)
    t_indices_valid, t_plus_1_indices_valid = zip(*end_pairs)
    return list(t_indices), list(t_plus_1_indices), list(t_indices_valid), list(t_plus_1_indices_valid)


boundaries = [6, 12, 19, 26, 32, 39, 44, 50, 53, 56, 59, 62]  # example boundaries
n_samples = 63
train_indices = make_pairs(boundaries, n_samples)
print(train_indices[2], train_indices[3])
