import random

result_true = []


def mutation_bag(chance, num_times):
    mut_obj = {"SNPs": 0, "dups": 0}
    for i in range(num_times):
        res = random.random() < chance
        if res == True:
            mut_obj["dups"] += 1
        else:
            mut_obj["SNPs"] += 1
    return mut_obj


# Example usage

print(random.sample(range(100), 1))
