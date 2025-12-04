import os
import itertools
import sys

from PIL import Image

#
# Aggregates basenji generated SAT files by selected generations
#

working_dir = "/Users/blake/Downloads/sat_output"

combo_dir = "loss_file_together"
os.system(f"mkdir -p {combo_dir}")
samples = ["line256rep6", "line256rep28", "line256rep31"]
for sample in samples:
    gens = list(map(lambda x: f"g{x}", [1, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300]))
    #
    seqs_to_test = set(map(lambda x: "_".join(x), list(itertools.product(gens, [sample]))))

    os.system(f"mkdir -p {combo_dir}/{sample}")

    for f in seqs_to_test:
        images = [
            Image.open(x) for x in [f"{working_dir}/{f}/plots/seq0_t1.png", f"{working_dir}/{f}/plots/seq1_t1.png"]
        ]
        images = [im.crop((500, 0, im.size[0] - 2000, im.size[1])) for im in images]
        widths, heights = zip(*(i.size for i in images))
        print(widths, heights)
        total_width = sum(widths)
        max_height = max(heights)

        new_im = Image.new("RGB", (total_width, max_height))

        x_offset = 0
        for im in images:
            new_im.paste(im, (x_offset, 0))
            x_offset += im.size[0]

        new_im.save(f"{combo_dir}/{sample}/{f}.png")
