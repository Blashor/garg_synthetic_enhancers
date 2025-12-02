#
# FIRE ce rank file to bedgraph
#
bg_dict = {}
with open("../data/GM12878_ce_rank.txt") as file:
    for line in file:
        print(line)
        cell = line.split("\t")
        chromA, coordsA = cell[0].split(":")
        startA, endA = map(int, coordsA.split("-"))
        chromB, coordsB = cell[1].split(":")
        startB, endB = map(int, coordsB.split("-"))
        score = float(cell[3])
        cA = f"{chromA}\t{startA}\t{endA}"
        cB = f"{chromB}\t{startB}\t{endB}"
        if cA not in bg_dict:
            bg_dict[cA] = score
        else:
            bg_dict[cA] = max(score, bg_dict[cA])
        if cB not in bg_dict:
            bg_dict[cB] = score
        else:
            bg_dict[cB] = max(score, bg_dict[cB])
bg_lines = []
for chrom in bg_dict:
    bg_lines.append(f"{chrom}\t{bg_dict[chrom]}\n")
with open("../data/GM12878_ce_rank.bedgraph", "w") as file:
    file.writelines(bg_lines)

# rgb(157, 36, 36)
"""

var fiber_window = document.querySelector("[clip-path='url(#fire_bed12bed_referenceFrame_0_guid_sytr_clip_rect)']");
console.log(fiber_window)
var paths = fiber_window.querySelectorAll("path")

for (var p of paths){
    p_attr = (p.getAttribute("d")+"").split(" ")
    if (p_attr[2] == p_attr[5]){
        p.setAttribute("d","")
    }
    
}
"""
