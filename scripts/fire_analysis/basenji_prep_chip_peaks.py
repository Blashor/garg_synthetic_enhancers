import urllib.request
import subprocess
import os
import asyncio
import json
import numpy
from matplotlib import pyplot as plt


async def async_cmd(command):
    # print(command)
    process = await asyncio.create_subprocess_shell(
        command,
        # stdin=asyncio.subprocess.PIPE,
        # stdout=asyncio.subprocess.PIPE,
        # stderr=asyncio.subprocess.PIPE,
    )
    await process.wait()


async def get_bigwig(SRX, itr, factorName, genome="mm10", directory="ChIP"):
    url = "https://chip-atlas.dbcls.jp/data/" + genome + "/eachData/bw/" + SRX + ".bw"
    await async_cmd("curl " + url + " -o " + directory + "/" + SRX + "_" + str(itr) + ".bw -s")
    print(directory + "/" + SRX + "_" + str(itr) + ".bw", factorName, sep="\t")
    return SRX


def load_experiments(experimentFile):
    with open(experimentFile) as json_file:
        json_data = json.load(json_file)
        return json_data


async def main():
    # stitch_peaks("FIRE_peaks")
    # chipAtlasDownloadHelper
    # "/home/bmt26/garg/roseBed/supplementary/tf_with_controls.json"
    factors = load_experiments("../roseBed/supplementary/tf_with_controls.json")
    limit = 100000
    itr = 0
    arg_arr = []
    promise_arr = []
    for factorName in factors:
        factor_itr = 0
        factor = factors[factorName]
        for experiment in factor:
            if experiment[1] != "NA":
                factor_itr += 1
                print(factorName)
                if factorName == "Klf4":
                    arg = (factorName, experiment, itr)
                    # path = "/home/bmt26/garg/basenji/basenji-master/data/ChIP"
                    promise_bw = get_bigwig(experiment[0], itr, factorName)
                    promise_arr.append(promise_bw)
                    # get_bigwig(SRX, itr, genome="mm10")
                    itr += 1
                    # EXPERIMENT[0] is the experiment, 1 is the supposed control
                    print(arg)
                    arg_arr.append(arg)
                    if factor_itr == 5:
                        pass  # more itrs
                        # break  # do this for only one instance

        if itr == limit:
            break
    await asyncio.gather(*promise_arr)
    print(arg_arr, "\n")
    print(len(arg_arr))


#
# calculate_cutoff(numpy.random.normal(0, 1, 100))
asyncio.run(main())
