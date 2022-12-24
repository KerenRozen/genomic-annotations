import pandas as pd

import json
import bisect
import time


class KeyList(object):
    # bisect doesn't accept a key function before 3.10,
    # so we build the key into our sequence.
    def __init__(self, l, key):
        self.l = l
        self.key = key
    def __len__(self):
        return len(self.l)
    def __getitem__(self, index):
        return self.key(self.l[index])



raw_path = 'C:/Users/keren/PycharmProjects/project/genes.csv'
unified_path = 'C:/Users/keren/PycharmProjects/project/unified.json'
test_samples = 'C:/Users/keren/PycharmProjects/project/example_read_file.tsv'

def iter_dataset():
    with open(raw_path, 'r') as f:
        next(f)
        for row in f:
            _, chrom, start, end = row.split()
            start, end = int(start), int(end)
            yield {"chrom": chrom, "start": start, "end": end}
    # for line in test:
    #     yield line


def unify_dataset():
    lines = sorted(iter_dataset(), key=lambda line: (line['chrom'], line["start"], -line["end"]))
    cur_line = lines[0]
    chrom_to_unified_lines = {}
    for line in lines[1:]:
        if cur_line["chrom"] != line["chrom"]:
            chrom_to_unified_lines.setdefault(cur_line["chrom"], []).append(cur_line)
            cur_line = line
        elif line["start"] <= cur_line["end"]:
            if line["end"] > cur_line["end"]:
                cur_line["end"] = line["end"]
        else:
            chrom_to_unified_lines.setdefault(cur_line["chrom"], []).append(cur_line)
            cur_line = line
    chrom_to_unified_lines.setdefault(cur_line["chrom"], []).append(cur_line)
    for i in range(1, 85):
        chrom_to_unified_lines.setdefault(str(i), [])

    with open(unified_path, 'w') as f:
        json.dump(chrom_to_unified_lines, f)

f = open(unified_path, 'r')
unified_index = json.load(f)

# def test_interval_overlaps(interval):
def test_interval_overlaps(chromosome, start_pos, end_pos):
    #with open(unified_path, 'r') as f:
     #   unified_index = json.load(f)
    #if interval["chrom"] == "100":
    if chromosome == 100:
        return False
    # chrom_index = unified_index[interval["chrom"]]
    chrom_index = unified_index[str(chromosome)]
    if chrom_index == []:
        return False
    #lower_bound_index = bisect.bisect_right(chrom_index, interval["start"] + 0.5, key=lambda x: x["start"]) - 1
    #lower_bound_index = bisect.bisect_right(KeyList(chrom_index, key=lambda x: x["start"]), interval["start"] + 0.5) - 1
    lower_bound_index = bisect.bisect_right(KeyList(chrom_index, key=lambda x: x["start"]), start_pos + 0.5) - 1
    #print(lower_bound_index)
    #return interval["start"] <= chrom_index[lower_bound_index]["end"] or (lower_bound_index < len(chrom_index) - 1 and interval["end"] >= chrom_index[lower_bound_index + 1]["start"])
    return start_pos <= chrom_index[lower_bound_index]["end"] or (lower_bound_index < len(chrom_index) - 1 and end_pos >= chrom_index[lower_bound_index + 1]["start"])



def tester():
    with open(test_samples, 'r') as f:
        next(f)
        i = 0
        while i < 100000:
            row = next(f)
        #for row in f:
            _, chrom, start, end = row.split()
            chrom, start, end = int(chrom),int(start), int(end)
            print(_)
            #print(test_interval_overlaps({"chrom":chrom, "start":start, "end":end}))
            print(test_interval_overlaps(chrom, start, end))
            i += 1
    f.close()


start_time = time.time()

tester()

print("--- %s seconds ---" % (time.time() - start_time))