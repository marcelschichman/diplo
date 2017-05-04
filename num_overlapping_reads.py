import sys

# reads_filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq"
# positions_filename = "/home/marcel/programming/data/positions.txt"
reads_filename = "/home/marcel/programming/diplo/tests/test_dataset.fasta"
positions_filename = "/home/marcel/programming/diplo/tests/positions.blasr"

# load reads
reads = []
name2id = {}
freads = open(reads_filename, "r")
# while True:
#     l1 = freads.readline()
#     if l1 == "":
#         break
#     l2 = freads.readline()
#     l3 = freads.readline()
#     l4 = freads.readline()

#     name2id[l1.strip()[1:]] = len(reads)
#     reads.append((l1.strip()[1:], len(l2.strip())))

while True:
    l1 = freads.readline()
    if l1 == "":
        break
    if l1[0] == ">":
        name2id[l1.strip()[1:]] = len(reads)
        reads.append((l1.strip()[1:], 0))

#print(name2id["m110618_035655_42142_c100158802555500000315044108071130_s1_p0/2958/0_4673"])

# load positions
positions = {}
fpos = open(positions_filename, "r")
while True:
    l = fpos.readline()
    if l == "":
        break
    cols = l.strip().split()
    # print(cols[0])
    # print(cols[0][:cols[0].rfind("/")])
    read_name = cols[0][:cols[0].rfind("/")]
    read_id = name2id[read_name]
    if not read_id in positions:
        positions[read_id] = []

    if cols[3] == '0':
        positions[read_id].append((int(cols[6]), int(cols[7])))
    else:
        genome_length = int(cols[8])
        positions[read_id].append((genome_length - int(cols[7]), genome_length - int(cols[6])))

    
# print("reads:", len(reads))
# print("read with single position:", len(reads_positions))

def overlap_size(read1, read2):
    overlap_start = max(read1[0], read2[0])
    overlap_end = min(read1[1], read2[1])
    return max(0, overlap_end - overlap_start)

def are_overlapping(pos1, pos2):
    for p1 in pos1:
        for p2 in pos2:
            if overlap_size(p1, p2) > 200:
                return True
    return False

def get_num_overlaps(idRead):
    read_name = reads[idRead][0]
    if not idRead in positions:
        return -1, []
    ref_read_pos = positions[idRead]
    
    num_overlaps = 0
    overlapping_reads_ids = []
    for idRead2, read_pos in positions.items():
        if idRead != idRead2:
            if are_overlapping(ref_read_pos, read_pos):
                # print(idRead2)
                num_overlaps += 1
                overlapping_reads_ids.append(idRead2)
    return num_overlaps, overlapping_reads_ids

#print("num overlaps:", get_num_overlaps(1001))

for i in range(len(reads)):
    num, overlaps = get_num_overlaps(i)
    if num != -1:
        print(i, " ".join([str(x) for x in overlaps]))


x = 4