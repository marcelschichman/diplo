
reference_fasta = open("/home/marcel/programming/data/ecoliK12_polished_assembly.fasta", "r")

genome_list = []

for line in reference_fasta:
    if line[0] == ">":
        continue
    genome_list.extend(list(line.strip()))

genome = "".join(genome_list)

read_pos_start = 3669256
read_pos_end = 3677922
read_length = read_pos_end - read_pos_start + 1

offset = 500
read_ref = genome[read_pos_start - offset:read_pos_end+1 + offset]

kmers_fasta = open("dobremiesta.fasta", "r")
kmers = []
while True:
    name = kmers_fasta.readline().strip()
    seq = kmers_fasta.readline().strip()    

    if name == "" or seq == "":
        break
    pos = int(name[8:])
    from_reference = int(name[4]) != 0
    reversed = int(name[6]) != 0
    if pos >= 0 and pos < read_length:
        kmers.append((name, pos, seq, from_reference, reversed))

def find_all(reference, substr):
    start = 0
    result = []
    while True:
        pos = reference.find(substr, start)
        if pos != -1:
            result.append(pos)
            start = pos + 1
        else:
            break
    return result

num_found = 0
num_found_multiple = 0
distances = []
distances_nfr_sense = []
distances_nfr_antisense = []

coverage = [0 for _ in range(read_length)]

for kmer in kmers:
    positions = [x - offset for x in find_all(read_ref, kmer[2])]
    if positions:
        num_found += 1
        distances.append(positions[0] - kmer[1])
        if kmer[3] == False:
            if kmer[4] == False:
                distances_nfr_sense.append(positions[0] - kmer[1])
            else:
                distances_nfr_antisense.append(positions[0] - kmer[1])
        #for i in range(positions[0], positions[0] + len(kmer[2])):
        for i in range(positions[0], positions[0] + 1):
            if i >= 0 and i < read_length:
                coverage[i] += 1
    if len(positions) > 1:
        num_found_multiple
    print (kmer[0], kmer[1], positions)

print("kmers:", len(kmers))
print("num found:", num_found)
print("num found multiple:", num_found_multiple)

distances_file = open("distances.txt", "w")
for distance in distances:
    distances_file.write(str(distance)+"\n")

distances_file = open("distances_sense.txt", "w")
for distance in distances_nfr_sense:
    distances_file.write(str(distance)+"\n")
print("average misplacement sense:", sum(distances_nfr_sense) / len(distances_nfr_sense))
print("average misplacement sense (absolute):", sum([abs(x) for x in distances_nfr_sense]) / len(distances_nfr_sense))

distances_file = open("distances_antisense.txt", "w")
for distance in distances_nfr_antisense:
    distances_file.write(str(distance)+"\n")
print("average misplacement antisense:", sum(distances_nfr_antisense) / len(distances_nfr_antisense))
print("average misplacement antisense (absolute):", sum([abs(x) for x in distances_nfr_antisense]) / len(distances_nfr_antisense))

coverage_file = open("coverage.txt", "w")
for i in range(len(coverage)):
    coverage_file.write(str(i)+"\t"+str(coverage[i])+"\n")

x = 5