import sys
def overlap_size(read1, read2):
    overlap_start = max(read1[0], read2[0])
    overlap_end = min(read1[1], read2[1])
    return max(0, overlap_end - overlap_start)

#pos_filename = "/home/marcel/programming/data/positions.txt"
pos_filename = "/home/marcel/programming/data/positions100.txt"
#reads_filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq"
reads_filename = "/home/marcel/programming/data/pacbio.100X.fastq"
pos_begin = 700000
pos_end = 800000
max_coverage = 50

if len(sys.argv) == 6:
    pos_filename = sys.argv[1]
    reads_filename = sys.argv[2]

    pos_begin = int(sys.argv[3])
    pos_end = int(sys.argv[4])
    max_coverage = int(sys.argv[5])

fpos = open(pos_filename)
freads = open(reads_filename)

reads_to_include = {}

while True:
    pos_line = fpos.readline()
    if pos_line == "":
        break
    cols = pos_line.strip().split()

    pos_in_genome = (0, 0)
    if cols[3] == '0':
        pos_in_genome = (int(cols[6]), int(cols[7]))
    else:
        genome_length = int(cols[8])
        pos_in_genome = (genome_length - int(cols[7]), genome_length - int(cols[6]))

    if overlap_size((pos_begin, pos_end), pos_in_genome) > 1000 and float(cols[5]) > 77:
        reads_to_include[cols[0][:cols[0].rfind("/")]] = 1

print("num reads:", len(reads_to_include), file=sys.stderr)
sum_bases = 0
lengths = []
while True:
    l1 = freads.readline()
    if l1 == "":
        break
    l2 = freads.readline().strip()
    l3 = freads.readline()
    l4 = freads.readline()

    if l1.strip()[1:] in reads_to_include:
        print(">" + l1.strip()[1:])
        sum_bases += len(l2)
        lengths.append(len(l2))
        for i in range(0, len(l2), 80):
            print(l2[i:min(i + 80, len(l2))])
    if sum_bases > (pos_end - pos_begin) * max_coverage:
        break

print("sum bases:", sum_bases, file=sys.stderr)
print("average coverage:", sum_bases / (pos_end - pos_begin), file=sys.stderr)
print("average read length:", sum_bases / len(lengths), file=sys.stderr)
print("median read length:", sorted(lengths)[len(lengths) // 2], file=sys.stderr)
print("number of reads:", len(lengths), file=sys.stderr)

x = 5