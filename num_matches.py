f = open("identical_kmer_count_distribution.txt", "r")

match_count = 0
for line in f:
    count, countCount = map(int, line.split())
    #print(count, " ", countCount)
    match_count += count * (count - 1) // 2 * countCount

print("total matches: ", match_count)