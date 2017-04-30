
def load_overlaps(filename):
    """
    asdf afdsa asdf
    """
    foverlaps = open(filename, "r")
    result = {}
    while True:
        l = foverlaps.readline()
        if l == "":
            break
        ints = [int(x) for x in l.strip().split()]
        result[ints[0]] = ints[1:]
    return result

true_overlaps = load_overlaps("true_overlaps.txt")
my_overlaps = load_overlaps("my_overlaps.txt")


total_correct = 0
total_false = 0
total_missing = 0
total_true = 0


for id_read in true_overlaps:
    if id_read not in my_overlaps:
        continue

    num_correct = 0
    num_false = 0
    num_missing = 0
    num_true = len(true_overlaps[id_read])
    
    for overlap in my_overlaps[id_read]:
        if overlap in true_overlaps[id_read]:
            num_correct += 1
        else:
            num_false += 1
    num_missing = num_true - num_correct
    print("%d\t%d\t%d\t%d" % (id_read, num_correct, num_false, num_missing))

    total_correct += num_correct
    total_false += num_false
    total_missing += num_missing
    total_true += num_true

print("total correct:", total_correct)
print("total false:", total_false)
print("total missing:", total_missing)
print("total true:", total_true)

X = 5