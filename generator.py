from random import randint

def randBase():
    bases = "ACGT"
    return bases[randint(0, 3)]

def randSequence(length):
    return "".join([randBase() for x in range(length)])

def copySubstr(from_, fromBegin, fromEnd, to, toBegin):
    toList = list(to)
    toList[toBegin:(toBegin + fromEnd - fromBegin)] = list(from_)[fromBegin:fromEnd]
    return "".join(toList)

seq1 = randSequence(300)
seq2 = randSequence(300)

seq2 = copySubstr(seq1, 50, 80, seq2, 80)
seq2 = copySubstr(seq1, 100, 130, seq2, 130)
seq2 = copySubstr(seq1, 140, 180, seq2, 170)

print("@prvy")
print(seq1)
print("+")
print("...")
print("@druhy")
print(seq2)
print("+")
print("...")