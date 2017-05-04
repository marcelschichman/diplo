import subprocess
import os
import time

start = time.time()
print("hello")
end = time.time()

def test_quality_and_time(tool_command, remove_overlap_graph, num_reads):
    if remove_overlap_graph:
        os.system("rm overlapgraph.txt")
        os.system(tool_command + " tests/test_dataset2.fasta --numReadsToFix 0")
    
    start = time.time()
    os.system(tool_command + " --numReadsToFix %d tests/test_dataset2.fasta  > tests/my_corrected.fasta" % (num_reads))
    end = time.time()

    blasr = subprocess.run(["blasr", "-bestn", "1", "tests/my_corrected.fasta", "/home/marcel/programming/data/ecoliK12_polished_assembly.fasta"], stdout=subprocess.PIPE)

    blasr_output = blasr.stdout.decode('ascii').strip()
    lines = blasr_output.split("\n")
    
    sum_quality = 0
    sum_length = 0
    original_sum_lengths = 29656

    for alig in lines:
        cols = alig.split()
        length = int(cols[10]) - int(cols[9])
        quality = float(cols[5])

        sum_length += length
        sum_quality += length * quality
    
    return ((end - start) / num_reads, sum_quality / sum_length, sum_length / original_sum_lengths * 100)

def test_quality_and_time_big(filename):
    blasr = subprocess.run(["blasr", "-bestn", "1", filename, "/home/marcel/programming/data/ecoliK12_polished_assembly.fasta", "-nproc", "4"], stdout=subprocess.PIPE)

    blasr_output = blasr.stdout.decode('ascii').strip()
    lines = blasr_output.split("\n")
    
    sum_quality = 0
    sum_length = 0
    original_sum_lengths = 5000000

    for alig in lines:
        cols = alig.split()
        length = int(cols[10]) - int(cols[9])
        quality = float(cols[5])

        sum_length += length
        sum_quality += length * quality
    
    return (sum_quality / sum_length, sum_length / original_sum_lengths * 100)

def test_seed_detection():
    kmerLengths = [13, 14, 15]
    minSumSeeds = [30, 40, 50, 60, 80]
    overlapFactors = [0.001, 0.005, 0.01]
    # kmerLengths = [13]
    # minSumSeeds = [40, 60]
    # overlapFactors = [0.001, 0.005]
    results = open("tests/seed_detection_results.txt", "w")

    for kmerLength in kmerLengths:
        results.write("k = %d\n" % (kmerLength))
        #results.write("")
        for overlapFactor in overlapFactors: 
            results.write("& %g" % (overlapFactor))
        results.write(" \\\\ \\hline \\hline\n")
        for minSum in minSumSeeds:
            results.write("%g" % (minSum))
            for overlapFactor in overlapFactors:
                score = test_quality_and_time("./main --overlapMinScore %d --overlapScoreRatio %f --kmerLength %d" % (minSum, overlapFactor, kmerLength), True, 10)
                results.write(" & %.2f\\%%; %.2f s; %.2f" % (score[1], score[0], score[2]))
            results.write(" \\\\ \\hline\n")
            results.flush()
            
def test_minimal_kmer_overlap():
    min_kmer_overlaps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    test_quality_and_time("./main", True, 10)

    results = open("tests/minimal_kmer_overlap_results.txt", "w")
    for overlap in min_kmer_overlaps:
        score = test_quality_and_time("./main --minKmerOverlap %d" % (overlap), False, 10)
        results.write("%d\t%.2f%%\t%.2f\t%.2f%%\n" % (overlap, score[1], score[0], score[2]))
        results.flush()
            
def test_max_edges_from_node():
    max_edges_from_node = [1, 2, 3, 5, 7, 10, 15, -1]

    test_quality_and_time("./main", True, 10)

    results = open("tests/max_edges_from_node_results.txt", "w")
    for num_edges in max_edges_from_node:
        score = test_quality_and_time("./main --maxEdgesFromNode %d" % (num_edges), False, 10)
        results.write("%d\t%.2f%%\t%.2f\t%.2f%%\n" % (num_edges, score[1], score[0], score[2]))
        results.flush()

            
def test_max_node_distance():
    node_distances = [20, 40, 60, 80, 100, 120, 140, 160]

    test_quality_and_time("./main", True, 10)

    results = open("tests/max_node_distance_results.txt", "w")
    for node_distance in node_distances:
        score = test_quality_and_time("./main --maxKmerDistance %d" % (node_distance), False, 10)
        results.write("%d\t%.2f%%\t%.2f\t%.2f%%\n" % (node_distance, score[1], score[0], score[2]))
        results.flush()

def test_max_distance_from_furthest_reaching():
    max_distances = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

    test_quality_and_time("./main", True, 10)

    results = open("tests/max_distance_from_furthest_reaching_results.txt", "w")
    for max_distance in max_distances:
        score = test_quality_and_time("./main --maxDistFromFurthest %d" % (max_distance), False, 10)
        results.write("%d\t%.2f%%\t%.2f\t%.2f%%\n" % (max_distance, score[1], score[0], score[2]))
        results.flush()

def test_local_match_rate():
    local_match_rates = [[15, 8], [15, 9], [15, 10], [20, 11], [20, 12], [20, 13], [20, 14], [30, 16], [30, 17], [30, 18]]

    test_quality_and_time("./main", True, 10)

    results = open("tests/local_match_rate.txt", "w")
    for local_match_rate in local_match_rates:
        score = test_quality_and_time("./main --minMatchesWindowSize %d --minMatches %d" % (local_match_rate[0], local_match_rate[1]), False, 10)
        results.write("%d & %d & %.2f%% & %.2f & %.2f%% \\ \hline\n" % (local_match_rate[0], local_match_rate[1], score[1], score[0], score[2]))
        results.flush()

def test_scoring():
    scorings = [[1, 2], [3, 5], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [1, 1]]

    test_quality_and_time("./main", True, 10)

    results = open("tests/scoring_results.txt", "w")
    for scoring in scorings:
        score = test_quality_and_time("./main --insertionPenalty %d --deletionPenalty %d --substitutionPenalty %d" % (scoring[0], scoring[1], scoring[1]), False, 10)
        results.write("%d & %d & %d & %.2f%% & %.2f & %.2f%% \\ \hline\n" % (scoring[0], scoring[1], scoring[1], score[1], score[0], score[2]))
        results.flush()

def test_jump_penalty():
    scorings = [[0, 0], [1/8, 0], [1/5, 0], [1/3, 0], [1, 0], [0, 1/10], [0, 1/8], [0, 1/5], [0, 1/3]]

    test_quality_and_time("./main", True, 10)

    results = open("tests/jump_penalty_results.txt", "w")
    for scoring in scorings:
        score = test_quality_and_time("./main --jumpSizeLinearPenalty %f --jumpSizeQuadraticPenalty %f" % (scoring[0], scoring[1]), False, 10)
        results.write("%f & %f & %.2f%% & %.2f & %.2f%% \\ \hline\n" % (scoring[0], scoring[1], score[1], score[0], score[2]))
        results.flush()


def test_max_kmer_misplacement():
    max_misplacements = [10, 20, 30, 40, 60, 80, 100, 120, 140]

    test_quality_and_time("./main", True, 10)

    results = open("tests/max_kmer_misplacement_results.txt", "w")
    for max_misplacement in max_misplacements:
        score = test_quality_and_time("./main --maxKmerMisplacement %d" % (max_misplacement), False, 10)
        results.write("%d\t%.2f%%\t%.2f\t%.2f%%\n" % (max_misplacement, score[1], score[0], score[2]))
        results.flush()

#test_seed_detection()
#test_minimal_kmer_overlap()
#test_max_edges_from_node()
#test_max_distance_from_furthest_reaching()
#test_max_node_distance()
#test_local_match_rate()
#test_scoring()
#test_max_kmer_misplacement()
print(test_quality_and_time_big("tests/final.fasta"))
print(test_quality_and_time_big("tests/my_corrected.fasta"))
print(test_quality_and_time_big("tests/my_corrected_2_runs.fasta"))