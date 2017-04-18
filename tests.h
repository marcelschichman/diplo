#pragma once
#include "matchfinder.h"


namespace Tests
{
    void AllocationOverhead();

    void KmerCountDistribution(const string& filename, int size = 13);

    long long MemoryConsumptionWithStructure(long long (*memorySize)(long long));

    void CompareDataStructures();

    void MatchFinderSimpleTest();

    void TestHashmapSize(const string& filename, int size);

    void FindOverlaps();

    void FindPath();
}