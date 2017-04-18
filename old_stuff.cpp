int main()
{
    string filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq";

    Tests::TestHashmapSize(filename, 13);
    return 0;


    MatchFinder mf(13);


    FASTQ fqReader(filename);
    int idSequence = 0;
    long long sumLengths = 0;
    while (fqReader)
    {
        Sequence s;
        fqReader >> s;
        mf.AddSequence(idSequence++, s);
        sumLengths += s.GetData().length();
        if (idSequence % 1000 == 0)
        {
            cout << "added: " << idSequence << endl;
        }

        if (idSequence > 20000)
            break;
    }
    //fqReader.close();
    cout << "num sequences: " << idSequence << endl;
    cout << "sum lengths: " << sumLengths << endl;
        
    map<int, vector<Match>> matches;
    mf.GetMatches(50, matches);

    int matchesCount = 0;

    for (auto matchesPair : matches)
    {
        if (matchesPair.second.size() < 10)
        {
            continue;
        }
        matchesCount += matchesPair.second.size();

        cout << "matches with sequence " << matchesPair.first << endl;
        for (auto match : matchesPair.second)
        {
            cout << match.pos1 << " " << match.pos2 << " " << (match.reversed ? "reversed" : "") << " " << match.length << endl;
        }
        cout << "extended matches with sequence " << matchesPair.first << endl;
        mf.ExtendMatches(matchesPair.second);
        for (auto match : matchesPair.second)
        {
            cout << match.pos1 << " " << match.pos2 << " " << (match.reversed ? "reversed" : "") << " " << match.length << endl;
        }
    }

    Sequence seq1;
    Sequence seq2;
    int seq1id = 50;
    int seq2id = 32;
    FASTQ fqReader2(filename);
    idSequence = 0;
    while (fqReader2)
    {
        Sequence s;
        fqReader2 >> s;

        if (idSequence == 50)
        {
            seq1 = s;
        }
        if (idSequence == 32)
        {
            seq2 = s;
        }

        idSequence++;
        if (idSequence > 100)
            break;
    }
    mf.ExtendMatches(matches[32]);
    for (auto match : matches[32])
    {
        cout << match.pos1 << " " << match.pos2 << " " << (match.reversed ? "reversed" : "") << " " << match.length << endl;
        pair<string, string> seeds = GetMatchingSequences(seq1, seq2, match);
        cout << seeds.first << endl << seeds.second << endl << endl;
        if (seeds.first != seeds.second)
        {
            cout << "megapipkos jemine" << endl;
        }
    }
    cout << "done" << endl;
}

string Reconstruction::FindPath(const Sequence& read, const vector<SequenceNode>& nodes)
{
    int endPos = read.GetData().length();

    int startNode = -1;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (nodes[i].isFromReference)
        {
            if (startNode == -1 || nodes[i].expectedPos < nodes[startNode].expectedPos)
            {
                startNode = i;
            }
        }
    }
    if (startNode == -1)
    {
        return "";
    }
    

    // [pos in read][node]
    vector<vector<Position>> positions(read.GetData().length() + 1, vector<Position>(nodes.size(), {-1, -1, 0}));

    struct ReconstructionStep
    {
        int distance;
        int posInRead;
        int node;
        int overlap;
        int previousNode;
    };

    auto compare = [](const ReconstructionStep& left, const ReconstructionStep& right)
    {
        return left.distance > right.distance || (!(right.distance > left.distance) && left.posInRead < right.posInRead);
    };

    // {edit distance, pos in read, node, overlap}
    priority_queue<ReconstructionStep, vector<ReconstructionStep>, decltype(compare)> steps(compare);

    steps.push({0, nodes[startNode].expectedPos + (int)nodes[startNode].sequence.length(), startNode, 0, -1});

    string result;

    long long iterationCounter = 0;

    while (!steps.empty())
    {
        cout << iterationCounter++ << " " << steps.size() << endl;
        auto step = steps.top();
        steps.pop();
        int &distance(step.distance);
        int &posInRead(step.posInRead);
        int &node(step.node);
        int &overlap(step.overlap);
        int &previousNode(step.previousNode);
        
        if (positions[posInRead][node].distance != -1)
        {
            continue;
        }

        positions[posInRead][node] = {distance, previousNode, overlap};

        if (posInRead >= endPos)
        {
            result = RecreatePath(positions, nodes, posInRead, node);
            return result;
        }

        for (auto& nextNode : nodes[node].overlaps)
        {
            // {alignment end pos in read, score}
            vector<pair<int, int>> alignmentScores;
            GetAlignmentScores(read, posInRead, nodes[nextNode.first], nextNode.second, alignmentScores);

            int stepLength = GetStepLength(nodes[nextNode.first], nextNode.second, posInRead);
            for (auto& alignment : alignmentScores)
            {
                int newDistance = distance + stepLength + alignment.second;
                steps.push({newDistance, alignment.first, nextNode.first, nextNode.second, node});
            }
        }
    }
    return "";
}

void Reconstruction::GetAlignmentScores(const Sequence& read, int readPos, const SequenceNode& node, int nodePos, vector<pair<int, int>>& scores)
{
    string readSeq = read.GetData().substr(readPos, min(read.GetData().length() - readPos, node.sequence.length() * 2));
    //string nodeSeq = node.sequence.substr(nodePos, node.sequence.length() - nodePos);
    string nodeSeq(node.sequence.begin() + nodePos, node.sequence.end());

    if (!alignmentMatrix.empty() && alignmentMatrix[0].size() < readSeq.length() + 1)
    {
        for (int i = 0; i < alignmentMatrix.size(); i++)
        {
            alignmentMatrix[i].resize(readSeq.length() + 1);
        }
    }

    if (alignmentMatrix.size() < nodeSeq.length() + 1)
    {
        alignmentMatrix.resize(nodeSeq.length() + 1, vector<int>(readSeq.length() + 1));
    }

    // fill trivial solutions (all inserts and all deletes)
    for (int i = 0; i <= nodeSeq.length(); i++)
    {
        alignmentMatrix[i][0] = i;
    }
    for (int i = 1; i <= readSeq.length(); i++)
    {
        alignmentMatrix[0][i] = i;
    }

    // fill alignment matrix
    for (int i = 1; i <= nodeSeq.length(); i++)
    {
        for (int j = 1; j <= readSeq.length(); j++)
        {
            int insertion = alignmentMatrix[i - 1][j] + scoring.insertion;
            int deletion = alignmentMatrix[i][j - 1] + scoring.deletion;
            int substitution = alignmentMatrix[i - 1][j - 1] + scoring.substitution * (nodeSeq[i - 1] != readSeq[j - 1]);
            alignmentMatrix[i][j] = min({insertion, deletion, substitution});
        }
    }

    scores.clear();
    for (int i = 0; i <= readSeq.length(); i++)
    {
        scores.push_back(make_pair(readPos + i, alignmentMatrix[nodeSeq.length()][i]));
    }
}

string Reconstruction::RecreatePath(const vector<vector<Position>>& positions, const vector<SequenceNode>& nodes, int posInRead, int node)
{
    string result;
    while (node >= 0)
    {
        const Position &pos(positions[posInRead][node]);
        string nextPiece = nodes[node].sequence.substr(pos.overlapWithPrevious, int(nodes[node].sequence.length()) - pos.overlapWithPrevious);
        reverse(nextPiece.begin(), nextPiece.end());
        result.append(nextPiece);
        posInRead = posInRead - nodes[node].sequence.length() + pos.overlapWithPrevious;
        node = pos.previousNode;
    }
    reverse(result.begin(), result.end());
    return result;
}

int Reconstruction::GetStepLength(const SequenceNode& node, int overlap, int posInRead)
{
    return (!node.isFromReference) * scoring.notFromReferencePenalty 
        + scoring.overlapPenalty(overlap, node.sequence.length()) 
        + scoring.misplacementPenalty(posInRead - node.expectedPos);
}

    void GetAlignmentScores();

void Tests::GetAlignmentScores()
{
    Scoring s;
    s.insertion = 1;
    s.deletion = 1;
    s.substitution = 1;
    
    Sequence read("ACTGCATGCGCTCGAGC");
    SequenceNode node1 = {"ATGACGC", 0, false, {}};

    Reconstruction r(s);
    vector<pair<int, int>> scores;
    r.GetAlignmentScores(read, 5, node1, 0, scores);
    for (auto& score : scores)
    {
        cout << score.first << ": " << score.second << endl;
    }
}