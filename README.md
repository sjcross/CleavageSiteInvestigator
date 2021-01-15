# CSI (Cleavage Site Identifier)

## Instructions
1. Install Python (tested with Python 3.7 via [Anaconda](https://www.anaconda.com/))
2. Install Microsoft Visual C++ 14.0 (Get it with [Microsoft C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/))
3. Install required libraries from command line
```
pip install biopython
pip install seaborn
pip install svgwrite
pip install tqdm
```
4. Open "csi.[]()py" file in "src" folder and edit filepaths on lines 13-39 to point to local files.  There needs to be a reference sequence, a cassete sequence and one ("OTHER") or two ("SANGER") test sequence files.  The mode is set on line 41 and must be either "Seqtype.OTHER" or "Seqtype.SANGER".
5. From command line, navigate to "src" folder and run "csi.[]()py"
```
python csi.py
```
6. Results should be shown almost instantly (example output below)
```
INPUT: Loading sequences from file
    Loading file "pUC19.fa"
        Reading as ".fasta" format
        Loaded 1 sequence(s)
    Loading file "Chloramphenicol Cassette overhang.fa"
        Reading as ".fasta" format
        Loaded 1 sequence(s)
    Loading file "Merge_test.fasta"
        Reading as ".fasta" format
        Loaded 475 sequence(s)

PROCESSING: Sequence(s)

RESULTS:
    Position:    427, 423
    Count:       163
    Type:        5' overhang
    Sequence:    5'...T↓CTAG A...3'
                 3'...A GATC↑T...5'

    Position:    2179, 2179
    Count:       155
    Type:        Blunt end
    Sequence:    5'...AGT↓ACT...3'
                 3'...TCA↑TGA...5'

    Position:    435, 439
    Count:       36
    Type:        3' overhang
    Sequence:    5'...C TGCA↓G...3'
                 3'...G↑ACGT C...5'

    Position:    427, 424
    Count:       34
    Type:        5' overhang
    Sequence:    5'...C↓TAG A...3'
                 3'...G ATC↑T...5'

    Position:    426, 423
    Count:       8
    Type:        5' overhang
    Sequence:    5'...T↓CTA G...3'
                 3'...A GAT↑C...5'


    Local dinucleotide frequencies:

    AA: 0
    AT: 0
    AG: 73
    AC: 1
    TA: 333
    TT: 0
    TG: 0
    TC: 383
    GA: 3
    GT: 0
    GG: 0
    GC: 0
    CA: 2
    CT: 49
    CG: 1
    CC: 1


    Local 5' nucleotide frequencies:

    A: 74
    T: 716
    G: 3
    C: 53


    Local 3' nucleotide frequencies:

    A: 338
    T: 49
    G: 74
    C: 385


    Completed with 52 errors (10.947368%)
```

## Required libraries
- [BioPython](https://biopython.org/)
- [Seaborn](https://seaborn.pydata.org/)
- [TQDM](https://github.com/tqdm/tqdm)