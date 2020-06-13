# SequenceAnalysis

## Instructions
1. Install Python (tested with Python 3.7 via [Anaconda](https://www.anaconda.com/))
2. Install required libraries from command line
```
pip install biopython
pip install seaborn
```
3. Open "main.[]()py" file in "src" folder and edit filepaths on lines 13-39 to point to local files.  There needs to be a reference sequence, a cassete sequence and one ("OTHER") or two ("SANGER") test sequence files.  The mode is set on line 41 and must be either "Seqtype.OTHER" or "Seqtype.SANGER".
4. From command line, navigate to "src" folder and run "main.[]()py"
```
python main.py
```
5. Results should be shown almost instantly (example output below)
```
Loading sequences from file
    Loading file "pRMA03+L2L2.dna"
        Reading as ".dna" format
    Loading file "CAT cassette as amplified by RA101 & 102 from pACYC184.dna"
        Reading as ".dna" format
    Loading file "45_044.ab1"
        Reading as ".ab1" format
        Found 2 instances of repeated sequence (loading first)
    Loading file "02_007.ab1"
        Reading as ".ab1" format
        Found 2 instances of repeated sequence (loading first)

Finding break in test sequence 1
    Finding cassette end in test sequence
        Best match for cassette end (TTCGCCAAgccggctttttt)
        Match score = 20.00 (quality 1.00)

    Finding cassette-adjacent sequence in reference
    Testing reference with 20 bases (1 matches found)
        Match score = 20.00 (quality 1.00)
        Reference break at position 412
        Reference break at sequence GGTAC | CCGCT

Finding break in test sequence 2
    Finding cassette end in test sequence
        Best match for cassette start RC (TTGGTGCCgccggctttttt)
        Match score = 16.00 (quality 0.80)

    Finding cassette-adjacent sequence in reference
    Reversing test target sequence
    Testing reference with 20 bases (1 matches found)
        Match score = 18.00 (quality 0.90)
        Reference break at position 408
        Reference break at sequence GCTCG | GTACC

Restriction site (3' overhang):
    5'...G GTAC↓C...3'
    3'...C↑CATG G...5'
```

## Required libraries
- [BioPython](https://biopython.org/)
