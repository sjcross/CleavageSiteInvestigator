# SequenceAnalysis

## Instructions
1. Install Python (tested with Python 3.7 via [Anaconda](https://www.anaconda.com/))
2. Install required libraries from command line
```
pip install biopython
```
3. Open "main.[]()py" file in "src" folder and edit filepaths on lines 9-16 to point to local files
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
    Loading file "48_041.ab1"
        Reading as ".ab1" format
        Found 2 instances of repeated sequence (loading first)

Finding cassette end in test sequence
    Best match for cassette start RC (TTGGTGCCgccggctttttt)
    Match score = 20.00 (quality 1.00)

Finding cassette-adjacent sequence in reference
    Testing reference with 20 bases (1 matches found)
    Match score = 20.00 (quality 1.00)
    Reference break at position 2396
    Reference break at sequence AcTCC | TCGAG
```

## Required libraries
- [BioPython](https://biopython.org/)