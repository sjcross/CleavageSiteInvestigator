<img src="./resources/logo.png">


# Features
- Run straight from command line
- Compatible with FASTA file format (.fa and .fasta)
- Determine top and bottom strand cleavage events
- Export results to .csv files
- Create visual event distributions as heatmaps and strand linkage plots


# Installation
1. Install Python (tested with Python 3.9.1)
2. Install required libraries ([BioPython](https://biopython.org/), [Seaborn](https://seaborn.pydata.org/), [SVGWrite](https://github.com/mozman/svgwrite) and [TQDM](https://github.com/tqdm/tqdm))
    - Either using Pip
    ```
    pip install biopython==1.79
    pip install tqdm==4.55.1
    pip install seaborn==0.11.1
    pip install svgwrite==1.4
    ```
    - Or using the provided Anaconda environment file ("csi.yml" in "resources" folder)
    ```
    conda env create -f csi.yml
    ```

# Usage
## Notes
- Example files for testing CSI can be downloaded from [TODO](TODO).  These files are:
  - "TODO" - Template sequence (must contain one sequence)
  - "TODO" - Cassette sequence (must contain one sequence)
  - "TODO" - Consensus sequence(s) (can contain multiple sequences)
- The above files are used throughout the following code demos.
- Each program (csi[]().py, heatmap[]().py and strandlinkageplot[]().py) can be run entirely from command line.  Full argument documentation is accessible using the `-h` (or `--help`) flag (e.g. `python csi.py -h`).


## Running CSI
CSI requires a minimum of three arguments, specifying paths to the cassette (`-c` or `--cass_path`), reference (`-r` or `--ref_path`) and test (`-t` or `--test_path`) files.
```
python csi.py 
```
2. Results should be shown almost instantly (example output below)
```

```

## Required libraries
- 
- 
- 
- 