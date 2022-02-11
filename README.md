# nodeSeqs
nodeSeqs.sh is a command-line script for extracting sequences from high-degree nodes within assembly graphs, using a [GFA file](http://gfa-spec.github.io/GFA-spec/GFA1.html) as input. In the language of GFA, this script looks for *segments* that have a large number of *links*. The creation of this script was motivated by [my question on Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/questions/18476), and my approach in addressing this problem is modeled on [this answer by Maximilian Press](https://bioinformatics.stackexchange.com/a/18528/3967). 

### Inputs

1. GFA file (`-g`, REQUIRED)
2. minimum number of links for a segment to be considered (`-d`, REQUIRED)
3. nucleotides of context to export (`-c`, default = 100)
4. *k*-mer size used during de Bruijn graph assembly (`-k`, default = 55)
5. minimum coverage for a segment to be considered (`-m`, default = 10)

For the full list of options, call the usage statement with `./nodeSeqs.sh -h`

    Usage: ./nodeSeqs.sh -g <file> -d <integer> [-o,-c,-k,-m,-t]
      -h    print this usage and exit
      -g    (REQUIRED) GFA file
      -d    (REQUIRED) minimum number of links (degree) needed to pull a segment (node)
      -o    output path (default: pwd)
      -c    nucleotide sequence context to extract from segment ends (default: 100)
      -k    k-mer size used in assembly (default for metaSPAdes: 55)
      -m    minimum coverage filter for segments with degree >= d (default: 10)
      -t    keep all temporary files (exclude -t to delete temp folder)

### Output

This script returns a single `.fasta` file with headers formatted as follows:

    >segment[segment ID]_degree=[number of links]_[orientation of overlap]_cov=[calculated coverage]

### Dependencies

This script should work on Unix systems with GNU Core Utilities and [GNU bc](https://www.gnu.org/software/bc/).

### Details




