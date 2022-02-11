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

This script is a wrapper for the `PWM()` and `matchPWM()` functions from the `Biostrings` package. (see [RDocumentation](https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/matchPWM))

`PWM()` is used to convert the input position-frequency matrix into a position-weight matrix using `type = "log2probratio"`. The parameters of the Dirichlet conjugate prior are adjusted to account for the GC content of the input fasta. To save time, the number of `matchPWM()` calls are reduced to one per motif by first concatenating the input fasta into a single sequence with contigs demarcated by non-IUPAC characters. 

There are three parameters that affect the runtime of the script.

1. The length of the input fasta
2. The number of motifs represented in the position-frequency table
3. The [information content](https://en.wikipedia.org/wiki/Position_weight_matrix#Information_content) of each motif. Motifs with lower information content will have lower maximum scores, and therefore more matches above the cutoff threshold.

Future implementations should further increase speed by parallelizing `matchPWM()` calls. 


