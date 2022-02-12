# nodeSeqs
`nodeSeqs.sh` is a command-line script for extracting sequences from high-degree nodes within assembly graphs, using a [GFA file](http://gfa-spec.github.io/GFA-spec/GFA1.html) as input. In the language of GFA, this script looks for *segments* that have a large number of *links*. The creation of this script was motivated by [my question on Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/questions/18476), and my approach to address this problem is modeled on [an answer by Maximilian Press](https://bioinformatics.stackexchange.com/a/18528/3967). 

## Inputs

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

## Output

This script returns a single `.fasta` file with headers formatted as follows:

    >segment[segment ID]_degree=[number of links]_[orientation of overlap]_cov=[calculated coverage]
    
Here are the first few lines of the example output [`segments.fasta`](https://github.com/acvill/nodeSeqs/blob/main/segments.fasta):

    >segment376898_degree=8_+_cov=69.91
    TGATATAATCCCCTTAGAGTAGACACATGAAATAATAAAATGGTGTAACTCTAAAGGGGGTTATTTTATGTCAAAGAAGAAATCACTTACAAGTGAAGAA
    >segment1437344_degree=7_-_cov=102.50
    GCCGATGTGGCTCAATTGGCAGAGCAGCTGATTTGTAATCAGCAGGTTATCGGTTCGAGTCCGATCATCGGCTT
    >segment8608999_degree=6_-_cov=22.66
    TTGACACACTCCCACGATTAAAATCGTGGGATTCTACTTCAACGAGGCCGCTGGCTG

## Dependencies

This script should work on Unix systems with Core Utilities and [GNU bc](https://www.gnu.org/software/bc/).

## Example

Download [`nodeSeqs.sh`](https://github.com/acvill/nodeSeqs/blob/main/nodeSeqs.sh) and [`assembly_graph_with_scaffolds.gfa.gz`](https://github.com/acvill/nodeSeqs/blob/main/assembly_graph_with_scaffolds.gfa.gz) and place them in the same directory. Decompress the GFA file and call `nodeSeqs.sh`, specifying a minimum degree of 5 links.

    cd /workdir/user/nodeSeqs
    wget https://github.com/acvill/nodeSeqs/raw/main/assembly_graph_with_scaffolds.gfa.gz
    wget https://raw.githubusercontent.com/acvill/nodeSeqs/main/nodeSeqs.sh
    gunzip assembly_graph_with_scaffolds.gfa.gz
    ./nodeSeqs.sh -g assembly_graph_with_scaffolds.gfa -d 5
    
Progress will print to the command line. Once finished, `segments.fasta` is written to the working directory. 

## Details

Link lines (beginning with `L`) are extracted from the GFA file, and segment occurrences are counted from both the `From` and `To` fields. To ensure that the final sequences capture the overlapping (linked) portions of segments, links are processed as segment-orientation pairs. Segment IDs with at least the specified minimum degree (`-d`) are written to `temp_XXXXX/segments.txt`, and the complementary segment lines (beginning with `S`) are pulled. *k*-mer counts (`KC` tag) are used to compute coverage for each segment using the approximation from [Gonnella & Kurtz 2016](https://dx.doi.org/10.7717%2Fpeerj.2681); rewritten:

![coverage estimation equation](https://user-images.githubusercontent.com/22378512/153694782-c890a32a-8863-4f64-b452-f598ca6d0447.png)

where *K*<sub>*S*</sub> is the *k*-mer count for a segment, *L*<sub>*S*</sub> is the length of a segment, and *k* is the *k*-mer size used in de Bruijn graph construction (`-k`). For (meta)SPAdes, graphs are built iteratively using an increasing *k*-mer size, so *k* is equal to the largest *k*-mer size used. Sequences are extracted corresponding to segments whose coverage is at least the specified minimum coverage (`-m`). Sequences longer than the specified context (`-c`) are truncated to the context length with respect to the link orientation. Including the `-t` flag will automatically remove the temporary folder when the program is finished, which holds intermediate files and `log.txt`. 

