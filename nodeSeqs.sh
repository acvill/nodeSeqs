#!/bin/bash

# https://bioinformatics.stackexchange.com/questions/18476

# This is a script to annotate sequences proximal to high-degree nodes in assembly graphs
## First, the GFA1 file (from metaSPAdes) is parsed to identify segments with degree >= d
## Then, links associated with high-degree nodes are pulled, and duplicate sequences are removed
## Finally, DIAMOND aligns sequences to DNA reference database
## The results can then be compared to HUMAnN3 output to see which genes are enriched at nodes in the assembly graph

# defaults
outdir=$(pwd)
context=100
kmer=55
mincov=10

# usage statement
usage() {
  printf "Usage: $0 -g <file> -d <integer> [-o,-c,-k,-m,-t]
      -h    print this usage and exit
      -g    (REQUIRED) GFA file
      -d    (REQUIRED) minimum number of links (degree) needed to pull a segment (node)
      -o    output path (default: pwd)
      -c    nucleotide sequence context to extract from segment ends (default: 100)
      -k    k-mer size used in assembly (default for metaSPAdes: 55)
      -m    minimum coverage filter for segments with degree >= d (default: 10)
      -t    keep all temporary files (exclude -t to delete temp folder)
"
  1>&2
  exit 1
}

# get inputs
while getopts "hg:d:o:c:k:m:t" opt; do
  case ${opt} in
    g) gfa=${OPTARG};;
    d) degree=${OPTARG};;
    o) outdir=${OPTARG};;
    c) context=${OPTARG};;
    k) kmer=${OPTARG};;
    m) mincov=${OPTARG};;
    t) keeptemp=1;;
    h) usage
    exit 0;;
    *) usage
    exit 0;;
  esac
done

if [[ -z ${gfa} ]] || \
   [[ -z ${degree} ]]; then
  usage
  exit 0
else
  printf "\n######################\n"
  printf "## input GFA file   -> ${gfa}\n"
  printf "## output directory -> ${outdir}\n"
  printf "## degree cutoff    -> ${degree} links\n"
  printf "## sequence context -> ${context} nt\n"
  printf "## k-mer length     -> ${kmer}\n"
  printf "## minimum coverage -> ${mincov}\n"
  printf "######################\n\n"
fi

# make tmp directory and initialize log file
tmp_dir=$(mktemp -d -t temp_XXXXX -p ${outdir})
log=${tmp_dir}/log.txt
printf "writing intermediate files to ${tmp_dir}\n" >> ${log}

# count occurrence of segment-orientation pairs in link lines
## FROM segment
grep -P "^L\t" ${gfa} | \
  cut -f2,3 \
  > ${tmp_dir}/seg_from.txt

## TO segment
grep -P "^L\t" ${gfa} | \
  cut -f4,5 \
  > ${tmp_dir}/seg_to.txt

## merge
cat ${tmp_dir}/seg_from.txt ${tmp_dir}/seg_to.txt | \
  sort | uniq -c | sort -nr | \
  sed 's/^\s*//' | sed 's/\s/\t/g' | \
  awk -v d=${degree} '($1 >= d)' \
  > ${tmp_dir}/segments.txt

# Progress bar function from https://stackoverflow.com/a/52581824/7976890
###
PROGRESS_BAR_WIDTH=30  # progress bar length in characters

draw_progress_bar() {
  # Arguments: current value, max value, unit of measurement (optional)
  local __value=$1
  local __max=$2
  local __unit=${3:-""}  # if unit is not supplied, do not display it

  # Calculate percentage
  if (( $__max < 1 )); then __max=1; fi  # anti zero division protection
  local __percentage=$(( 100 - ($__max*100 - $__value*100) / $__max ))

  # Rescale the bar according to the progress bar width
  local __num_bar=$(( $__percentage * $PROGRESS_BAR_WIDTH / 100 ))

  # Draw progress bar
  printf "["
  for b in $(seq 1 $__num_bar); do printf "#"; done
  for s in $(seq 1 $(( $PROGRESS_BAR_WIDTH - $__num_bar ))); do printf " "; done
  printf "] $__percentage%% ($__value / $__max $__unit)\r"
}
###

# for each high-degree segment, pull sequence from segment lines
totsegs=$(wc -l < ${tmp_dir}/segments.txt)
line=0
while read segment; do

  ## separate fields for each segment line
  deg=$(echo "$segment" | awk -F$'\t' '{print $1}')
  seg=$(echo "$segment" | awk -F$'\t' '{print $2}')
  ori=$(echo "$segment" | awk -F$'\t' '{print $3}')
  sqn=$(grep -P "^S\t${seg}\t" ${gfa} | cut -f3)
  kci=$(grep -P "^S\t${seg}\t" ${gfa} | awk '{for(i=4;i<=NF;i++){if($i~/^KC/){a=$i}} print a}' | sed 's/KC:i://')
  sqL=$(printf ${sqn} | wc -c)
  cov=$(echo "scale=2; ${kci} / (${sqL} - ${kmer} + 1)" | bc)

  ## extract segment sequence given orientation and context
  ## add sequence and header to output fasta if coverage >= mincov
  if (( $(echo "${cov} < ${mincov}" | bc -l) )); then
    printf "Segment ${seg} has coverage below the set threshold (${cov} < ${mincov})\n" >> ${log}
  elif [[ ${sqL} -le ${context} ]]; then
    tsq=${sqn}
    header=">segment${seg}_degree=${deg}_${ori}_cov=${cov}"
    printf "${header}\n${tsq}\n" >> ${outdir}/segments.fasta
    printf "Segment ${seg} is <= context length (${sqL} <= ${context}) -- exporting full segment\n" >> ${log}
  elif [[ ${ori} = "+" ]]; then
    tsq=${sqn:0:${context}}
    header=">segment${seg}_degree=${deg}_${ori}_cov=${cov}"
    printf "${header}\n${tsq}\n" >> ${outdir}/segments.fasta
  elif [[ ${ori} = "-" ]]; then
    tsq=${sqn: -${context}}
    header=">segment${seg}_degree=${deg}_${ori}_cov=${cov}"
    printf "${header}\n${tsq}\n" >> ${outdir}/segments.fasta
  else
    printf "ERROR: link orientation not found for segment ${seg}, aborting\n" | tee -a ${log}
    exit 1
  fi

  ((line+=1))
  draw_progress_bar $line $totsegs "segments"

done < ${tmp_dir}/segments.txt

if [[ ${keeptemp} -ne 1 ]]; then
  rm -rf $tmp_dir
fi

printf "\nDone!\n"
