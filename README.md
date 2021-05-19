# SPPsystemCodonTranscoder
Codon transcoder to obtain ACA-less sequence for Single-Protein-Production system

# Requirement
* Biopython
* Numpy
* Pandas

# Usage
`python spp_transcoder.py [Input, fasta] [Output, fasta]`

Example:
`python spp_transcoder.py lysozyme.fasta lysozyme_acaless.fasta`

## Important
The target DNA sequence must begin without shifting of codon-reading frame.
