# 
  
A method for creating minimum unique k-mer wiggle tracks from a set of fasta sequences, both left- and right-anchored as well as k-mer coverage and mappability wiggle tracks.
  
### Building the scripts

Before running any of the scripts, first build them using `make`. This will clone two additional git repositories that are dependencies. All scripts are build in the `bin` directory.

### Running minUniqueKmer or minUniqueBismapKmer

In order to create the minimum unique k-mer length tracks, you can run one of the provided shell scripts. For the standard minimum unique k-mer length tracks, run

```
./find_minUniqueKmer.sh <genome file (Fasta format)> [threads]
```

For minimum unique k-mer length tracks for a bisulfite treated genome, run

```
./find_minUniqueBismapKmer.sh <genome file (Fasta format)> [threads]
```

Output files will include:

- `<genome file>.refrev`  or `<genome file>.refbismap` A concatenated character file of all sequences and their reverse complements from the original fasta file (with or without bisulfite conversion)
- `<genome file>.sa` A binary suffix array file
- `<genome file>.mul.wig` A left-anchored minimum unique k-mer wiggle file
- `<genome file>.mur.wig` A right-anchored minimum unique k-mer wiggle file

### Program details

Each sequence denoted by the minimum unique value, starting from its anchored base and extending left or right for the number of bases equal to the value of the minimum unique k-mer, left-anchored (MUL) or minimum unique k-mer, right-anchored (MUR), respectively, represents the shortest unique sequence starting from the anchoring base. This considers both forward and reverse-compliment sequences for uniqueness. It does not mean that there aren't shorter unique sequences contained in a given MUL or MUR sequence, merely that conditioned on the anchoring base and direction, the MUL or MUR represents the shortest unique sequence. 

While `minUniqueKmer` calculates a value for every base position, unique k-mer values that overlap N's in the sequence or run over the bounds of chromosomes are invalid and consitute missing data in the wiggle files.

This is also one difference that produces slightly different results from a traditional k-mer counter. `minUniqueKmer` does not count unique palindromic forward-reverse complement sequences as unique. This was a conscious choice. While the sequence itself only occurs at a single position, it is direction-agnostic and therefore not truly unique.

### K-mer coverage and uniquely mappable regions

There are two additional information tracks that can be produced using the minimum unique k-mer length. The first is a k-mer coverage track. This finds the percent of a specific length k-mer that cover a given base and are unique. This is useful for finding the mappability of a region for a certain read length.

```
bin/meanKmerCoverage [-t#THREADS] <kmer length> <mul wiggle file>
```

Output is a wiggle track written to stdout. This track can be used to create a bed file containing regions that containing all bases that are mappable with at least one unique k-mer of a given length.

```
bin/kmerUniqueMapping <kmer length> <kmer coverage wiggle file>
```

Again, the output is written to stdout.

### References

The suffix array is constructed using a lightweight library from Yuta Mori and can be found <a href="https://github.com/y-256/libdivsufsort">here</a>.

The thread pool management is performed using a library from Vitaliy Vitsentiy and can be found <a href="https://github.com/vit-vit/CTPL">here</a>.

The reference concatenation and longest common prefix functions were adapted from the preprocessing steps of Sapling, a faster suffix array query method by Melanie Kirsche. Sapling can be found <a href="https://github.com/mkirsche/sapling">here</a>.
> Bioinformatics, Volume 37, Issue 6, 15 March 2021, Pages 744–749, https://doi.org/10.1093/bioinformatics/btaa911
