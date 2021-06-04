# 
  
A method for creating minimum unique k-mer wiggle tracks from a set of fasta sequences, both left- and right-anchored.
  
### Running minUniqueKmer:

The first time minUniqueKmer is run, the needed programs will be built if they do not currently exist.

```
./find_minUniqueKmer.sh <genome file (Fasta format)>
```

Output files will include:

- `<genome file>.refrev` A concatenated character file of all sequences and their reverse complements from the original fasta file
- `<genome file>.sa` A binary suffix array file
- `<genome file>.mul.wig` A left-anchored minimum unique k-mer wiggle file
- `<genome file>.mur.wig` A right-anchored minimum unique k-mer wiggle file
  
### Program details

Each sequence denoted by the minimum unique value, starting from its anchored base and extending left or right for the number of bases equal to the value of the minimum unique k-mer, left-anchored (MUL) or minimum unique k-mer, right-anchored (MUR), respectively, represents the shortest unique sequence starting from the anchoring base. This considers both forward and reverse-compliment sequences for uniqueness. It does not mean that there aren't shorter unique sequences contained in a given MUL or MUR sequence, merely that conditioned on the anchoring base and direction, the MUL or MUR represents the shortest unique sequence. 

While `minUniqueKmer` calculates a value for every base position, unique k-mer values that overlap N's in the sequence or run over the bounds of chromosomes are invalid and consitute missing data in the wiggle files.

This is also one difference that produces slightly different results from a traditional k-mer counter. `minUniqueKmer` does not count unique palindromic forward-reverse complement sequences as unique. This was a conscious choice. While the sequence itself only occurs at a single position, it is direction-agnostic and therefore not truly unique.

### References

The suffix array is constructed using a lightweight library from Yuta Mori and can be found <a href="https://github.com/y-256/libdivsufsort">here</a>.

The reference concatenation and longest common prefix functions were adapted from the preprocessing steps of Sapling, a faster suffix array query method by Melanie Kirsche. Sapling can be found <a href="https://github.com/mkirsche/sapling">here</a>.
> Bioinformatics, Volume 37, Issue 6, 15 March 2021, Pages 744–749, https://doi.org/10.1093/bioinformatics/btaa911
