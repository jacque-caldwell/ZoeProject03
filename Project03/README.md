# Introduction
Gibbs sampling is an MCMC approach to identify enrichments. Here, we will implement a method to identify motifs from a set of regions.

Probability in this case will be robability of choosing position in DNA _m_ in (A_sample/sum(A_sub_l) for possitions l in DNA_sub_i

## Notes and Important Considerations:

Note: I have also added a function to motif_ops.py This will calculate
the information content of your motifs. This is useful to observe the
progression of your Gibbs sampler as well as a measure of convergence.
You can use this function as IC = pfm_ic(pfm). You should expect a
slow increase of IC until it plateaus, such as in the plot below from
your lecture slides:

We will need to score each sequence with a PWM using the score_kmer() or
score_sequence() functions.  You will need to investigate the help
documentation and libraries to identify how best to use these functions.
These sites are often not strand-specific, and so both scores on the
negative as well as positive strands should be considered.

To select a random sequence, use random.randint() or
numpy.random.randint()

To select a new position $m$ (as defined below) use random.choices() or
numpy.random.choice()

# Pseudocode
Put pseudocode in this box:

From diagram shown in class 9/24/2025:

```
Set of unaligned input sequences
   read of the file --> list of sequences that have promoter (AGGAGG) in them;
   and any promoter seq of GFF file with type "CDS"; these will be 50 bp long

  Randomly assign initial motif positions in each sequence


  iterate:
    Choose Sequence at Random

    Create Motif and background model all other sequences (create PWM(PFM(sequ, score)))

    Score all possible motifs in chosen sequence PWM.  loop over "40" looking at 10bp -> loop over score_kmer (PWM?). save somewhere....   want motif with the "highest frequency"
```

# Functions that we were given:
##data_readers.py
* Class GffEntry:
    __init__(self, args)
    __str__(self)
    __len__(self)
    __eq__(self, other)
    __lt__(self, other)

* def get_gff((str)gff_file)
* def get_fasta((str)file)

## motif_ops.py
* def build_pfm(sequences: List[str], length: int) -> np.ndarray:
* def build_pwm(pfm: np.ndarray) -> np.ndarray:
* def score_kmer(seq: str, pwm: np.ndarray) -> float:
* def pfm_ic(pfm: np.ndarray) -> float:

##seq_ops.py
* def reverse_complement(seq):
* def get_seq(seq, start, end, strand, size):
        size (int): how far upstream to get extra sequence (default: 50)
# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
