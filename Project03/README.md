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
## Function to create possible motifs
```
create_possible_motif(full sequence, k-length)
```
## Function to handle fuzzy differences
```
fuzzy_diff(x,y)
{
  return(abs(x-y)) > 0.000001
}
```
## GibbsMotifFinder

```
GibbsMotifFinder(seqs, k-length) DNA - list of strings?, k-length = 10 initially
  
  # Initialize lists to store possible motifs and background (possible motifs - 1)
  poss_motif - list of k-length motifs from sequences
  background - (possible motifs - 1)
  
  
  Loop through seqs
    CALL on create_poss_motifs (based on seqs) (one random from each seqs[]) to get k-length motifs of the seqs
    APPEND to poss_motif
      
  #Set counter for j and IC
  j=counter_ic=0
       
  while loop (j==10000 or counter_ic == 100). # convergence conditions
      N = randomly index to pick one poss_motif
  
      new_seq = seqs[N]  # full sequence for randomly selected motif
  
      background = new set of seqs that doesn’t contain the randomly selected motif 
      bgPWM = Calculate PWM(PFM(new_set))
  
      Create the second for loop within 
       for in in length(new_seq-k)
         # sliding window of k-mer bp over our 50bp
             k_motif = new_seq[i:i+k]
             new_k-mer_score = score_k-mer(seq, PWM)
          
      	if (new_kmer > than old_kmer)
                  old_kmer = new_kmer
                  motif = kMotif
         pfm=build_pfm(motif)
         new_ic = pfm_ic(pfm)   (need to start the counter_ic = 0) 
   
         if (fuzzy_diff (new_ic,old_ic))
  	update old_ic = new_ic
              counter_ic = 0
          else
               counter_ic+=1
          j+=1
   returning (pfm)
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
