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
## Driver
```
READ files 
INITIALIZE list for sequences

LOOP through files to get promoter sequences if gff_entry is a 'CDS'
  IF 'AGGAGG' is in the promoter sequence
    APPEND promoter sequence to list of sequences
    APPEND complementary sequence to list of sequences using reverse_complement(promoter sequence) 
```

## Support Functions for GibbsMotifFinder
### Function to create k-length motifs of the full length sequences
```
FUNCTION create_possible_motif(full sequence, k-length)
  N <- random integer between 0 to the sequence - k-length 
  RETURN k-length sequence starting at the Nth position

```
### Function to handle fuzzy differences between k-mer scores
```
FUNCTION fuzzy_diff(x: first value, 
                    y: second value, 
                    epsilon: tolerance between x and y)
  
  IF the absolute difference between x and y is less than epsilon
    RETURN true
  ELSE 
    RETURN false
```
### Function to normalize and log2 scores 
```
FUNCTION normalized_scores(poss_scores: list of k-mer scores)
  exp_scores <- CALCULATE the exp2 of the poss_scores to handle negative log2 values
  sum_exp_scores <- CALCULATE the sum of the exp_scores 
  probabilities <- CALCULATE the probability of each score in poss_scores by dividing the exp_scores and the sum_exp_scores
  
  RETURN probabilities
```
### Function to choose motif using probablistic selection 
```
FUNCTION choose_motif(motif_list: list of kmers from the selected sequence, 
                      weights: scores of the kmers from the selected sequence based on the background PWM)
  
  IF motif_list exists
    IF the length of the motif_list and scores are the same 
      picked_item <- SELECT random motif from list based on their scores using random.choices()
      RETURN the index of the picked_item
    ELSE
      inform user that the motif_list and scores are not the same length
  ELSE
    inform the user that there were no motifs found
```

## GibbsMotifFinder
```
GibbsMotifFinder(seqs: list of full length sequences, 
                 k: length of desired motifs)
  
  # Initialize lists to store possible motifs and background (possible motifs - 1)
  all_motifs <- list of k-length motifs from sequences
  background <- list of motifs except motifi for PWM
  pfm <- array to store PFMs
  old_ic <- comparison for IC, start with 0 then update in the loop
  j <- loop counter
  counter_ic <- counter for IC (if IC stays the same, it should update, else reset to 0)
  
  
  
  FOR each seq in seqs
    CALL on create_poss_motifs to get k-length motifs of the seq
    APPEND to all_motifs
      
       
  WHILE (j<10000 or counter_ic < 100). # convergence conditions
      N <- random index to pick one motif from all_motifs
      new_seq <- seqs[N]  # full sequence for randomly selected motif
      background <- all motifs except randomly selected at index N
      bgPWM <- Calculate PWM for the background list
      
      poss_motifs <- INITIALIZE list of all possible motifs within the new_seq 
      poss_motif_scores <- INITIALIZE list of scores for each possible motif within the new_seq
  
      CREATE sliding window for k bps of the full sequence 
        k_motif <- k-length sequence starting from index i
        APPEND k_motif to poss_motifs
        
        kmer_score <- CALCULATE score of k_motif using score_kmer() and the bgPWM
        APPEND kmer_score to poss_motif_scores
      
      weights <- CALCULATE normalized scores of poss_motif_scores using normalized_scores() to handle negative scores
      selected_motif <- CALL choose_motif(poss_motifs, weights) to get probablistically selected motif
      UPDATE index N of all_motifs using probablistically selected motif
      
      pfm <- BUILD pfm of the updated all_motif list with k values
      current_ic <- IC of the pfm
      
      IF the difference of the current_ic and old_ic are minimal
        UPDATE counter_ic by 1
      ELSE
        RESET counter_ic to 0
      
      UPDATE old_ic to current_ic
      UPDATE j by 1

   RETURN array of pfm after convergence
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
* After multiple brainstorming sessions, we were able to successfully understand the processes of Gibb's Sampling.

# Struggles
* We initially struggled to understand how to incorporate the k-mer scores and figure out how Gibb's Sampling converges.
* We had issues with tabbing in this project. Although our settings are set to 4 spaces for tabs, when programming, there were only 2 spaces, which caused a lot of indentation errors when running. 

# Personal Reflections
## Group Leader - Zoe Chow
Group leader's reflection on the project

## Other member - Jacqueline Caldwell
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
