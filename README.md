> OVERVIEW:

The Greedy Motif Search algorithm is a heuristic approach to motif discovery, which aims to identify recurring patterns (motifs) in a set of DNA sequences. It works by iteratively selecting the best k-mers from each sequence based on a probabilistic scoring method.


> ALGORITHM WORKFLOW:

1- Initialize Motifs:
Select the first k-mer from the first DNA sequence as the initial motif.
This is an assumption that it might be part of the true motif.

2- Iterate Through Sequences

3- For each remaining sequence, use the current motifs to build a profile matrix (probability distribution of nucleotides at each position).

4- Select the Profile-most probable k-mer from the sequence using this profile.

5- Update the motif set with the new k-mer.

6- Scoring & Optimization:
Compute the score of the motifs (based on similarity to the consensus sequence).
Keep track of the best motif set found so far.

7- Return the Best Motifs found after all iterations.


> FUNCTIONS AND THEIR ROLES:

Count(motifs): Computes a count matrix for each nucleotide at each position.

Consensus(motifs): Derives the consensus string (most frequent nucleotide per position).

Score(motifs): Evaluates how well motifs match the consensus.

Profile(motifs): Generates a probability matrix from the count matrix.

Pr(text, profile): Computes the probability of a k-mer occurring based on the given profile.

ProfileMostProbableKmer(text, k, profile): Finds the highest scoring k-mer in a sequence based on the profile.

GreedyMotifSearch(DNA, k, t): Implements the Greedy Motif Search algorithm using the above functions a subroutines.


> WHY GREEDY MOTIF SEARCH?

It offers fast heuristic approach for motif discovery and is more efficient than brute force methods.


> LIMITATION:

Not guaranteed to find the optimal motif, since it relies on an initial assumption. It is thus sensitive to the first k-mer selection, which may affect final results.
