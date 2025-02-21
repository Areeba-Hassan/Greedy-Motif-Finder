def Count(motifs):
    count = {}  # initializing the count dictionary
    n = len(motifs[0])  # assuming that all motifs have the same length
    for nucl in "ACGT":  # creating an empty list with each nucleotide as the key
        count[nucl] = []
        for i in range(n):  # initiating the list of each nucl key with a zero
            count[nucl].append(0)

    l = len(motifs)
    for m in range(l):  # iterating over each motif in motifs
        for i in range(n):  # interating over each nucleotide of l motif
            nucl = motifs[m][i]  # accessing the i-th position of m-th motif in motifs
            count[nucl][i] += 1  # adding one to the value list of count where the key = the i-th index of the current motif
    return count

#creating a consensus string containing the the most popular nucl at each position of the motifs according to the scoring matrix
def Consensus(Motifs):
    l = len(Motifs[0]) #to iterate over the consensus string
    counts = Count(Motifs) #getting the score matrix
    consensus = "" #will be equal to len(Motifs[0]) given that all sequences in Motifs have the same len
    for i in range (l): #looping over every nucl in a ceratin motif of len l-- doing this first because we need to iterate column-wise now, instead of row-wise
        temp = 0
        freq = "" #will store the most popular nucl at each position
        for nucl in "ACGT": #iterating through each nucl to check it's score for a certain position i
            if counts[nucl][i] > temp:
                temp = counts[nucl][i] #store the greater score
                freq = nucl #store the nucl with the greatest score
        consensus = consensus + freq #concatenate the most popular nucl to consensus at position [i]
    return consensus

def Score(Motifs):
    con = Consensus(Motifs)
    l = len(Motifs[0])
    s = 0 #initializing score
    for motif in Motifs: #iterating each sequence in motifs
        for i in range (l): #checking each nucl of each sequence in motifs
            if motif[i] != con[i]: #comparing motifs to consensus string
                s += 1 #incrementing score in case of mismatch between motif and string
    return s

#creating a profile function to generate the profile of motif matrix storing the probabilities of count over total number of outcomes
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0]) #number of rows in motif matrix
    profile = Count(Motifs) #using count as a subroutine
    for nucl in profile.keys(): #interating over each key in the profile dictionary
        for i in range(k): #iterating over each nucleotide in the i-th value of profile dictionary
            profile[nucl][i] = (profile[nucl][i])/t #reassigning probability to the i-th element at the nucl value of profile
    return profile

#making a probability profile matrix
def Pr(text, profile):
    p = 1.0
    for i in range(len(text)):
        p *= profile[text[i]][i] #calculate probability score of a kmer based on an existing probability profile
    return p
807

def ProfileMostProbableKmer (text, k, profile): #gets the kmer from text with highest probability score based on existing probability profile
    #simply put, this function fetches the best fitting kmer to add to out motif set
    r = len(text)-k+1
    mp = ""
    max_p = -1

    for i in range(r):
        kmer = text[i:i + k]
        score = Pr(kmer, profile) #calculating probability score of a certain kmer
        if score > max_p: #comparing the probability score of current kmer to previos maximum score
            max_p = score #if the current score exceeds the previous one, reassing the max score
            mp = kmer #also reassign the best kmer in that case
    return mp

def GreedyMotifSearch(DNA, k, t):
    n = len(DNA[0])
    best_motifs = [seq[:k] for seq in DNA]

    for i in range(n - k + 1):
        motifs = [DNA[0][i:i + k]]
        for j in range(1, t):
            profile = Profile(motifs)
            motifs.append(ProfileMostProbableKmer(DNA[j], k, profile))

        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs

#test variable
DNA = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# setting t equal to the number of strings in Dna
t = len(DNA)
k=15
# Calling GreedyMotifSearch()
Motifs = GreedyMotifSearch(DNA, k, t)
print(Motifs)
print(Score(Motifs))