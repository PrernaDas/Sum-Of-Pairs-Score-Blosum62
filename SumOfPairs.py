from Bio.SubsMat import MatrixInfo # importing Package Bio of Biopython
                                   # Package Bio has a submodule Bio.SubsMat
                                   # Bio.SubsMat has a module called MatrixInfo


seq1 = 'ATTDEWKKQRKDSHKEVERRRRENINTAINVLSDLLPVRESSKAAILACAAEYIQKLKETDEANIEKWTLQKLLSEQNASQLASANEKLQEELGNAYKEIEYMKRVLRK----------'
seq2 = 'HGSEEWHRQRRENHKEVERKRRESINTGIRELARLIPTTDTNKAQILQRAVEYIKRLKENENNNIEKWTLEKLLTEQAVSELSASNEKLKHELESAYREIEQLKRGKK-----------'
seq3 = 'TGSTAWKQQRKESHKEVERRRRQNINTAIEKLSDLLPVKETSKAAILSRAAEYIQKMKETETANIEKWTLQKLLGEQQVSSLTSANDKLEQELSKAYKNLQELKKKLKEAGIEDPTEEE'


blosum = MatrixInfo.blosum62 # MatrixInfo has a number of matrices, including blosum62
# blosum62 in Bio.SubsMat.MatrixInfo is a dictionary (biopython.org/DIST/docs/_api_161/Bio.SubsMat.MatrixInfo-module.html#blosum62
# the blosum62 dict is arranged as {('A', 'A'): 4,
#                                   ('B', 'A'): -2,
#                                   ('B', 'B'): 4,
#                                   ('B', 'C'): -3
#                                   ('B', 'D'): 4,
#                                   ......}

def cal_pairwise_score(seq1, seq2, matrix):
    """
    This function constructs pairs taking one element from seq1, and another element from the corresponding position in the second sequence.
    It then searches for the pair in the blosum62 matrix
    and assigns the value/score associted with that match pair to a score variable
    """
    score = 0
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-': # disregarding pairs in which either or one of the element has '-'
            pair = (seq1[i], seq2[i])
            if pair not in blosum: # the blosum pairs are arranged as ('B', 'A'), ('B', 'C'), so in order to look up for the value of ('A','B'), we need to reverse it ('B', 'A) 
                reverse_pair = tuple(reversed(pair))
                score += blosum[reverse_pair] # score = score + blosum[reverse_pair]
            else:
                score += blosum[pair] # score = score + blosum[pair]
    return score # gives the cumulative score for all the pairs evaluated in range of len(seq1), i.e i takes values from 0 uptil len(seq1)
            
score1_2 = cal_pairwise_score(seq1, seq2, blosum) # using the cal_pairwise_score() with seq1, seq2, and blosum62 matrix as the argumen

score1_3 = cal_pairwise_score(seq1, seq3, blosum) # score for seq1 and seq3

score2_3 = cal_pairwise_score(seq2, seq3, blosum) # score for seq2 and seq3

SP_SCORE = score1_2 + score1_3 + score2_3 # Sum of Pairs score for the multisequence alignment

print " The SP score of the given alignment is: ", SP_SCORE
