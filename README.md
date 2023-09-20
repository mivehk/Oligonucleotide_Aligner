*<h1> Oligonucleotide_Aligner</h1>*</br>
C++ &amp; Python</br>

Deoxyribonucleic acid (DNA) is a macromolecule polymer comprised of double-helical strands of phosphoric acid alternating with deoxyribose sugar on the backbone with appendage bonds to acidic nucleotides like base rinds of a ladder. These four nitrogenous bases form bonds between one molecule of purine (Adenine or Guanine) and one cyclohexane ring of pyrimidine (cytosine or Thymine). Therefore, the Adenine nucleotide can only form hydrogen bonds with the Thymine nucleotide, while the Guanine base can only bond with the Cytosine base. The central dogma of biology theorized that DNA in the nucleus of our cells transcribes to ribonucleic acids (RNA), then within cytoplasmic ribosomes, the messenger RNA molecules translate the chain of instructions to amino acids that are building blocks of proteins.</br> </br> 
   
To understand the importance of nucleotide sequence alignment, we can talk about the inheritance of two alleles of homologous genes, one paternal and the other maternal, which primary point mutation in their nucleotides during gametic meiosis can cause genetic diseases. Moreover, Epigenetic changes (e.g., methylation of cytosine or histone modifications) can change the order of these nucleotides for different trait expressions. Consequently, post-translational alterations in the structure and functionality of the proteins can contribute to the pathogenesis of diseases like Cancer or Parkinson's. </br></br> 

On the other hand, when individuals have infections, pathogens flood in the bloodstream and penetrate cells' cytoplasm with endocytosis to hijack cellular structures. Some evolved DNA viruses infiltrate the chain of nucleotides inside the nucleus to multiply their genetic code. Integration of viral DNA with the host genome can disrupt regulatory signal pathways and cause continuous cellular growth and differentiation. Consequently, scientists adopted the interspaced palindromic repeats observed in the immunity of procaryotic cells to synthesize complementary RNA sequences of viral DNA as a guided RNA that can be attached to CAS9 enzyme to regulate the expression of a gene or disrupt a viral DNA. For example, Hepatocellular carcinoma(HCC) as one of the most common types of liver cancer (Muflikhah & Santoso, 2017), has disclosed similar DNA sequences in hepatocytes that are carried in the inner core protein shell of Hepatitis-B antigen (HBcAg).</br></br> 

The literature review embodies the showcases of dynamic programming algorithms for aligning nucleotide sequences. Global alignment algorithms that map the entire length of sequences for similarity by filling a scoring matrix and traceback scores-penalties to return possibilities with the highest scores. Local alignment algorithms look for similarity in sequences' specific motifs and consider a pointing system for match/mismatch nucleotides with penalties for gap opening and extension. Because BLAST, the most well-known local sequence alignment tool maintained by the National Institute of Health (NIH), is developed in C and C++, I decided to investigate the domain of pairwise sequence comparison by practical application of C++ style API design that can engage both features of global alignment and local alignment methods. Currently, the Needleman-Wunsch (NW) is one of the well-known global alignment algorithms used in software tools like GENEIOUS and CLUSTALW; in contrast, local alignment tools like BLOSUM62, BLOSUM50 and PAM250 use variations of Smith-Waterman algorithm.</br></br> 

Alignment against the human genome with three billion base pairs (or other species with long strands of nucleotides) demands indexing methods like FM-index used by Burrow-Wheeler Aligner (BWA) tool and typically demands acceleration with Graphical Processing Units (GPU) or other parallel computing methods, so I developed an algorithm to align and anchor identical sites of oligonucleotides with M and N base pairs (Ezz El-Din Rashed et al., 2021, p. 109522) for enhancing the asymptotic runtime of O(MN). I propose a combination of seed-extend technique (Bayat et al., 2019, p. 1) and dynamic programming for the quantitative measurement of the longest pairwise identical segments with pinpointing their respectful locations and engaging well-known matrix scores to draw the suboptimal alignments that reflect the structural variations (Zhao et al., 2013, p. 2).</br></br> 



```

(example 1):

Please, Enter full name of Target sequence (Longer Fasta file): 
seq5
Please, Enter full name of Query sequence (Shorter Fasta file): 
seq3

Top Sequence is called seq5 and Bottom sequence is seq3
The best Identical Site is 42.6471% at multiple loci: 
 in locus 507
suboptimal tracked on 473
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
Target Sequence - Chain One: TTCCAAGGAACAGTGTGGCCAAGGCCTTTCGTTCCGCAATGCATGTTGGAAATAGTAGTTCTTTCCCTCCACCTCCCAACAATCCTTTTATTTACCTAAACTGGAGACCTCCATTAGGGCGGAAAGAGTGGGGTAATGGGACCTCTTCTTAAGACTGCTTTGGACACTATCTTACGCTGATATTCAGGCCTCAGGTGGCGATTCTGACCTTGGTACAGCAATTACTGTGACGTAATAAGCCGCAACTGGAAGCGTAGAGGCGAGAGGGCGGGCGCTTTACGGCGAACTCAGGTAGAATTCTTCCTTTTCCGTCTCTTTCTTTTTATGTCACCAGGGGAGGACTGGGTGGCCAACCCAGAGCCCCGAGAGATGCTAGGCTCTTTCTGTCCCGCCCTTCCTCTGACTGTGTCTTGATTTCCTATTCTGAGAGGCTATTGCTCAGCGGTTTCCGTGGCAACAGTAAAGCGTGGGAATTACAGATAAATTAAAACTGTGGAACCCCTTT
Target Sequence - Chain Two Included Anchor: CTCGGCTGCCGCCAAGGTGTTCGGTCCTTCCGAGGAAGCTAAGGCCGCGTTGGGGTGAGACCCTCACTTCATCCGGTGAGTAGCACCGCGTCCG
============================================================================================================= 
                             Query Sequence: CGAGAATGTGCCCTAGGAGGACGGGACCCGCCTCCGAGAATGTGCCCGAGGAGGACGGGACCCGCTTC

====================================Showing Top Candidate Alignment====================================== 
Line 1 include Query nucleotides 1 till 35
C T C G G C T G C C G C C A A G G T G T T C G G T C C T T C C G A G G 
| : : | : : | | : : : | | : | | | : | : : | | | : : | : : : | : : : : 
C G A G A A T G T G C C C T A G G A G G A C G G G A C C C G C C T C C 
Line 2 include Query nucleotides 36 till 68
A A G C T A A G G C C G C G T T G G G G T G A G A C C C T C A C T 
: | | : : : : : | | | : : : : : : | | : : | : | | | | | : | : : : 
G A G A A T G T G C C C G A G G A G G A C G G G A C C C G C T T C 
The Identical Site is 42.6471%

====================== Needleman-Wunsch for suboptimal Alignment===================================================== 
Secondary Sequence 1: ATTACAGATAAATTAAAACTGTGGAACCCCTTTCCTCGGCTGCCGCCAAGGTGTTCGGTCCTTCCGAGGAAGCTAAGGCCGCGTTGGGGTGAGACCCTCACTTCATCCGGTGAGTAGCACCGCGTCCG CCTTCCGAGGAAGCT       P       P AGTGTGGCCAAGGCCTTTCGTTCCGCAATGCATGTTGGAAATAG                 @�@��  GCG;  ���>                  P       P ---------  -----A-----------------C-- ----G---------- --A---------  -----A-----------------C-- -GTA-A-------------------C--  --A--A----A--T--C--T--T--       ------------------------C-- CTT--G--T--GTCCG G      TTAA--A--G--GTCC--      P CGGTGAGTAGCACCGCGTCCG G      TTCCAAGGAACAGTGTGGCCAAGGCCTTTCGTTCCGCAATGCATGTTGGAAATAGT
==================================================================================== 
Secondary Sequence 2: ------C-G-AG--A--ATGT----GC---------CC--TAGG---A-------G----G---------A------C---G-G------GA-C------------------C-------C---G--C-----------------CTC------C-GA---------------------------G------------------------------------------------AA----------T-----G----------------------------TG-----------C------C--------------------------------------C---G---A----------------------------G--------------G-----A----------------------------G--------G--A------CG-------------------------G----------G----A-CCC----------------------------------------G------CT--------------------------TC-----------------------------
==================================================================================== 

```
