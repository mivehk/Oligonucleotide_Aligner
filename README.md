*<h1> Oligonucleotide_Aligner</h1>*</br>
C++ &amp; Python</br>

Deoxyribonucleic acid (DNA) is a macromolecule polymer comprised of double-helical strands of phosphoric acid alternating with deoxyribose sugar on the backbone with appendage bonds to acidic nucleotides like base rinds of a ladder. These four nitrogenous bases form bonds between one molecule of purine (Adenine or Guanine) and one cyclohexane ring of pyrimidine (cytosine or Thymine). Therefore, the Adenine nucleotide can only form hydrogen bonds with the Thymine nucleotide, while the Guanine base can only bond with the Cytosine base. The central dogma of biology theorized that DNA in the nucleus of our cells transcribes to ribonucleic acids (RNA), then within cytoplasmic ribosomes, the messenger RNA molecules translate the chain of instructions to amino acids that are building blocks of proteins.</br> </br>
   
To understand the importance of nucleotide sequence alignment, we can talk about the inheritance of two alleles of homologous genes, one paternal and the other maternal, which primary point mutation in their nucleotides during gametic meiosis can cause genetic diseases. Moreover, Epigenetic changes (e.g., methylation of cytosine or histone modifications) can change the order of these nucleotides for different trait expressions. Consequently, post-translational alterations in the structure and functionality of the proteins can contribute to the pathogenesis of diseases like Cancer or Parkinson's. </br></br>

On the other hand, when individuals have infections, pathogens flood in the bloodstream and penetrate cells' cytoplasm with endocytosis to hijack cellular structures. Some evolved DNA viruses infiltrate the chain of nucleotides inside the nucleus to multiply their genetic code. Integration of viral DNA with the host genome can disrupt regulatory signal pathways and cause continuous cellular growth and differentiation. Consequently, scientists adopted the interspaced palindromic repeats observed in the immunity of procaryotic cells to synthesize complementary RNA sequences of viral DNA as a guided RNA that can be attached to CAS9 enzyme to regulate the expression of a gene or disrupt a viral DNA. For example, Hepatocellular carcinoma(HCC) as one of the most common types of liver cancer (Muflikhah & Santoso, 2017), has disclosed similar DNA sequences in hepatocytes that are carried in the inner core protein shell of Hepatitis-B antigen (HBcAg).</br></br>

The literature review embodies the showcases of dynamic programming algorithms for aligning nucleotide sequences. Global alignment algorithms that map the entire length of sequences for similarity by filling a scoring matrix and traceback scores-penalties to return possibilities with the highest scores. Local alignment algorithms look for similarity in sequences' specific motifs and consider a pointing system for match/mismatch nucleotides with penalties for gap opening and extension. Because BLAST, the most well-known local sequence alignment tool maintained by the National Institute of Health (NIH), is developed in C and C++, I decided to investigate the domain of pairwise sequence comparison by practical application of C++ style API design that can engage both features of global alignment and local alignment methods. Currently, the Needleman-Wunsch (NW) is one of the well-known global alignment algorithms used in software tools like GENEIOUS and CLUSTALW; in contrast, local alignment tools like BLOSUM62, BLOSUM50 and PAM250 use variations of Smith-Waterman algorithm.</br></br>

Alignment against the human genome with three billion base pairs (or other species with long strands of nucleotides) demands indexing methods like FM-index used by Burrow-Wheeler Aligner (BWA) tool and typically demands acceleration with Graphical Processing Units (GPU) or other parallel computing methods, so developed algorithm anchors identical base pair sites between pair of oligonucleotides with M and N base pairs to reduce asymptotic runtime of O(MN). I am considering combining Seed-Extend and dynamic programming techniques for a Quantitative measurement of the longest pairwise identical segments (seeds or k-mers), with their respectful locations that also engage well-known matrix scores to draw the suboptimal alignment reflected by structural variations.</br></br>



```

(Example one):

Top Sequence is called seq4 and Bottom sequence is seq5
The best Identical Site is at 246th nc which is 30.6667%
102
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
Target Sequence - Chain One: ATAAGAGACCACAAGCGACCCGCAGGGCCAGACGTTCTTCGCCGAGAGTCGTCGGGGTTTCCTGCTTCAACAGTGCTTGGACGGAACCCGGCGCTCGTTCCCCACCCCGGCCGGCCGCCCATAGCCAGCCCTCCGTCACCTCTTCACCGCACCCTCGGACTGCCCCAAGGCCCCCGCCGCCGCTCCAGCGCCGCGCAGCCACCGCCGCCGCCGCCGCCTCTCCTTAGTCGCCGCCATGACGACC
Target Sequence - Chain Two Included Anchor: CGTCCACCTCGCAGGTGCGCCAGAACTACCACCAGGACTCAGAGGCCGCCATCAACCGCCAGATCAACCTGGAGCTCTACGCCTCCTACGTTTACCTGTCCATGTCTTACTACTTTGACCGCGATGATGTGGCTTTGAAGAACTTTGCCAAATACTTTCTTCACCAATCTCATGAGGAGAGGGAACATGCTGAGAAACTGATGAAGCTGCAGAACCAACGAGGTGGCCGAATCTTCCTTCAGGATATCAAGAAACCAGACTGTGATGACTGGGAGAGCGGGCTGAATGCAATGGAGTGTGCATTACATTTGGAAAAAAATGTGAATCAGTCACTACTGGAACTGCACAAACTGGCCACTGACAAAAATGACCCCCATTTGTGTGACTTCATTGAGACACATTACCTGAATGAGCAGGTGAAAGCCATCAAAGAATTGGGTGACCACGTGACCAACTTGCGCAAGATGGGAGCGCCCGAATCTGGCTTGGCGGAATATCTCTTTGACAAGCACACCCTGGGAGACAGTGATAATGAAAGCTAAGCCTCGGGCTAATTTCCCCATAGCCGTGGGGTGACTTCCCTGGTCACCAAGGCAGTGCATGCATGTTGGGGTTTCCTTTACCTTTTCTATAAGTTGTACCAAAACATCCACTTAAGTTCTTTGATTTGTACCATTCCTTCAAATAAAGAAATTTGGTACCCAGGTGTTGTCTTTGAGGTCTTGGGATGAATCAGAAATCTATCCAGGCTATCTTCCAGATTCCTTAAGTGCCGTTGT
=============================================================================================================
Query Sequence: TTCCAAGGAACAGTGTGGCCAAGGCCTTTCGTTCCGCAATGCATGTTGGAAATAGTAGTTCTTTCCCTCCACCTCCCAACAATCCTTTTATTTACCTAAACTGGAGACCTCCATTAGGGCGGAAAGAGTGGGGTAATGGGACCTCTTCTTAAGACTGCTTTGGACACTATCTTACGCTGATATTCAGGCCTCAGGTGGCGATTCTGACCTTGGTACAGCAATTACTGTGACGTAATAAGCCGCAACTGGAAGCGTAGAGGCGAGAGGGCGGGCGCTTTACGGCGAACTCAGGTAGAATTCTTCCTTTTCCGTCTCTTTCTTTTTATGTCACCAGGGGAGGACTGGGTGGCCAACCCAGAGCCCCGAGAGATGCTAGGCTCTTTCTGTCCCGCCCTTCCTCTGACTGTGTCTTGATTTCCTATTCTGAGAGGCTATTGCTCAGCGGTTTCCGTGGCAACAGTAAAGCGTGGGAATTACAGATAAATTAAAACTGTGGAACCCCTTTCCTCGGCTGCCGCCAAGGTGTTCGGTCCTTCCGAGGAAGCTAAGGCCGCGTTGGGGTGAGACCCTCACTTCATCCGGTGAGTAGCACCGCGTCCG

======================Showing Top Candidate Alignment===================================================== 
Line 1 include Query nucleotides 1 till 35
C G T C C A C C T C G C A G G T G C G C C A G A A C T A C C A C C A G 
: : : | : | : : : : : : : : | | | : : | : | | : : | | : : | : : : : : 
T T C C A A G G A A C A G T G T G G C C A A G G C C T T T C G T T C C 
Line 2 include Query nucleotides 36 till 70
G A C T C A G A G G C C G C C A T C A A C C G C C A G A T C A A C C T 
| : : : : : : | : | : : | : : | : : | : : : | : : : : : | | : : : | : 
G C A A T G C A T G T T G G A A A T A G T A G T T C T T T C C C T C C 
Line 3 include Query nucleotides 71 till 105
G G A G C T C T A C G C C T C C T A C G T T T A C C T G T C C A T G T 
: : : : | : | : | | : : : : | : | : : : | | | | | | | : : : | : : | : 
A C C T C C C A A C A A T C C T T T T A T T T A C C T A A A C T G G A 
Line 4 include Query nucleotides 106 till 140
C T T A C T A C T T T G A C C G C G A T G A T G T G G C T T T G A A G 
: : : : : : : : | | : | : : | | : : | : | | : : : | | : | : : : : : | 
G A C C T C C A T T A G G G C G G A A A G A G T G G G G T A A T G G G 
Line 5 include Query nucleotides 141 till 175
A A C T T T G C C A A A T A C T T T C T T C A C C A A T C T C A T G A 
| : | | : | : | : : | | : | | | : : : | | : : : | | : | : | | : | : : 
A C C T C T T C T T A A G A C T G C T T T G G A C A C T A T C T T A C 
Line 6 include Query nucleotides 176 till 210
G G A G A G G G A A C A T G C T G A G A A A C T G A T G A A G C T G C 
| : : | | : : : : : : : : : | | : | | : : : : : | | | : : : | : : : : 
G C T G A T A T T C A G G C C T C A G G T G G C G A T T C T G A C C T 
Line 7 include Query nucleotides 211 till 245
A G A A C C A A C G A G G T G G C C G A A T C T T C C T T C A G G A T 
: | : : : | | : | : | : : : : : : : | | : : : : : : : : : | : | : | : 
T G G T A C A G C A A T T A C T G T G A C G T A A T A A G C C G C A A 
Line 8 include Query nucleotides 246 till 280
A T C A A G A A A C C A G A C T G T G A T G A C T G G G A G A G C G G 
: | : : | : : : : : : : : : : : | : | | : | : | : | | : : : : : : : : 
C T G G A A G C G T A G A G G C G A G A G G G C G G G C G C T T T A C 
Line 9 include Query nucleotides 281 till 315
G C T G A A T G C A A T G G A G T G T G C A T T A C A T T T G G A A A 
| : : | | | : : | | : : : : : : : : | : : : : : : : : | : : | : : : : 
G G C G A A C T C A G G T A G A A T T C T T C C T T T T C C G T C T C 
Line 10 include Query nucleotides 316 till 350
A A A A T G T G A A T C A G T C A C T A C T G G A A C T G C A C A A A 
: : : : | : | : : | | : : : : | : : : : : : : | : | | | | : : : : : : 
T T T C T T T T T A T G T C A C C A G G G G A G G A C T G G G T G G C 
Line 11 include Query nucleotides 351 till 385
C T G G C C A C T G A C A A A A A T G A C C C C C A T T T G T G T G A 
| : : : | | | : : | : | : : : | : : | | : : | : : : : : | : | : | : : 
C A A C C C A G A G C C C C G A G A G A T G C T A G G C T C T T T C T 
Line 12 include Query nucleotides 386 till 420
C T T C A T T G A G A C A C A T T A C C T G A A T G A G C A G G T G A 
: | : | : : : : : : : | : : : | : | | : : : : : : : : | : : : : : : : 
G T C C C G C C C T T C C T C T G A C T G T G T C T T G A T T T C C T 
Line 13 include Query nucleotides 421 till 455
A A G C C A T C A A A G A A T T G G G T G A C C A C G T G A C C A A C 
| : : | : : : : | : : : : | | | | : : : : : | : : : : | : : : : : : | 
A T T C T G A G A G G C T A T T G C T C A G C G G T T T C C G T G G C 
Line 14 include Query nucleotides 456 till 490
T T G C G C A A G A T G G G A G C G C C C G A A T C T G G C T T G G C 
: : : : | : | | : : : | : | : | : : : : : : | : : : : : : : | : : : : 
A A C A G T A A A G C G T G G G A A T T A C A G A T A A A T T A A A A 
Line 15 include Query nucleotides 491 till 525
G G A A T A T C T C T T T G A C A A G C A C A C C C T G G G A G A C A 
: : : : : : : : : | : : | : : | : : : : : | : : | | : : : : | | : : : 
C T G T G G A A C C C C T T T C C T C G G C T G C C G C C A A G G T G 
Line 16 include Query nucleotides 526 till 560
G T G A T A A T G A A A G C T A A G C C T C G G G C T A A T T T C C C 
: | : : : : : : : : : : | : : : | : : | | : : | | | : : : : | | : : : 
T T C G G T C C T T C C G A G G A A G C T A A G G C C G C G T T G G G 
Line 17 include Query nucleotides 561 till 595
C A T A G C C G T G G G G T G A C T T C C C T G G T C A C C A A G G C 
: : : | | : | : : : : : : | : : : | : | : : | | : : : | : | | : : | | 
G T G A G A C C C T C A C T T C A T C C G G T G A G T A G C A C C G C 
Line 18 include Query nucleotides 596 till 600
A G T G C 
: : : : : 
G T C C G 
The Identical Site is 30.6667%

====================== Needleman-Wunsch for suboptimal Alignment===================================================== 
Secondary Sequence 1: CCACCCC-GGCCGGCCGCCCATAGCCAGCCCTCCGTCACCTCTTCACCGCACCCTCGGACTGCCCCAAGGCCCCCGCCGCCGCTCCAGCGCCGCGCAGCCACCGCCGCCGCCGCCGCCTCTCCTTAGTCGCCGCCATGACGACCG-CGTCCACCTCGCAGGTGCGCCAGAACTACCACCAGGACTCAGAGGCCGCCATCAACCGCCAGATCAACCTGGAGCTCTACGCCTCCTACGTTTACCTGTCCATGTCTTACTACTTTGACCGCGATGATGTGGCTTTGAAGAACTTTGCCAAATACTTTCTTCACCAATCTCATGAGGAGAG-GGAACATGCTGAGAAACTGATGAAGCTGC-AGAACCAACGA-GGTGGC-CGAATCTTCCTTCAGGATATCA-AGAAACCAGACTGT-GATGACTGGGAGAGCGGGCTGAATGCAATGGAGTGTGCATTACATTTGGAAAAAAATGTGAATCAGTCACTACTGGAACTGCACAAACTGGCCACTGACAAAA-ATGACCCCCATTTGTGTGACTTCATTGAGACACATTACCTGAATGAGCAGGTGAAAGCCATCAAAGAATTGGGTG-ACCACGTGACCAACTTGCGCAAGATGGGAGCGCCCGAATCTGGCTTGGCGGAATATCTCTTTGACAAGCACACCCTGGGAGACAGTGATAATGAAAGCTAAGCCTCGGGCTAATTTCCCCATAGCCGTGGGGTGACTTCCCTGGTCACCAAGGCAGTG-CATGCATGTTGGGGTTTCCTTTACCTTTTCTATAAGTTGTACCAAAACATCCACTTAAGTTCTT-TGATTTGTACCATTCCTTCAAATAAAGAAATTTGGTACCCAG-GTGTTGT-CTTT-GAGGTCTTGGGATGAATCAGAAATCTATCCAGGCTATCTTCCAGATTCCTTAAGTGCCGTTGT                                                                                                     
==================================================================================== 
Secondary Sequence 2: TT--CCAAGGAAC------A-GTG-T-----GG-CC-A---AGG-C--CTTTC--G-T---TCCGC--A-----ATGC--ATG--------------------T-TG---GAA-A-----TAG---T-AG----T--T----CTTTCCCT--CCACCTCCCA-ACAATC-C----TT-T-------T-ATTT-A---C--CT-AAAC--TG-----GAG-A-----C-C-TCCA------TTAGGG-C---G----GAAA--G-AG----T---G-GGG--TAAT-G-G-G---A------CCTCT----TCT-TAAGACTGC---TTTGGA-CAC----TATCTTACGCTGATA-TTCAGGCCTCAGGTGGCGATTCTG-ACCTT---G-G-TA-CAGCAA-TTA-C-T--GTGACG--T----AA-TA---AGC-C---G---C-A-ACT-G--GAAG--CG-T-AGA--G-GCGAGA----GGGC-GG-G-CGCTTTA-CGG-C---G--AA--CTC--AGGTAGA---A---TTC-TTC-CT-TT-TC-C-GT--CT---C-TT----TC-TTTT-TATGT-CAC--C-AGGGGA-GGA---C--T--------G-GGTGGC--C--A-AC-CCA---G-AGC-CC-CGAGAGAT-GCTA-GGCTCT-TTCTGTCCCGCCC---T------TCCT--C-TGAC-TGTGTC-TT--G-----A-TT-TCCTAT----T--CTG-AGA--GGCTA---T-TGCT---CA-GCGGTT--T-C-CG---TG---GCAAC-AG-TAA-AG--CGTGGGAATTACAGA-TAAA----TTAAA-A-C-TGT-G-G-AA-C--CCCT-T-TCCTCGGCTG-C-C-GC--CAAG--GTGTT--CGGTCC-T-TC-C--GA-GG--AA-G---------C-TA-AGGCCGCGT--TGGGG-T----G--A-GA--CC---CTCA-CTT-C----A------TCC--GG------T----GAGTAGCA-CCGCGTC-C-G------
==================================================================================== 

```