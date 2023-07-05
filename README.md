# Oligonucleotide_Aligner
C++ &amp; Python

This tool will provide two different C++ builts to compare oligonucleotides sequence alignment with scoring system. One envision will allow you to enter the sequence on the prompt & other accept sequence files (w/wo) headers which will be aligned to the length of shorter sequence.

The python code was translate the file reader envision to python intrepeter.

Top Sequence is called seq1 and Bottom sequence is seq2
oligonucleotides one till 34
C G A G A A T G T G C C C G A G G A G G A C G G G A C C C G C T T C 
| | | | | | | | | | | | | | | | | | | | | |   | | | | | | | |   | | 
C G A G A A T G T G C C C G A G G A G G A C A G G A C C C G C C T C 
The score is 0.9411764705882353
 
 
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Top Sequence is called seq2 and Bottom sequence is seq3
oligonucleotides one till 34
C G A G A A T G T G C C C G A G G A G G A C A G G A C C C G C C T C 
| | | | | | | | | | | | |   | | | | | | | |   | | | | | | | | | | | 
C G A G A A T G T G C C C T A G G A G G A C G G G A C C C G C C T C 
The score is 0.9411764705882353

 
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Top Sequence is called seq1 and Bottom sequence is seq4
oligonucleotides one till 34
A T A A G A G A C C A C A A G C G A C C C G C A G G G C C A G A C G 
    |     |           |         | |             |     | |           
C G A G A A T G T G C C C G A G G A G G A C G G G A C C C G C T T C 
The score is 0.23529411764705882
