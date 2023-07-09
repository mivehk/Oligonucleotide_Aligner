*<h1> Oligonucleotide_Aligner</h1>*</br>
C++ &amp; Python</br>

This tool will provide two different C++ builts to compare oligonucleotides sequence alignment with scoring system. One envision will allow you to enter the sequence on the prompt & other accept sequence files (w/wo) headers which will be aligned to the length of shorter sequence.</br>

The python code was translate the file reader envision to python intrepeter.</br>
```
Top Sequence is called seq1 and Bottom sequence is seq2
oligonucleotides one till 34
C G A G A A T G T G C C C G A G G A G G A C G G G A C C C G C T T C
| | | | | | | | | | | | | | | | | | | | | |   | | | | | | | |   | |
C G A G A A T G T G C C C G A G G A G G A C A G G A C C C G C C T C 
The score is 0.9411764705882353
```
 
==================================================</br>
```
Top Sequence is called seq2 and Bottom sequence is seq3
oligonucleotides one till 34
C G A G A A T G T G C C C G A G G A G G A C A G G A C C C G C C T C 
| | | | | | | | | | | | |   | | | | | | | |   | | | | | | | | | | | 
C G A G A A T G T G C C C T A G G A G G A C G G G A C C C G C C T C 
The score is 0.9411764705882353
```
 
======================================================</br>
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Top Sequence is called seq5 and Bottom sequence is seq4
Line 1 include nucleotides 1 till 35
T T C C A A G G A A C A G T G T G G C C A A G G C C T T T C G T T C C 
  |       | |               |   |   | |                     |         
A T A A G A G A C C A C A A G C G A C C C G C A G G G C C A G A C G T 

Line 2 include nucleotides 36 till 70
G C A A T G C A T G T T G G A A A T A G T A G T T C T T T C C C T C C 
  |       | |           |             |       | | |   |   |           
T C T T C G C C G A G A G T C G T C G G G G T T T C C T G C T T C A A 

Line 3 include nucleotides 71 till 105
A C C T C C C A A C A A T C C T T T T A T T T A C C T A A A C T G G A 
      |   |         |                             |         |         
C A G T G C T T G G A C G G A A C C C G G C G C T C G T T C C C C A C 

Line 4 include nucleotides 106 till 140
G A C C T C C A T T A G G G C G G A A A G A G T G G G G T A A T G G G 
    |     | |         |     |     |       | |                         
C C C G G C C G G C C G C C C A T A G C C A G C C C T C C G T C A C C 

Line 5 include nucleotides 141 till 175
A C C T C T T C T T A A G A C T G C T T T G G A C A C T A T C T T A C 
  |   | |     |     |           |       | |     |           |       | 
T C T T C A C C G C A C C C T C G G A C T G C C C C A A G G C C C C C 

Line 6 include nucleotides 176 till 210
G C T G A T A T T C A G G C C T C A G G T G G C G A T T C T G A C C T 
| |   |         | |     | |     |     |     | |   |             | |   
G C C G C C G C T C C A G C G C C G C G C A G C C A C C G C C G C C G 

Line 7 include nucleotides 211 till 245
T G G T A C A G C A A T T A C T G T G A C G T A A T A A G C C G C A A 
    |                           | |     |                 |     |     
C C G C C G C C T C T C C T T A G T C G C C G C C A T G A C G A C C G 

Line 8 include nucleotides 246 till 280
C T G G A A G C G T A G A G G C G A G A G G G C G G G C G C T T T A C 
|         |   |         | | |   |   |       |             |       |   
C G T C C A C C T C G C A G G T G C G C C A G A A C T A C C A C C A G 

Line 9 include nucleotides 281 till 315
G G C G A A C T C A G G T A G A A T T C T T C C T T T T C C G T C T C 
|   |     |                   |               |           |     |     
G A C T C A G A G G C C G C C A T C A A C C G C C A G A T C A A C C T 

Line 10 include nucleotides 316 till 350
T T T C T T T T T A T G T C A C C A G G G G A G G A C T G G G T G G C 
          |   |               |   |   |                           |   
G G A G C T C T A C G C C T C C T A C G T T T A C C T G T C C A T G T 

Line 11 include nucleotides 351 till 385
C A A C C C A G A G C C C C G A G A G A T G C T A G G C T C T T T C T 
|       |   |             |                       | | | |   |         
C T T A C T A C T T T G A C C G C G A T G A T G T G G C T T T G A A G 

Line 12 include nucleotides 386 till 420
G T C C C G C C C T T C C T C T G A C T G T G T C T T G A T T T C C T 
    |         | |           | |     | |         |         |           
A A C T T T G C C A A A T A C T T T C T T C A C C A A T C T C A T G A 

Line 13 include nucleotides 421 till 455
A T T C T G A G A G G C T A T T G C T C A G C G G T T T C C G T G G C 
          |   | |       |     | |       |   |   |   |       |     | | 
G G A G A G G G A A C A T G C T G A G A A A C T G A T G A A G C T G C 

Line 14 include nucleotides 456 till 490
A A C A G T A A A G C G T G G G A A T T A C A G A T A A A T T A A A A 
|     |     | |   |   |     | |         |                         |   
A G A A C C A A C G A G G T G G C C G A A T C T T C C T T C A G G A T 

Line 15 include nucleotides 491 till 525
C T G T G G A A C C C C T T T C C T C G G C T G C C G C C A A G G T G 
  |       | | |   | |             |                 |       | |     | 
A T C A A G A A A C C A G A C T G T G A T G A C T G G G A G A G C G G 

Line 16 include nucleotides 526 till 560
T T C G G T C C T T C C G A G G A A G C T A A G G C C G C G T T G G G 
      |                 |     |           |       |                   
G C T G A A T G C A A T G G A G T G T G C A T T A C A T T T G G A A A 

Line 17 include nucleotides 561 till 595
G T G A G A C C C T C A C T T C A T C C G G T G A G T A G C A C C G C 
      |                     | | |             | |       | | | |       
A A A A T G T G A A T C A G T C A C T A C T G G A A C T G C A C A A A 

Line 18 include nucleotides 596 till 600
G T C C G 
  |       
C T G G C 
The score is 0.25666666666666665

```
