#import os
from osa import compare_sequences 

#current_width = os.get_terminal_size().columns
#print(current_width)
#os.environ['COLUMNS'] = str(1000)

compare_sequences('seq1','seq2')
compare_sequences('seq1','seq3')
compare_sequences('seq2','seq3')
compare_sequences('seq4','seq1')
compare_sequences('seq1','seq4')
compare_sequences('seq4','seq5')
compare_sequences('seq5','seq4')
#compare_sequences('macGenMD5.txt','macfileMD5.txt')
#compare_sequences('centOSfileMD5.txt','centOSGenMD5.txt')
