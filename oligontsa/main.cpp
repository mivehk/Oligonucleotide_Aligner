#include <iostream>
#include "osa.h"
#include <string>
using namespace std;

/*bool isSeqNameValid (string& seqname) {
    for (char c: seqname) {
        //char c;
        if (isalnum(c) || c =='_' ){
            return true;
        }
    }
    return false;
}*/


int main() {

    string name1;
    string name2;

    cout<<"Please, Enter full name of Subject sequence (Longer Fasta file): "<<endl;
    cin>> name1;
/*    if (!isSeqNameValid(name1)) {
        cerr<<"Error:  Invalid file name"<<endl;
        return 1;
    }else{
        name1 += ".fasta" ;
    }
    cout<<endl;*/
    cout<<"Please, Enter full name of Query sequence (Shorter Fasta file): "<<endl;
    cin>> name2;
/*    if (!isSeqNameValid(name2)) {
        cerr<<"Error: Invalid file name"<<endl;
        return 1;
    } else{
        name2 += ".fasta" ;
    }*/

    compare_sequences(name1,name2);

    //compare_sequences("seq2","seq1");
    //compare_sequences("seq1","seq3");
    //compare_sequences("seq2","seq3");
    //compare_sequences("seq4","seq1");
    //compare_sequences("seq1","seq4");
    //compare_sequences("seq4","seq5");
    //compare_sequences("seq5","seq4");

    return 0;

}
