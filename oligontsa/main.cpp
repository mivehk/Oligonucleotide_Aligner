#include <iostream>
#include "osa.h"
#include <string>
using namespace std;

int main() {

    compare_sequences("seq1","seq2");
    compare_sequences("seq1","seq3");
    compare_sequences("seq2","seq3");
    compare_sequences("seq4","seq1");
    compare_sequences("seq1","seq4");
    compare_sequences("seq4","seq5");
    compare_sequences("seq5","seq4");

    return 0;

}
