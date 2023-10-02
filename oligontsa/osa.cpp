
//
// Created by Kayvan Mivehnejad on 6/19/23.
//
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "osa.h"
#include <vector>
#include <map>
#include <stack>
#include "printAlign.h"
#include <queue>

using namespace std;

const string abspath="/Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/";

int count_lines(ifstream& file){
    int lc =0;
    string line;
    while(getline(file, line)){
        lc++;
    }
    return lc;
}

int findMax(int x , int y, int z){
    return max(max(x,y), z);
}

//Penalty and score values for each char used for pairwise sequences comparison
int isBPMatch (char a , char b){
    return(a==b)? 1 : -1;
}

const int MATCH_SCORE = 2;
const int MISMATCH_SCORE = -1;
const int GAP_PENALTY = -1;

int swScore(char a, char b) {
    return (a == b) ? MATCH_SCORE : MISMATCH_SCORE;
}

void compare_sequences(string sequ1 ,string sequ2 ) {

    cout << endl;
    //scrupulous words toward others has consequences! Don't!
    string firstline1;
    string firstline2;
    string lineignore1;
    string lineignore2;
    string fing1;
    string fing2;

    cout << "Top Sequence is called " << sequ1 << " and Bottom sequence is " << sequ2 << endl;
    ifstream istrm1;
    ifstream istrm2;
    //ofstream out_stream1;
    //ofstream out_stream2;
    istrm1.open(abspath + sequ1);
    if (istrm1.fail()) {
        cout << "input of first file failed. \n";
        exit(1);
    };
    istrm2.open(abspath + sequ2);
    if (istrm2.fail()) {
        cout << "input of second file failed. \n";
        exit(1);
    };

    //cout << count_lines(istrm1) <<endl;
    istrm1.clear();
    istrm1.seekg(0);
    //cout << count_lines(istrm2) <<endl;
    istrm2.clear();
    istrm2.seekg(0);

    if (getline(istrm1, firstline1, '\n')) {
        lineignore1 = firstline1;
        //cout<<lineignore1 <<endl;
    }
    if (lineignore1[0] == '>') {
        istrm1.clear();
        istrm1.seekg(0);
        getline(istrm1, firstline1, '\n');
        if (getline(istrm1, firstline1, '\n')) {
            fing1 = firstline1;
        }
    } else {
        istrm1.clear();
        istrm1.seekg(0);
        if (getline(istrm1, firstline1, '\n')) {
            fing1 = firstline1;
        }
    }

    if (getline(istrm2, firstline2, '\n')) {
        lineignore2 = firstline2;
        //cout<<lineignore2 <<endl;
    }

    if (lineignore2[0] == '>') {
        istrm1.clear();
        istrm2.seekg(0);
        getline(istrm2, firstline2, '\n');
        if (getline(istrm2, firstline2, '\n')) {
            fing2 = firstline2;
        }
    } else {
        istrm1.clear();
        istrm2.seekg(0);
        if (getline(istrm2, firstline2, '\n')) {
            fing2 = firstline2;
        }
    }

    int l1 = fing1.length(); //target
    int l2 = fing2.length(); //query
    float d1 = floor(static_cast<float>(l1) / 35);
    int r1 = l1 % 35;
    float d2 = floor(static_cast<float>(l2) / 35);
    int r2 = l2 % 35;
    int M = 1;
    int G = 0;
    int doc = 35;
    int iam = 0;


    cout<<l1 <<endl;
    cout<<l2 <<endl;
    cout<<d1 <<endl;
    cout<<d2 <<endl;
    cout<<r1 <<endl;
    cout<<r2 <<endl;


    /*
     substr hop = 0, 1, ...
     Index position = 0, 1, ...
     K count = 1, 2, ...
     (Max_ind increments only if score g8t max and it is one plus hops)
    */

    //=======================started conditional printing==============================

    string temp1 = fing1;
    string temp2 = fing2;
    //string temp1 = alignedfing1;
    //string temp2 = alignedfing2;

    //=======================================================
    //A condition for when oligonucleotides with equal length are compared for similarity. (No-SV)
    if (l1 == l2 and d1 < 1) {
        float ps1 = 0;
        float ss1 = 0;
        cout << "oligonucleotides one till " << l1 << endl;
        //cout << "yek \n";
        for (int i = 0; i < l2; i++) {
            cout << fing1[i] << " ";
        }
        cout << endl;
        for (int k = 0; k < l1; k++) {
            if (fing1[k] == fing2[k]) {
                cout << "| ";
                ps1 = ps1 + 1;
            } else {
                cout << ": ";
                ss1++;
            }
        }
        cout << endl;
        for (int j = 0; j < l2; j++) {
            cout << fing2[j] << " ";
        }
        cout << endl;
        cout << "The Identical Site is " << (static_cast<float>(ps1 / l1)) * 100 << "%" << endl;

        cout << "===============Showing Alignment Using Smith-Waterman============================================="<<endl;
        string eqswfing1 = fing1;
        string eqswfing2 = fing2;

        std::vector<std::vector<int> > eqswmatrix(l2 + 1, std::vector<int>(l1 + 1, 0));

        int eqmax_score = 0;
        int eqmax_i = 0;
        int eqmax_j = 0;

        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                int match = eqswmatrix[i - 1][j - 1] + swScore(eqswfing1[i - 1], eqswfing2[j - 1]);
                int gap_seq1 = eqswmatrix[i - 1][j] + GAP_PENALTY;
                int gap_seq2 = eqswmatrix[i][j - 1] + GAP_PENALTY;
                eqswmatrix[i][j] = findMax(match, gap_seq1, gap_seq2); //std::max({0, match, gap_seq1, gap_seq2});

                if (eqswmatrix[i][j] > eqmax_score) {
                    eqmax_score = eqswmatrix[i][j];
                    eqmax_i = i;
                    eqmax_j = j;
                }
            }
        }

        std::string eqswaligned_seq1 = "";
        std::string eqswaligned_seq2 = "";

        int eqswi = eqmax_i; //these two integer need to be swapped? think about why swi = max_j & swj = max_i
        int eqswj = eqmax_j;

        while (eqswi > 0 && eqswj > 0 && eqswmatrix[eqswi][eqswj] != 0) {
            if (eqswmatrix[eqswi][eqswj] == eqswmatrix[eqswi - 1][eqswj - 1] + swScore(eqswfing1[eqswi], eqswfing2[eqswj])) {
                eqswaligned_seq1 = eqswfing1[eqswi - 1] + eqswaligned_seq1;
                eqswaligned_seq2 = eqswfing2[eqswj - 1] + eqswaligned_seq2;
                eqswi--;
                eqswj--;
            } else if (eqswmatrix[eqswi][eqswj] == eqswmatrix[eqswi - 1][eqswj] + GAP_PENALTY) {
                eqswaligned_seq1 = eqswfing1[eqswi - 1] + eqswaligned_seq1;
                eqswaligned_seq2 = "-" + eqswaligned_seq2;
                eqswi--;
            } else {
                eqswaligned_seq1 = "-" + eqswaligned_seq1;
                eqswaligned_seq2 = eqswfing2[eqswj - 1] + eqswaligned_seq2;
                eqswj--;
            }
        }
        std::cout << "Suboptimal Subject Sequence 1: " << eqswaligned_seq1 << std::endl;
        std::cout << "Suboptimal  Query Sequence  2: " << eqswaligned_seq2 << std::endl;
    }

    if (l2 == l1 and d1 >= 1 and r1 > 0) {
        string eqopt_align;
        print_Alignment(fing1 , fing2);

        //==================Needleman-Wunsch implementation for equal length reads ==================
        //Steps: initialization, matrix filling, backtracing and score generation.
        //number of rows and colomns are plus one for when we initialize first row/column.
        vector<vector<int> > eqmatrix( l2+1, vector<int>(l1+1, 0));

        //initialized first col
        for(int i=0; i<=l2; i++){
            eqmatrix[i][0] = i * -1;
        }
        //initializing first row
        for(int j=0; j<= l1; j++){
            eqmatrix[0][j] = j * -1;
        }
        //post initialization matrix filling with insertion-score/deletion-penalty of -1 and match score one & miss penalty -1
        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                int scoreDiagonal = eqmatrix[i - 1][j - 1] + isBPMatch(fing1[i-1], fing2[j-1]);
                int scoreUp = eqmatrix[i - 1][j] - 1; //insertion
                int scoreLeft = eqmatrix[i][j - 1] - 1; //deletion
                eqmatrix[i][j] = findMax(scoreDiagonal, scoreUp, scoreLeft);
            }
        }

        //l2(alignedfing2) is query, and l1(alignedl1) is target
        int i = l2; //n
        int j = l1;  //m

        //string kmer_fing1 = fing1.substr(max4_ind4-1 ,l2);
        //string kmer_fing2 = fing2;
        //string kmer_fing3 = fing2;
//        while ( ! subopt_count4.empty()) {
//            cout << subopt_count4.top() << " kayvon "<<endl;
//            subopt_count4.pop();
//        }
        //string kmer_fing4 = fing1.substr(subopt_count4.top()-1);
        string eqnwfing1 = fing1;
       // string nwfing3 = fing1.substr(0, max4_ind4-1);
        string eqnwfing2 = fing2;
        //string nwfing4 = fing1.substr((subopt_count4.top())-1  ); //looking for index which is suboptimal location minus one

        //cout << "Target Sequence - Chain One: " << nwfing3 << endl;
        //cout << "Target Sequence - Chain Two Included Anchor: " << nwfing1 << endl;
        //cout<<"============================================================================================================= \n";
        //cout << "Query Sequence: " << nwfing2 << endl <<endl;
       // cout<<"======================Showing Top Candidate Alignment===================================================== \n";

        string eqalignedfing1 = "";
        string eqalignedfing2 = "";

        //backtracing
        while (i > 0 || j > 0) {
            //while (i > 0 and j > 0) {
            //while (i == j and  i > 0 ){
            if (i > 0 && eqmatrix[i][j] == eqmatrix[i - 1][j] - 1) {
                eqalignedfing2  = eqnwfing2[i - 1] + eqalignedfing2;
                eqalignedfing1 = "-" + eqalignedfing1;
                i--;
            } else if (j > 0 && eqmatrix[i][j] == eqmatrix[i][j - 1] - 1) {
                eqalignedfing2  = "-" + eqalignedfing2;
                eqalignedfing1 = eqnwfing1[j - 1] + eqalignedfing1;
                j--;
            } else {
                eqalignedfing1  = eqnwfing1[j - 1] + eqalignedfing1;
                eqalignedfing2 = eqnwfing2[i - 1] + eqalignedfing2;
                i--;
                j--;
            }
        }
        //}

        ////finding kmer from optimal and suboptimal alignments
        ////string kmer_fing1 = fing1.substr(max4_ind4-1 );
        //// string kmer_fing2 = fing2;
        //// string kmer_fing3 = fing2;
        //// string kmer_fing4 = fing1.substr((subopt_count4.top())-1);
        int eqkmer_opt_five;
        //int eqkmer_subopt_five;
        int eqkmer_opt_four;
        //int eqkmer_subopt_four;
        int eqcount_opt4;
        int eqcount_subopt4;
        ///=======at this point i want to show the kmer detected on primary alignment
        for(int jm=0; jm < eqnwfing2.length(); jm++){
            if(eqnwfing1[jm] == eqnwfing2[jm]){
                eqopt_align = eqopt_align + '|';
            }else{
                eqopt_align = eqopt_align + ':';
            }
        }
        eqkmer_opt_five = eqopt_align.find("|||||");
        eqkmer_opt_four = eqopt_align.find("||||");
        if ( eqkmer_opt_five != std::string::npos) {
            eqcount_opt4 = 5;
            cout<<"K-mer found on optimal alignment is: "<< eqnwfing1.substr((eqkmer_opt_five), eqcount_opt4) <<endl;
        } else if(eqkmer_opt_four != std::string::npos and eqkmer_opt_five == std::string::npos) {
            eqcount_opt4 = 4;
            cout<<"K-mer found on optimal alignment is: "<< eqnwfing1.substr((eqkmer_opt_four), eqcount_opt4) <<endl;
        } else if(eqkmer_opt_four == std::string::npos and eqkmer_opt_five == std::string::npos) {
            cout<<"No K-mer found on optimal sequence!"<<endl;
        }

        cout<<"====================== Needleman-Wunsch for suboptimal Alignment===================================================== \n";
        cout << "Suboptimal Subject Sequence: " << eqalignedfing1 << endl;
        cout<<"==================================================================================== \n";
        cout << "Suboptimal Query  Sequence : " << eqalignedfing2 << endl;
        cout<<"==================================================================================== \n";
    }

    //========================================================
    //A condition for when l1 is larger and l2 which is a read larger than 35bp (investigate structural variants e.g., INDEL)
    //Assists the highest optimal score to NW for best and secondary scores.
    map <int, float> score4;
    float max4 = std::numeric_limits<float>::min();
    int max4_ind4;
    float ps4 ;
    float ss4 ;
    vector <int> ST4 ; //vector records indices of all loci with max identities
    bool multi4;
    priority_queue<float> subopt_count4; ////record values of sorted scores
    string opt_align4;
    string subopt_align4;
;
    if(l1>l2 and d2 >= 1 and r2>0 ){
        //for (int k=0; k<(l1-1) ; k++){
        for (int k=0; k<=(l1-l2); k++){
                if (k!= 0) {temp1 = temp1.substr(1);}
                ps4 = 0;
                ss4 = 0;
                for(int h=0; h < l2 ; h++) {
                    if (temp1[h] == temp2[h]) {
                        ps4++;
                    } else {
                        ss4++;
                    }
                }
                //cout<< "Line " <<d2+1 <<" include nucleotides "<< doc+1<<" till "<<l2 << endl;
                //cout << "The score is "<<static_cast<float>(ps4/l2) <<endl;
                score4.insert(make_pair(k+1, static_cast<float>(ps4/l2)));
                //cout<<static_cast<float>(ps4/l2)<<endl;
            }
       // }
        //cout << score4.size() <<endl;
        ////subopt_count4 = 1;
        ST4.push_back(0);
//        for ( auto &pair : score4) {
//            cout<< pair.second <<" at "<< pair.first << endl;
//        }
        for ( auto &pair : score4) {
            if (pair.second > max4) {
                subopt_count4.push(max4); //stack holding indices of subs
                ST4.clear();
                max4 = pair.second;
                max4_ind4 = pair.first;
                ST4.push_back(max4_ind4); //vector holding indices of maxes
            } else if ( pair.second == max4 ){
                //subopt_count4.push(pair.first);
                multi4 = true;
                max4 = pair.second;
                ST4.push_back(pair.first);
            } else if( pair.second < max4 and pair.second > subopt_count4.top()){
                subopt_count4.push(pair.second);
            }
        }
        if (max4 != std::numeric_limits<int>::min() and multi4 != true ) {
            cout << "The best Identical Site is at "<< max4_ind4 << "th nc which hits " <<max4*100<<"%"<<endl;
        } else if (max4 != std::numeric_limits<int>::min() and multi4 == true ) {
            cout << "The best Identical Site is "<< max4*100 <<"% at multiple loci: "<<endl;
            for ( int i=0 ; i<ST4.size() ; i+=1){
                cout<<" in locus " << ST4[i]  <<endl;
            }
        } else {
            cout << "The map is empty." << std::endl;
        }
        //finding the locus from subopt value
        int locus4;
        for (auto &pair : score4){
            if (pair.second == subopt_count4.top()){
                locus4 = pair.first;
            }
        }

        cout<<"suboptimal sequence score is "<<(subopt_count4.top())*100<<"% in locus "<<locus4<<endl;
        cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n";
        //==================Needleman-Wunsch implementation
        //Steps: initialization, matrix filling, backtracing and score generation.
        //number of rows and colomns are plus one for when we initialize first row/column.
        vector<vector<int> > matrix( l2+1, vector<int>(l1+1, 0));

        //initialized first col
        for(int i=0; i<=l2; i++){
            matrix[i][0] = i * -1;
        }
        //initializing first row
        for(int j=0; j<= l1; j++){
            matrix[0][j] = j * -1;
        }
        //post initialization matrix filling with insertion-score/deletion-penalty of -1 and match score one & miss penalty -1
        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                int scoreDiagonal = matrix[i - 1][j - 1] + isBPMatch(fing1[i-1], fing2[j-1]);
                int scoreUp = matrix[i - 1][j] - 1; //insertion
                int scoreLeft = matrix[i][j - 1] - 1; //deletion
                matrix[i][j] = findMax(scoreDiagonal, scoreUp, scoreLeft);
            }
        }

        //l2(alignedfing2) is query, and l1(alignedl1) is target
        int i = l2; //n
        int j = l1;  //m

        string kmer_fing1 = fing1.substr(max4_ind4-1 ,l2);
        string kmer_fing2 = fing2;
        string kmer_fing3 = fing2;
//        while ( ! subopt_count4.empty()) {
//            cout << subopt_count4.top() << " kayvon "<<endl;
//            subopt_count4.pop();
//        }
        string kmer_fing4 = fing1.substr(subopt_count4.top()-1);
        string nwfing1 = fing1.substr(max4_ind4-1 ,l2);
        string nwfing3 = fing1.substr(0, max4_ind4-1);
        string nwfing2 = fing2;
        //string nwfing4 = fing1.substr((subopt_count4.top())-1  ); //looking for index which is suboptimal location minus one

        cout << "Target Sequence - Chain One: " << nwfing3 << endl;
        cout << "Target Sequence - Chain Two Included Anchor: " << nwfing1 << endl;
        cout<<"============================================================================================================= \n";
        cout << "Query Sequence: " << nwfing2 << endl <<endl;
        cout<<"======================Showing Top Candidate Alignment===================================================== \n";

        string alignedfing1 = "";
        string alignedfing2 = "";

        //backtracing
        while (i > 0 || j > 0) {
            //while (i > 0 and j > 0) {
            //while (i == j and  i > 0 ){
            if (i > 0 && matrix[i][j] == matrix[i - 1][j] - 1) {
                alignedfing2  = nwfing2[i - 1] + alignedfing2;
                alignedfing1 = "-" + alignedfing1;
                i--;
            } else if (j > 0 && matrix[i][j] == matrix[i][j - 1] - 1) {
                alignedfing2  = "-" + alignedfing2;
                alignedfing1 = nwfing1[j - 1] + alignedfing1;
                j--;
            } else {
                alignedfing1  = nwfing1[j - 1] + alignedfing1;
                alignedfing2 = nwfing2[i - 1] + alignedfing2;
                i--;
                j--;
            }
        }
        //}

        ////finding kmer from optimal and suboptimal alignments
        ////string kmer_fing1 = fing1.substr(max4_ind4-1 );
        //// string kmer_fing2 = fing2;
        //// string kmer_fing3 = fing2;
        //// string kmer_fing4 = fing1.substr((subopt_count4.top())-1);
        int kmer_opt_five;
        int kmer_subopt_five;
        int kmer_opt_four;
        int kmer_subopt_four;
        int count_opt4;
        int count_subopt4;
        ///=======at this point i want to show the kmer detected on primary alignment
        for(int jm=0; jm < kmer_fing2.length(); jm++){
            if(kmer_fing1[jm] == kmer_fing2[jm]){
                opt_align4 = opt_align4 + '|';
            }else{
                opt_align4 = opt_align4 + ':';
            }
        }
        kmer_opt_five = opt_align4.find("|||||");
        kmer_opt_four = opt_align4.find("||||");
        if ( kmer_opt_five != std::string::npos) {
            count_opt4 = 5;
            cout<<"K-mer found on optimal alignment is: "<< kmer_fing1.substr((kmer_opt_five), count_opt4) <<endl;
        } else if(kmer_opt_four != std::string::npos and kmer_opt_five == std::string::npos) {
                count_opt4 = 4;
                cout<<"K-mer found on optimal alignment is: "<< kmer_fing1.substr((kmer_opt_four), count_opt4) <<endl;
        } else if(kmer_opt_four == std::string::npos and kmer_opt_five == std::string::npos) {
            cout<<"No K-mer found on optimal sequence!"<<endl;
        }

        ///=======at this point i want to show the kmer detected on secondary alignment
       /* for(int km=0; km < kmer_fing3.length(); km++){
            if(kmer_fing4[km] == kmer_fing3[km]){
                subopt_align4 = subopt_align4 + '|';
            }else{
                subopt_align4 = subopt_align4 + ':';
            }
        }
        kmer_subopt_five = subopt_align4.find("|||||");
        kmer_subopt_four = subopt_align4.find("||||");
        if ( kmer_subopt_five != std::string::npos) {
            count_subopt4 = 5;
            cout<<"K-mer found on suboptimal alignment is: "<< kmer_fing4.substr((kmer_subopt_five), count_subopt4) <<endl;
        } else if(kmer_subopt_four != std::string::npos and kmer_subopt_five == std::string::npos) {
            count_subopt4 = 4;
            cout<<"K-mer found on suboptimal alignment is: "<< kmer_fing4.substr((kmer_subopt_four), count_subopt4) <<endl;
        } else if(kmer_subopt_four == std::string::npos and kmer_subopt_five == std::string::npos) {
            cout<<"No K-mer found on suboptimal sequence!"<<endl;
        }
        cout<<endl;*/
        print_Alignment(nwfing1 , nwfing2); cout<<endl;
        cout<<"====================== Needleman-Wunsch for suboptimal Alignment===================================================== \n";
        cout << "Suboptimal Subject Sequence: " << alignedfing1 << endl;
        cout<<"==================================================================================== \n";
        cout << "Suboptimal Query  Sequence : " << alignedfing2 << endl;
        cout<<"==================================================================================== \n";
    }
    cout<<endl;

    //===================================================
    //A condition for when l1 is larger and l2 is a seed smaller than 35bp(e.g., k-mers)
    map<int, float> score2;
    float max2 = std::numeric_limits<float>::min();
    int max2_ind2 = 0;
    vector<int> ST2;
    bool multi2;
    priority_queue<float> subopt_count2;
    string opt_align2;
    string suboptalign2;
    if (l1 > l2 and d2 < 1) {
        for (int k = 0; k <= (l1 - l2); k++) {
            //Hopping nucleotides in each k iteration along nucleic acid to find similarity and score.
            if (k != 0) { temp1 = temp1.substr(1); }
            float ps2 = 0;
            float ss2 = 0;
            //dropped fs-printing from here
            for (int m = 0; m < l2; m++) {
                if (temp1[m] == temp2[m]) {
                    //cout<<"| ";
                    ps2++;
                } else {
                    //cout<<"  ";
                    ss2++;
                }
            }
            //dropped ss-printing from here
            score2.insert(make_pair(k + 1, static_cast<float>(ps2 / l2)));
            //cout<<static_cast<float>(ps2/l2)<<endl;
        }
        ST2.push_back(0);
        for (const auto &pair: score2) {
            if (pair.second > max2) {
                subopt_count2.push(max2);
                ST2.clear();
                max2 = pair.second;
                max2_ind2 = pair.first;
                ST2.push_back(max2_ind2);
            } else if (pair.second == max2) {
                multi2 = true;
                max2 = pair.second;
                ST2.push_back(pair.first);
            } else if(pair.second <max2 and pair.second > subopt_count2.top()){
                subopt_count2.push(pair.second);
            }
            //cout<<pair.first<<endl;
        }
        //index 0 in seq is nc one which is at zero hop, so k one is at first bit
        if (max2 != std::numeric_limits<int>::min() and multi2 != true) {
            cout << "The best Identical Site is at " << max2_ind2 << "th nc which hits " << max2 * 100 << "%" << endl;
        } else if (max2 != std::numeric_limits<int>::min() and multi2 == true) {
            cout << "The best Identical Site is " << max2 * 100 << "% at multiple loci: " << endl;
            for (int i = 0; i < ST2.size(); i += 1) {
                cout << " in locus " << ST2[i] << endl;
            }
            //cout << "with score value of " <<max2 <<endl;
        } else {
            cout << "The map is empty." << std::endl;
        }
        //finding the locus from subopt value
        int locus2;
        for (auto &pair : score2){
            if (pair.second == subopt_count2.top()){
                locus2 = pair.first;
            }
        }
        cout<<"suboptimal sequence score is "<<(subopt_count2.top())*100<<"% in locus "<<locus2<<endl;
        cout << "===============Showing Alignment Using Smith-Waterman============================================="<<endl;
        string swfing1 = fing1.substr(max2_ind2 - 1);
        string swfing3 = fing1.substr(max2_ind2 - 1);
        string swfing2 = fing2;
        string swfing4 = fing2;
        std::vector<std::vector<int> > swmatrix(l2 + 1, std::vector<int>(l1 + 1, 0));

        int max_score = 0;
        int max_i = 0;
        int max_j = 0;

        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                int match = swmatrix[i - 1][j - 1] + swScore(swfing1[i - 1], swfing2[j - 1]);
                int gap_seq1 = swmatrix[i - 1][j] + GAP_PENALTY;
                int gap_seq2 = swmatrix[i][j - 1] + GAP_PENALTY;
                swmatrix[i][j] = findMax(match, gap_seq1, gap_seq2); //std::max({0, match, gap_seq1, gap_seq2});

                if (swmatrix[i][j] > max_score) {
                    max_score = swmatrix[i][j];
                    max_i = i;
                    max_j = j;
                }
            }
        }

        std::string swaligned_seq1 = "";
        std::string swaligned_seq2 = "";

        int swi = max_i; //these two integer need to be swapped? think about why swi = max_j & swj = max_i
        int swj = max_j;

        while (swi > 0 && swj > 0 && swmatrix[swi][swj] != 0) {
            if (swmatrix[swi][swj] == swmatrix[swi - 1][swj - 1] + swScore(swfing1[swi], swfing2[swj])) {
                swaligned_seq1 = swfing1[swi - 1] + swaligned_seq1;
                swaligned_seq2 = swfing2[swj - 1] + swaligned_seq2;
                swi--;
                swj--;
            } else if (swmatrix[swi][swj] == swmatrix[swi - 1][swj] + GAP_PENALTY) {
                swaligned_seq1 = swfing1[swi - 1] + swaligned_seq1;
                swaligned_seq2 = "-" + swaligned_seq2;
                swi--;
            } else {
                swaligned_seq1 = "-" + swaligned_seq1;
                swaligned_seq2 = swfing2[swj - 1] + swaligned_seq2;
                swj--;
            }
        }

        // }
        if (l1 > l2 and d2 < 1) {
            //fing1 = swfing1.substr(max2_ind2 - 1);
            float ps2 = 0;
            float ss2 = 0;
            cout << "Query Oligonucleotides one till " << l2 << endl;
            //cout << "do \n";
            //for (int KM: ST){
            //int xy = ST2[0];
            int g = swfing4.length();
            //cout<<x<<endl;
            //for( int M = 0 ; M < ST.size() ; M++){
            for (int i = 0; i < g; i++) {
                cout << swfing3[i] << " ";
            }
            cout << endl;
            for (int k = 0; k < g; k++) {
                if (swfing3[k] == swfing4[k]) {
                    cout << "| ";
                    ps2++;
                } else {
                    cout << ": ";
                    ss2++;
                }
            }
            cout << endl;
            for (int j = 0; j < g; j++) {
                cout << swfing4[j] << " ";
            }
            //}
            //}
            cout << endl;
            //cout << "The calculated score is "<<static_cast<float>(ps2/l2) <<endl;
            ST2.clear();
        }

        std::cout << "Suboptimal Subject Sequence 1: " << swaligned_seq1 << std::endl;
        std::cout << "Suboptimal  Query Sequence  2: " << swaligned_seq2 << std::endl;
    }

    istrm1.close();
    istrm2.close();
}