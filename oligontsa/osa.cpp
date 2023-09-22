
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

    /*
    cout<<l1 <<endl;
    cout<<l2 <<endl;
    cout<<d1 <<endl;
    cout<<d2 <<endl;
    cout<<r1 <<endl;
    cout<<r2 <<endl;
    */

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
    //first and sixth condition for when oligonucleotides with equal length are compared for similarity. (No-SV)
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
    }
    if (l2 == l1 and d1 >= 1 and r1 > 0) {
        float ps6 = 0;
        float ss6 = 0;
        //cout<<"shish \n";
        for (int i = 0; i < d1; i++) {
            doc = 35 * M;
            iam = 35 * G;
            cout << "Line " << i + 1 << " include Query nucleotides " << iam + 1 << " till " << doc << endl;
            for (int j = iam; j < doc; j++) {
                cout << fing1[j] << " ";
            }
            cout << endl;
            for (int k = iam; k < doc; k++) {
                if (fing1[k] == fing2[k]) {
                    cout << "| ";
                    ps6++;
                } else {
                    cout << ": ";
                    ss6++;
                }
            }
            cout << endl;
            for (int a = iam; a < doc; a++) {
                cout << fing2[a] << " ";
            }
            cout << endl;
            M = M + 1;
            G = G + 1;
        }
        cout << "Line " << d1 + 1 << " include Query nucleotides " << doc + 1 << " till " << l2 << endl;
        for (int i = doc; i < l1; i++) {
            cout << fing1[i] << " ";
        }
        cout << endl;
        for (int i = doc; i < l1; i++) {
            if (fing1[i] == fing2[i]) {
                cout << "| ";
                ps6++;
            } else {
                cout << "  ";
                ss6++;
            }
        }
        cout << endl;
        for (int i = doc; i < l1; i++) {
            cout << fing2[i] << " ";
        }
        cout << endl;
        cout << "The Identical Site is " << (static_cast<float>(ps6 / l2)) * 100 << "%" << endl;
    }
    //===================================================
    //second condition for when l1 is larger and l2 is a seed smaller than 35bp(e.g., k-mers)
    map<int, float> score2;
    float max2 = std::numeric_limits<float>::min();
    int max2_ind2 = 0;
    vector<int> ST2;
    bool multi2;
    if (l1 > l2 and d2 < 1) {
        for (int k = 0; k < (l1 - 1); k++) {
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
        for (const auto &pair: score2) {
            if (pair.second > max2) {
                ST2.clear();
                max2 = pair.second;
                max2_ind2 = pair.first;
                ST2.push_back(max2_ind2);
            } else if (pair.second == max2) {
                multi2 = true;
                max2 = pair.second;
                ST2.push_back(pair.first);
            }
            //cout<<pair.first<<endl;
        }
        //index 0 in seq is nc one which is at zero hop, so k one is at first bit
        if (max2 != std::numeric_limits<int>::min() and multi2 != true) {
            cout << "The best Identical Site is at " << max2_ind2 << "th nc which is " << max2 * 100 << "%" << endl;
        } else if (max2 != std::numeric_limits<int>::min() and multi2 == true) {
            cout << "The best Identical Site is " << max2 * 100 << "% at multiple loci: " << endl;
            for (int i = 0; i < ST2.size(); i += 1) {
                cout << " in locus " << ST2[i] << endl;
            }
            //cout << "with score value of " <<max2 <<endl;
        } else {
            cout << "The map is empty." << std::endl;
        }

        cout << "===============Showing Alignment Using Smith-Waterman==============================================="
             << endl;

        string swfing1 = fing1.substr(max2_ind2 - 1);
        string swfing3 = fing1.substr(max2_ind2 - 1);
        string swfing2 = fing2;
        string swfing4 = fing2;
        std::vector<std::vector<int>> swmatrix(l2 + 1, std::vector<int>(l1 + 1, 0));

        int max_score = 0;
        int max_i = 0;
        int max_j = 0;

        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                int match = swmatrix[i - 1][j - 1] + swScore(swfing1[i - 1], swfing2[j - 1]);
                int gap_seq1 = swmatrix[i - 1][j] + GAP_PENALTY;
                int gap_seq2 = swmatrix[i][j - 1] + GAP_PENALTY;
                swmatrix[i][j] = std::max({0, match, gap_seq1, gap_seq2});

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
        std::cout << "Suboptimal Target Sequence 1: " << swaligned_seq1 << std::endl;
        std::cout << "Suboptimal Query Sequence  2: " << swaligned_seq2 << std::endl;
     }
    //========================================================
    //fourth condition for when l1 is larger and l2 which is a read larger than 35bp (investigate structural variants e.g., INDEL)
    //Assists the highest optimal score to NW for best and secondary scores.
    map <int, float> score4;
    float max4 = std::numeric_limits<float>::min();
    int max4_ind4;
    float ps4 ;
    float ss4 ;
    vector <int> ST4 ;
    bool multi4;
    int subopt_count;
    if(l1>l2 and d2 >= 1 and r2>0 ){
        for (int k=0; k<(l1-1) ; k++){
            if (k!= 0) {temp1 = temp1.substr(1);}
            ps4 = 0;
            ss4 = 0;
                for(int h=0; h < l2 ; h++) {
                    if (temp1[h] == temp2[h]) {
                        //cout <<"| ";
                        ps4++;
                    } else {
                        // cout<<"  ";
                        ss4++;
                    }
                }
            //cout<< "Line " <<d2+1 <<" include nucleotides "<< doc+1<<" till "<<l2 << endl;
            //cout<<endl;
            //cout << "The score is "<<static_cast<float>(ps4/l2) <<endl;
            score4.insert(make_pair(k+1, static_cast<float>(ps4/l2)));
            //cout<<static_cast<float>(ps4/l2)<<endl;
        }
        //cout << ss4 <<endl;
        //cout<< ps4 << endl;
        //float max = 0.0;
        subopt_count = 1;
        ST4.push_back(0);
        for ( auto &pair : score4) {
            if (pair.second > max4) {
                subopt_count =  ST4[0];
                //cout<<subopt_count<<endl;
                ST4.clear();
                max4 = pair.second;
                max4_ind4 = pair.first;
                ST4.push_back(max4_ind4);
            } else if ( pair.second == max4 ){
                subopt_count =  pair.first;
                multi4 = true;
                max4 = pair.second;
                ST4.push_back(pair.first);
            }
            //cout<<pair.first<<endl;
        }
        if (max4 != std::numeric_limits<int>::min() and multi4 != true ) {
            cout << "The best Identical Site is at "<< max4_ind4 << "th nc which is " <<max4*100<<"%"<<endl;
        } else if (max4 != std::numeric_limits<int>::min() and multi4 == true ) {
            cout << "The best Identical Site is "<< max4*100 <<"% at multiple loci: "<<endl;
            for ( int i=0 ; i<ST4.size() ; i+=1){
                cout<<" in locus " << ST4[i]  <<endl;
            }
        } else {
            cout << "The map is empty." << std::endl;
        }
    }
    cout<<subopt_count<<endl;
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

        string nwfing1 = fing1.substr(max4_ind4-1 );
        string nwfing3 = fing1.substr(0, max4_ind4-2);
        string nwfing2 = fing2;
        string nwfing4 = fing1.substr(subopt_count-1); //looking for index which is suboptimal location minus one

    cout << "Target Sequence - Chain One: " << nwfing3 << endl;
    cout << "Target Sequence - Chain Two Included Anchor: " << nwfing1 << endl;
    cout<<"============================================================================================================= \n";
    cout << "Query Sequence: " << nwfing2 << endl <<endl;
    cout<<"======================Showing Top Candidate Alignment===================================================== \n";

        string alignedfing4 = "";
        string alignedfing2 = "";

        //backtracing
        while (i > 0 || j > 0) {
        //while (i > 0 and j > 0) {
                //while (i == j and  i > 0 ){
                if (i > 0 && matrix[i][j] == matrix[i - 1][j] - 1) {
                    alignedfing2  = nwfing2[i - 1] + alignedfing2;
                    alignedfing4 = "-" + alignedfing4;
                    i--;
                } else if (j > 0 && matrix[i][j] == matrix[i][j - 1] - 1) {
                    alignedfing2  = "-" + alignedfing2;
                    alignedfing4 = nwfing4[j - 1] + alignedfing4;
                    j--;
                } else {
                    alignedfing4  = nwfing4[j - 1] + alignedfing4;
                    alignedfing2 = nwfing2[i - 1] + alignedfing2;
                    i--;
                    j--;
                }
            }
        //}

    if(l1>l2 and d2 >= 1 and r2>0){
        float ps4 = 0;
        float ss4 = 0;
        //fing1= fing1.substr(max4_ind4-1 );
        for( int i=0; i<d2; i++){
            doc = 35 * M;
            iam = 35 * G;
            cout<< "Line "<< i+1 <<" include Query nucleotides "<< iam+1 << " till " << doc <<endl ;
            for(int k=iam; k<doc ; k++){
                cout << nwfing1[k] << " ";
            }
            cout <<endl;
            for(int k=iam; k<doc ; k++){
                if(nwfing1[k] == nwfing2[k]){
                    cout <<"| ";
                    ps4++;
                }else{
                    cout<<": ";
                    ss4++;
                }
            }
            cout << endl;
            for(int k=iam; k<doc ; k++){
                cout << nwfing2[k] << " ";
            }
            cout <<endl;
            M = M + 1;
            G = G + 1;
        }
        cout<< "Line " <<d2+1 <<" include Query nucleotides "<< doc+1<<" till "<<l2 << endl;
        for(int i=doc ; i<l2; i++){
            cout<< nwfing1[i] << " ";
        }
        cout<<endl;
        for(int i=doc ; i<l2; i++){
            if(nwfing1[i] == nwfing2[i]) {
                cout << "| ";
                ps4++;
            }else{
                cout<< ": ";
                ss4++;
            }
        }
        cout<<endl;
        for(int i=doc ; i<l2; i++){
            cout<< nwfing2[i] << " ";
        }
        cout<<endl;
        cout << "The Identical Site is "<<(static_cast<float>(ps4/l2))*100 <<"%"<<endl;
    }
    cout<<endl;
    cout<<"====================== Needleman-Wunsch for suboptimal Alignment===================================================== \n";
        cout << "Secondary Sequence 1: " << alignedfing4 << endl;
    cout<<"==================================================================================== \n";
        cout << "Secondary Sequence 2: " << alignedfing2 << endl;
    cout<<"==================================================================================== \n";

    istrm1.close();
    istrm2.close();
}