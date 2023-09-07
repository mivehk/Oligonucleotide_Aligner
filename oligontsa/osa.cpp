
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
    istrm1.open(abspath+sequ1);
    if (istrm1.fail()){ cout<<"input of first file failed. \n"; exit (1);};
    istrm2.open(abspath+sequ2);
    if (istrm2.fail()){ cout<<"input of second file failed. \n"; exit (1);};

    //cout << count_lines(istrm1) <<endl;
    istrm1.clear();
    istrm1.seekg(0);
    //cout << count_lines(istrm2) <<endl;
    istrm2.clear();
    istrm2.seekg(0);


    if( getline(istrm1 ,firstline1 , '\n' )){
        lineignore1 = firstline1;
        //cout<<lineignore1 <<endl;
    }
    if(lineignore1[0] == '>'){
        istrm1.clear();
        istrm1.seekg(0);
        getline(istrm1,firstline1, '\n' );
        if(getline(istrm1, firstline1, '\n')){
           fing1 = firstline1;
        }
    } else{
        istrm1.clear();
        istrm1.seekg(0);
        if(getline(istrm1, firstline1, '\n')){
            fing1 = firstline1;
        }
    }

    if( getline(istrm2 ,firstline2 , '\n' )){
        lineignore2 = firstline2;
        //cout<<lineignore2 <<endl;
    }

    if(lineignore2[0] == '>'){
        istrm1.clear();
        istrm2.seekg(0);
        getline(istrm2,firstline2 , '\n');
        if(getline(istrm2, firstline2, '\n')){
            fing2 = firstline2;
        }
    } else{
        istrm1.clear();
        istrm2.seekg(0);
        if(getline(istrm2, firstline2, '\n')){
            fing2 = firstline2;
        }
    }

    int l1 = fing1.length();
    int l2 = fing2.length();
    double d1 = floor(static_cast<double>(l1) / 35);
    int r1 = l1 % 35;
    double d2 = floor(static_cast<double>(l2) / 35);
    int r2 = l2 % 35;
    int M = 1;
    int G = 0;
    int doc = 35 ;
    int iam = 0;


   /* cout<<l1 <<endl;
    cout<<l2 <<endl;
    cout<<d1 <<endl;
    cout<<d2 <<endl;
    cout<<r1 <<endl;
    cout<<r2 <<endl;*/

   string fing3= fing1 ;
   string fing4 = fing2;

   /*
    substr hop = 0, 1, ... 10
    Index position = 0, 1, ... 10
    K count = 1, 2, ... 11
   */

    //started conditional printing==============================
    //first and sixth condition for equal oligo length
    if(l1 == l2 and d1<1){
        double ps1=0 ;
        double ss1=0 ;
        cout<<"oligonucleotides one till " <<l1 <<endl;
        //cout << "yek \n";
        for(int i=0; i<l2 ; i++){
           cout <<fing1[i] << " ";
        }
        cout<<endl;
        for(int k=0 ; k<l1; k++ ){
            if(fing1[k]==fing2[k]){
                cout<<"| ";
                ps1= ps1+1;
            }else{
                cout<<"  ";
                ss1++;
            }
        }
      cout<< endl;
        for(int j=0 ;j<l2; j++){
            cout <<fing2[j] << " ";
        }
        cout<< endl;
        cout<<"The score is "<< static_cast<double>(ps1/l1) << endl;
    }
    if(l2 == l1 and d1>1 and r1>0){
        double ps6=0;
        double ss6=0;
        //cout<<"shish \n";
        for(int i=0; i<d1; i++){
            doc = 35 * M;
            iam = 35 * G;
            cout << "Line "<< i+1 <<" include nucleotides "<<iam+1 <<" till "<<doc <<endl;
            for(int j=iam; j<doc ; j++){
                cout<< fing1[j] <<" ";
            }
            cout<< endl;
            for(int k=iam; k<doc; k++) {
                if (fing1[k] == fing2[k]) {
                    cout << "| ";
                    ps6++;
                } else {
                    cout << "  ";
                    ss6++;
                }
            }
            cout<<endl;
            for(int a=iam; a<doc; a++){
                cout<<fing2[a] << " ";
            }
            cout<<endl;
            M = M + 1;
            G = G + 1;
        }
        cout <<"Line "<< d1+1 <<" include nucleotides "<<doc+1<<" till "<<l2 <<endl;
        for(int i=doc ; i<l1; i++){
            cout<< fing1[i] << " ";
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            if(fing1[i] == fing2[i]) {
                cout << "| ";
                ps6++;
            }else{
                cout<< "  ";
                ss6++;
            }
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            cout<< fing2[i] << " ";
        }
        cout<<endl;
        cout << "The score is "<<static_cast<double>(ps6/l2) <<endl;
    }
    //===================================================
    //second condition for when l1 is larger and l2 is short
    map<int, double> score2;
    double max2 = std::numeric_limits<double>::min();
    int max2_ind2 = 0;
    if(l1>l2 and d2<1 ){
        for (int k=1; k<(l1-1); k++){
            //Hopping nucleotides in each k iteration along nucleic acid to find similarity and score.
            if (k !=1) {fing3 = fing3.substr(1);}
            double ps2 = 0;
            double ss2 = 0;
            //dropped fs-printing from here
            for(int m=0 ; m<l2; m++ ){
                if(fing3[m]==fing4[m]){
                    //cout<<"| ";
                    ps2 ++;
                }else{
                    //cout<<"  ";
                    ss2 ++;
                }
            }
            //cout<< endl;
            //dropped ss-printing from here
            score2.insert(make_pair(k, static_cast<double>(ps2/l2)));
            cout<<static_cast<double>(ps2/l2)<<endl;
            //cout<<ps2<<endl;
            //cout<<ss2<<endl;
        }
        for (const auto& pair : score2) {
            if (pair.second > max2) {
                max2 = pair.second;
                max2_ind2 = pair.first;
            }
            //cout<<pair.first<<endl;
        }
        //index 0 in seq is nc one which is at zero hop, so k one is at first bit
        if (max2 != std::numeric_limits<int>::min()) {
            cout << "The best score is at "<< max2_ind2-1 << "th nc which is " <<max2 <<endl;
        } else {
            cout << "The map is empty." << std::endl;
        }
    }
    if(l1>l2 and d2<1){
        fing1= fing1.substr(max2_ind2-1 );
        double ps2 = 0;
        double ss2 = 0;
        cout<<"oligonucleotides one till " <<l2 <<endl;
        //cout << "do \n";
        for(int i=0; i<l2 ; i++){
            cout <<fing1[i] << " ";
        }
        cout<<endl;
        for(int k=0 ; k<l2; k++ ){
            if(fing1[k]==fing2[k]){
                cout<<"| ";
                ps2++;
            }else{
                cout<<"  ";
                ss2++;
            }
        }
        cout<< endl;
        for(int j=0 ;j<l2; j++){
            cout <<fing2[j] << " ";
        }
        cout<< endl;
        cout << "The calculated score is "<<static_cast<double>(ps2/l2) <<endl;
     }
    //========================================================
    //fourth condition for when l1 is larger and l2 is not short
    map <int, double> score4;
    double max4 = std::numeric_limits<double>::min();
    int max4_ind4 = 0;
    if(l1>l2 and d2 >= 1 and r2>0 ){
        for (int k=1; k<(l1-1) ; k++){
            if (k!= 1) {fing3= fing3.substr(1);}
            double ps4;
            double ss4;
                for(int m=0; m<l2 ; m++){
                    if(fing3[m] == fing4[m]){
                        //cout <<"| ";
                        ps4++;
                    }else{
                       // cout<<"  ";
                        ss4++;
                    }
                }
            //cout<< "Line " <<d2+1 <<" include nucleotides "<< doc+1<<" till "<<l2 << endl;
            //cout<<endl;
            //cout << "The score is "<<static_cast<double>(ps4/l2) <<endl;
            score4.insert(make_pair(k, static_cast<double>(ps4/l2)));
            //cout<<static_cast<double>(ps4/l2)<<endl;
        }
        //double max = 0.0;

        for (const auto& pair : score4) {
            if (pair.second > max4) {
                max4 = pair.second;
                max4_ind4 = pair.first;
            }
        }
        //cout << "The score is "<< max <<endl;
        cout << "The best score is at "<< max4_ind4-1 << "th nc which is " <<max4 <<endl;
    }
    if(l1>l2 and d2 >= 1 and r2>0){
         double ps4;
         double ss4;
         fing1= fing1.substr(max4_ind4-1 );
         //cout<<"chohaar \n";
         for( int i=0; i<d2; i++){
             doc = 35 * M;
             iam = 35 * G;
             cout<< "Line "<< i+1 <<" include nucleotides "<< iam+1 << " till " << doc <<endl ;
             for(int k=iam; k<doc ; k++){
                 cout << fing1[k] << " ";
             }
             cout <<endl;
             for(int k=iam; k<doc ; k++){
                 if(fing1[k] ==fing2[k]){
                     cout <<"| ";
                     ps4++;
                 }else{
                     cout<<"  ";
                     ss4++;
                 }
             }
             cout << endl;
             for(int k=iam; k<doc ; k++){
                 cout << fing2[k] << " ";
             }
             cout <<endl;
             M = M + 1;
             G = G + 1;
         }
         cout<< "Line " <<d2+1 <<" include nucleotides "<< doc+1<<" till "<<l2 << endl;
         for(int i=doc ; i<l2; i++){
             cout<< fing1[i] << " ";
         }
         cout<<endl;
         for(int i=doc ; i<l2; i++){
             if(fing1[i] == fing2[i]) {
                 cout << "| ";
                 ps4++;
             }else{
                 cout<< "  ";
                 ss4++;
             }
         }
         cout<<endl;
         for(int i=doc ; i<l2; i++){
             cout<< fing2[i] << " ";
         }
         cout<<endl;
         cout << "The score is "<<static_cast<double>(ps4/l2) <<endl;
     }
    //========================================================
    //Third condition for when l2 is larger and l1 is short
    map <int, double> score3;
    double max3 = std::numeric_limits<double>::min();
    int max3_ind3 = 0;
    if(l1<l2 and d1<1){
        for (int k=1; k<(l2-1) ; k++) {
           if (k != 1) {fing4 = fing4.substr(1);}
            double ps3 = 0;
            double ss3 = 0;
            //cout<<"oligonucleotides one till " <<l1 <<endl;
            //fs is not needed
            for(int k=0 ; k<l1; k++ ){
                if(fing3[k]==fing4[k]){
                    //cout<<"| ";
                    ps3 ++;
                }else{
                    //cout<<"  ";
                    ss3 ++;
                }
            }
            //cout<< endl;
           //ss is not needed
            //cout << "The score is "<<static_cast<double>(ps3/l1) <<endl;
            score3.insert(make_pair(k, static_cast<double>(ps3/l1)));
        }
        for (const auto& pair : score3) {
            if (pair.second > max3) {
                max3 = pair.second;
                max3_ind3 = pair.first;
            }
        }
        cout << "The best score is at "<< max3_ind3-1 << "th hop which is " << max3 <<endl;
    }
    if(l1<l2 and d1<1){
        fing2 = fing2.substr(max3_ind3-1);
        double ps3 = 0;
        double ss3 = 0;
        cout<<"oligonucleotides one till " <<l1 <<endl;
        //cout << "se \n";
        for(int i=0; i<l1 ; i++){
            cout <<fing1[i] << " ";
        }
        cout<<endl;
        for(int k=0 ; k<l1; k++ ){
            if(fing1[k]==fing2[k]){
                cout<<"| ";
                ps3++;
            }else{
                cout<<"  ";
                ss3++;
            }
        }
        cout<< endl;
        for(int j=0 ;j<l1; j++){
            cout <<fing2[j] << " ";
        }
        cout<< endl;
        cout << "The score is "<<static_cast<double>(ps3/l1) <<endl;
    }
    //========================================================
    //Fifth condition for when l2 is larger and l1 is not short
    map <int, double> score5;
    double max5 = 0.0; //std::numeric_limits<double>::min();
    int max5_ind5 = 0;
    if(l1<l2 and d1 >= 1 and r1>0 ){
        for (int k=1; k<(l2-1) ; k++) {
            if (k != 1 ) {fing4 = fing4.substr(1);}
            double ps5 =0;
            double ss5 =0;
            //cout<<"panj \n";
            //for( int i=0; i<d1; i++){
               // doc = 35 * M;
                //iam = 35 * G;
                //cout<< "Line "<< i+1 <<" include nucleotides "<< iam+1 << " till " << doc <<endl ;
                //cout <<endl;
                //for(int k=iam; k<doc ; k++){
                for (int m=0 ; m <l1; m++){
                    if(fing3[m] ==fing4[m]){
                        //cout <<"| ";
                        ps5++;
                    }else{
                        //cout<<"  ";
                        ss5++;
                    }
                }
                //cout << endl;
                //cout <<endl;
               // M = M + 1;
                //G = G + 1;
           // }
            //cout<< "Line " <<d1+1 <<" include nucleotides "<< doc+1<<" till "<<l1 << endl;
            //for(int i=doc ; i<l1; i++){
                //cout<< fing1[i] << " ";
            //}
           // cout<<endl;
           // for(int i=doc ; i<l1; i++){
               // if(fing1[i] == fing2[i]) {
                    //cout << "| ";
                  //  ps5++;
             //   }else{
                   // cout<< "  ";
                  //  ss5++;
              //  }
           // }
            //cout<<endl;
           // for(int i=doc ; i<l1; i++){
                //cout<< fing2[i] << " ";
           // }
            //cout<<endl;
            //cout << "The score is "<<static_cast<double>(ps5/l1) <<endl;
            score5.insert(make_pair(k, static_cast<double>(ps5/l1)));
        }
        for (const auto& pair : score5) {
            if (pair.second > max5) {
                max5 = pair.second;
                max5_ind5 = pair.first;
            }
        }
        cout << "The best score is at "<< max5_ind5-1 << "th hop which is " <<max5 <<endl;
    }
    if(l2>l1 and d1 >= 1 and r1>0){
        fing2= fing2.substr(max5_ind5-1 );
        double ps5 =0;
        double ss5 =0;
        //cout<<"panj \n";
        for( int i=0; i<d1; i++){
            doc = 35 * M;
            iam = 35 * G;
            cout<< "Line "<< i+1 <<" include nucleotides "<< iam+1 << " till " << doc <<endl ;
            for(int k=iam; k<doc ; k++){
                cout << fing1[k] << " ";
            }
            cout <<endl;
            for(int k=iam; k<doc ; k++){
                if(fing1[k] ==fing2[k]){
                    cout <<"| ";
                    ps5++;
                }else{
                    cout<<"  ";
                    ss5++;
                }
            }
            cout << endl;
            for(int k=iam; k<doc ; k++){
                cout << fing2[k] << " ";
            }
            cout <<endl;
            M = M + 1;
            G = G + 1;
        }
        cout<< "Line " <<d1+1 <<" include nucleotides "<< doc+1<<" till "<<l1 << endl;
        for(int i=doc ; i<l1; i++){
            cout<< fing1[i] << " ";
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            if(fing1[i] == fing2[i]) {
                cout << "| ";
                ps5++;
            }else{
                cout<< "  ";
                ss5++;
            }
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            cout<< fing2[i] << " ";
        }
        cout<<endl;
        cout << "The score is "<<static_cast<double>(ps5/l1) <<endl;
    }
    cout<<endl;
    cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n";
    istrm1.close();
    istrm2.close();
}