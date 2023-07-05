
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


/*    cout<<l1 <<endl;
    cout<<l2 <<endl;
    cout<<d1 <<endl;
    cout<<d2 <<endl;
    cout<<r1 <<endl;
    cout<<r2 <<endl;*//*
*/
    if(l1 == l2 and d1 < 1){
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
    if(l1 > l2 and d2 < 1){
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
        cout << "The score is "<<static_cast<double>(ps2/l2) <<endl;
    }
    if(l1 < l2 and d1 < 1){
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
    if(l1>l2 and d2 >=1 and r2>0){
        double ps4;
        double ss4;
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
    if (l2 > l1 and d1 >=1 and r1>0){
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
    if (l2 == l1 and d1 >1 and r1>0) {
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
    cout<<endl;
    cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n";

    istrm1.close();
    istrm2.close();

}
