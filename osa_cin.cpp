//
// Created by kayvan mivehnejad on 6/21/23.
//

#include <iostream>
#include <string>
#include <cmath>
#include "osa_cin.h"
#include <cstdlib>
#include <vector>
using namespace std;

void compare_sequences(string sequ1 ,string sequ2 ) {

    cout << endl;
    cout << "First entered sequence alignment against second entered sequence" << endl;

    int l1 = sequ1.length();
    int l2 = sequ2.length();
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
    cout<<r2 <<endl;*/
    if(l1 == l2 and d1 < 1){
        double ps1=0;
        double ss1=0;
        cout<<"oligonucleotides one till " <<l1 <<endl;
        //cout << "yek \n";
        for(int i=0; i<l2 ; i++){
            cout <<sequ1[i] << " ";
        }
        cout<<endl;
        for(int k=0 ; k<l1; k++ ){
            if(sequ1[k]==sequ2[k]){
                cout<<"| ";
                ps1++;
            }else{
                cout<<"  ";
                ss1++;
            }
        }
        cout<< endl;
        for(int j=0 ;j<l2; j++){
            cout <<sequ2[j] << " ";
        }
        cout<< endl;
        cout<<"The score is "<< static_cast<double>(ps1/l1) << endl;
    }
    if(l1 > l2 and d2 < 1){
        double ps2 = 0;
        double ss2=0;
        cout<<"oligonucleotides one till " << l2 <<endl;
        //cout << "do \n";
        for(int i=0; i<l2 ; i++){
            cout <<sequ1[i] << " ";
        }
        cout<<endl;
        for(int k=0 ; k<l2; k++ ){
            if(sequ1[k]==sequ2[k]){
                cout<<"| ";
                ps2++;
            }else{
                cout<<"  ";
                ss2++;
            }
        }
        cout<< endl;
        for(int j=0 ;j<l2; j++){
            cout <<sequ2[j] << " ";
        }
        cout<< endl;
        cout<<"The score is "<< static_cast<double>(ps2/l2) << endl;
    }
    if(l1 < l2 and d1 < 1){
        double ps3=0;
        double ss3=0;
        cout<<"oligonucleotides one till " <<l1 <<endl;
        //cout << "se \n";
        for(int i=0; i<l1 ; i++){
            cout <<sequ1[i] << " ";
        }
        cout<<endl;
        for(int k=0 ; k<l1; k++ ){
            if(sequ1[k]==sequ2[k]){
                cout<<"| ";
                ps3++;
            }else{
                cout<<"  ";
                ss3++;
            }
        }
        cout<< endl;
        for(int j=0 ;j<l1; j++){
            cout <<sequ2[j] << " ";
        }
        cout<< endl;
        cout<<"The score is "<< static_cast<double>(ps3/l1) << endl;
    }
    if(l1>l2 and d2 >=1 and r2>0){
        double ps4 =0;
        double ss4=0;
        //cout<<"chohaar \n";
        for( int i=0; i<d2; i++){
            doc = 35 * M;
            iam = 35 * G;
            cout<< "Line "<< i+1 <<" include nucleotides "<< iam+1 << " till " << doc <<endl ;
            for(int k=iam; k<doc ; k++){
                cout << sequ1[k] << " ";
            }
            cout <<endl;
            for(int k=iam; k<doc ; k++){
                if(sequ1[k] ==sequ2[k]){
                    cout <<"| ";
                    ps4++;
                }else{
                    cout<<"  ";
                    ss4++;
                }
            }
            cout << endl;
            for(int k=iam; k<doc ; k++){
                cout << sequ2[k] << " ";
            }
            cout <<endl;
            M = M + 1;
            G = G + 1;
        }
        cout<< "Line " <<d2+1 <<" include nucleotides "<< doc+1<<" till "<<l2 << endl;
        for(int i=doc ; i<l2; i++){
            cout<< sequ1[i] << " ";
        }
        cout<<endl;
        for(int i=doc ; i<l2; i++){
            if(sequ1[i] == sequ2[i]) {
                cout << "| ";
                ps4++;
            }else{
                cout<< "  ";
                ss4++;
            }
        }
        cout<<endl;
        for(int i=doc ; i<l2; i++){
            cout<< sequ2[i] << " ";
        }
        cout<<endl;
        cout<<"The score is "<< static_cast<double>(ps4/l2) << endl;
    }
    if (l2 > l1 and d1 >=1 and r1>0){
        double ps5=0;
        double ss5=0;
        //cout<<"panj \n";
        for( int i=0; i<d1; i++){
            doc = 35 * M;
            iam = 35 * G;
            cout<< "Line "<< i+1 <<" include nucleotides "<< iam+1 << " till " << doc <<endl ;
            for(int k=iam; k<doc ; k++){
                cout << sequ1[k] << " ";
            }
            cout <<endl;
            for(int k=iam; k<doc ; k++){
                if(sequ1[k] ==sequ2[k]){
                    cout <<"| ";
                    ps5++;
                }else{
                    cout<<"  ";
                    ss5++;
                }
            }
            cout << endl;
            for(int k=iam; k<doc ; k++){
                cout << sequ2[k] << " ";
            }
            cout <<endl;
            M = M + 1;
            G = G + 1;
        }
        cout<< "Line " <<d1+1 <<" include nucleotides "<< doc+1<<" till "<<l1 << endl;
        for(int i=doc ; i<l1; i++){
            cout<< sequ1[i] << " ";
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            if(sequ1[i] == sequ2[i]) {
                cout << "| ";
                ps5++;
            }else{
                cout<< "  ";
                ss5++;
            }
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            cout<< sequ2[i] << " ";
        }
        cout<<endl;
        cout<<"The score is "<< static_cast<double>(ps5/l1) << endl;
    }
    if (l2 == l1 and d1 >1 and r1>0) {
        double ps6 =0;
        double ss6 =0;
        //cout<<"shish \n";
        for(int i=0; i<d1; i++){
            doc = 35 * M;
            iam = 35 * G;
            cout << "Line "<< i+1 <<" include nucleotides "<<iam+1 <<" till "<<doc <<endl;
            for(int j=iam; j<doc ; j++){
                cout<< sequ1[j] <<" ";
            }
            cout<< endl;
            for(int k=iam; k<doc; k++) {
                if (sequ1[k] == sequ2[k]) {
                    cout << "| ";
                    ps6++;
                } else {
                    cout << "  ";
                    ss6++;
                }
            }
            cout<<endl;
            for(int a=iam; a<doc; a++){
                cout<<sequ2[a] << " ";
            }
            cout<<endl;
            M = M + 1;
            G = G + 1;
        }
        cout <<"Line "<< d1+1 <<" include nucleotides "<<doc+1<<" till "<<l2 <<endl;
        for(int i=doc ; i<l1; i++){
            cout<< sequ1[i] << " ";
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            if(sequ1[i] == sequ2[i]) {
                cout << "| ";
                ps6++;
            }else{
                cout<< "  ";
                ss6++;
            }
        }
        cout<<endl;
        for(int i=doc ; i<l1; i++){
            cout<< sequ2[i] << " ";
        }
        cout<<"The score is "<< static_cast<double>(ps6/l1) << endl;
        cout<<endl;
    }
    cout<<endl;
    cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n";

}