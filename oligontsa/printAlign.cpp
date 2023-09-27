//
// Created by kayvan mivehnejad on 9/26/23.
//

#include <iostream>
#include "printAlign.h"
#include <cmath>

using namespace std;

void print_Alignment(string fing1 , string fing2){

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

    float ps6 = 0;
    float ss6 = 0;

    for (int i = 0; i < d2; i++) {
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
    cout << "Line " << d2 + 1 << " include Query nucleotides " << doc + 1 << " till " << l2 << endl;
    for (int i = doc; i < l2; i++) {
        cout << fing1[i] << " ";
    }
    cout << endl;
    for (int i = doc; i < l2; i++) {
        if (fing1[i] == fing2[i]) {
            cout << "| ";
            ps6++;
        } else {
            cout << ": ";
            ss6++;
        }
    }
    cout << endl;
    for (int i = doc; i < l2; i++) {
        cout << fing2[i] << " ";
    }
    cout << endl;
    cout << "The Identical Site is " << (static_cast<float>(ps6 / l2)) * 100 << "%" << endl;
}

