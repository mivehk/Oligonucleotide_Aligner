#include <iostream>
#include "osa_cin.h"
#include <string>
using namespace std;

int main() {
    string s1;
    //char s2[1025];
    string s2;
    cout<< "Please enter your first sequence: ";
    //cin >> s1;
    getline( cin, s1 ,'\n' );
/*    while (!s1.empty() && s1.back() =='\n' ){
        cin.ignore(1);
    }*/
    //cout<< s1<<endl;
    int l1 = s1.length();
    //cout<< l1<<endl;
    cout<< "Please enter your second sequence: ";
    //cin >> s2;
    //cin.getline(s2 ,1025 , '\n');
    getline( cin, s2 );
    /*The line "cin.getline(s1,101,'\n')" uses cin member function of c-string "getline()" reads into a character array "char[]" until \n or 101 char,
     * while the line "getline(cin,s2,'\n')" uses global getline() function of standard string class (<string> header) to read into a string object.*/
    while (!s2.empty() && s2.back() =='\n' ){
        cin.ignore(1);
    }
    //cout<< s2<<endl;
    //int l2 = strlen(s2);
    int l2 = s2.length();
    //cout<< l2<<endl;

    compare_sequences(s1,s2);
    return 0;
}
