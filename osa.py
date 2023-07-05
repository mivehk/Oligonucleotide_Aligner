import os
import time
os.chdir('/Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK04/HW4')
os.environ['COLUMNS'] = str(1000)

print()


def compare_sequences(sequ1: str, sequ2: str):
    '''scrupulous words toward others has consequences! Don't! '''
    print(f"Top Sequence is called {sequ1} and Bottom sequence is {sequ2}")
    with open(sequ1, 'r') as data1, open(sequ2, 'r') as data2:
        # print(sum(1 for i in data1))
        # print(sum(1 for i in data2))
        a1 = data1.readline()
        # print(a1)
        if a1[0] == '>':
            data1.seek(0)
            ignore = data1.readline()
            # data1.seek(0)
            s1 = data1.read().rstrip('\n')
            '''it is important to strip any special character like "\n" otherwise len return extra char '''
            # print(s1)
        else:
            data1.seek(0)
            s1 = data1.read().rstrip('\n')

        a2 = data2.readline()
        # print(a2)
        if a2[0] == '>':
            data2.seek(0)
            ignore = data2.readline()
            s2 = data2.read().rstrip('\n')
            # print(s2)
        else:
            data2.seek(0)
            s2 = data2.read().rstrip('\n')
        l1 = len(s1)
        l2 = len(s2)
        d1 = l1 // 35
        r1 = l1 % 35
        d2 = l2 // 35
        r2 = l2 % 35
        M = 1
        G = 0
        # print(l1)
        # print(l2)
        # print(d1)
        # print(d2)
        # print(r1)
        # print(r2)
        if l1 == l2 and d1 < 1:
            ps1 = 0
            ss1 = 0
            print(f"oligonucleotides one till {l1}")
            # print('yek')
            for i in range(0, l2):
                # ali.append(s1[i])
                print(s1[i], end=' ')
            # sys.stdout.reconfigure(wrap=0)
            print()
            for k in range(0, l1):
                if s1[k] == s2[k]:
                    print("|", end=' ')
                    ps1 += 1
                else:
                    print(" ", end=' ')
                    ss1 += 1
            print()
            for j in s2:
                print(j, end=' ')
            print()
            print(f"The score is {ps1/l1}")
            print(' ', end='\r')
        if l1 > l2 and d2 < 1:
            ps2 = 0
            ss2 = 0
            # print('do')
            print(f"oligonucleotides one till {l2}")
            for i in range(0, l2):
                # ali.append(s1[i])
                print(s1[i], end=' ')
            # sys.stdout.reconfigure(wrap=0)
            print()
            for k in range(0, l2):
                if s1[k] == s2[k]:
                    print("|", end=' ')
                    ps2 += 1
                else:
                    print(" ", end=' ')
                    ss2 += 1
            print()
            for j in s2:
                print(j, end=' ')
            print()
            print(f"The score is {ps2/l2}")
            print(' ', end='\r')
        if l2 > l1 and d1 < 1:
            print(f"oligonucleotides one till {l1}")
            ps3 = 0
            ss3 = 0
            # print('se')
            for i in range(0, l1):
                # ali.append(s1[i])
                print(s2[i], end=' ')
            # sys.stdout.reconfigure(wrap=0)
            print()
            for k in range(0, l1):
                if s1[k] == s2[k]:
                    print("|", end=' ')
                    ps3 += 1
                else:
                    print(" ", end=' ')
                    ss3 += 1
            print()
            for j in s1:
                print(j, end=' ')
            print()
            print(f"The score is {ps3/l1}")
            print(' ', end='\r')
        if l1 > l2 and d2 > 0:
            ps4 = 0
            ss4 = 0
            # print('char')
            for i in range(0, d2):
                doc = 35 * M
                iam = 35 * G
                # M=1
                print(f"Line {i+1} include nucleotides {iam+1} till {doc}")
                for j in range(iam, doc):
                    print(s1[j], end=' ')
                print()
                for k in range(iam, doc):
                    if s1[k] == s2[k]:
                        print('|', end=' ')
                        ps4 += 1
                    else:
                        print(' ', end=' ')
                        ss4 += 1
                print()
                for a in range(iam, doc):
                    print(s2[a], end=' ')
                print()
                print()
                M = M + 1
                G = G + 1
            print(f"Line {i+2} include nucleotides {doc+1} till {l2}")
            for i in range(doc, l2):
                print(s1[i], end=' ')
            print()
            for i in range(doc, l2):
                if s1[i] == s2[i]:
                    print('|', end=' ')
                    ps4 += 1
                else:
                    print(' ', end=' ')
                    ss4 += 1
            print()
            for i in range(doc, l2):
                print(s2[i], end=' ')
            print()
            print(f"The score is {ps4/l2}")
            print()
        if l2 > l1 and d1 > 0:
            ps5 = 0
            ss5 = 0
            # print('panj')
            for i in range(0, d1):
                doc = 35 * M
                iam = 35 * G
                # M=1
                print(f"Line {i+1} include nucleotides {iam+1} till {doc}")
                for j in range(iam, doc):
                    print(s1[j], end=' ')
                print()
                for k in range(iam, doc):
                    if s1[k] == s2[k]:
                        print('|', end=' ')
                        ps5 += 1
                    else:
                        print(' ', end=' ')
                        ss5 += 1
                print()
                for a in range(iam, doc):
                    print(s2[a], end=' ')
                print()
                print()
                M = M + 1
                G = G + 1
            print(f"Line {i+2} include nucleotides {doc+1} till {l1}")
            # print(doc)
            # print(l1)
            for i in range(doc, l1):
                print(s1[i], end=' ')
            print()
            # for u in range(doc, l1):
            #     print(s2[u] , end='-')
            for i in range(doc, l1):
                if s1[i] == s2[i]:
                    print('|', end=' ')
                    ps5 += 1
                else:
                    print(' ', end=' ')
                    ss5 += 1
            print()
            for i in range(doc, l1):
                print(s2[i], end=' ')
            print()
            print(f"The score is {ps5/l1}")
            print()
        if l2 == l1 and d1 > 0:
            ps6 = 0
            ss6 = 0
            # print('shish')
            for i in range(0, d1):
                doc = 35 * M
                iam = 35 * G
                # M=1
                print(f"Line {i+1} include nucleotides {iam+1} till {doc}")
                for j in range(iam, doc):
                    print(s1[j], end=' ')
                print()
                for k in range(iam, doc):
                    if s1[k] == s2[k]:
                        print('|', end=' ')
                        ps6 += 1
                    else:
                        print(' ', end=' ')
                        ss6 += 1
                print()
                for a in range(iam, doc):
                    print(s2[a], end=' ')
                print()
                print()
                M = M + 1
                G = G + 1
            print(f"Line {i+2} include nucleotides {doc+1} till {l2}")
            for i in range(doc, l1):
                print(s1[i], end=' ')
            print()
            for i in range(doc, l1):
                if s1[i] == s2[i]:
                    print('|', end=' ')
                    ps6 += 1
                else:
                    print(' ', end=' ')
                    ss6 += 1
            print()
            for i in range(doc, l1):
                print(s2[i], end=' ')
                print()
            print(f"The score is {ps6/l1}")
            print()
        print()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
