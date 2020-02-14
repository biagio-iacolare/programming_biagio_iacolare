import input_data
seq1=input_data.seq1
seq2=input_data.seq2
BLOSUM52=input_data.BLOSUM52
def matrices_F(seq1,seq2,dictionary,d):
    F=[]
    P=[]
    for i in range(len(seq2)+1):
        line=[]
        line.append(i*d)
        move=[]
        if i==0:
            move.append("0")
        else:
            move.append("u")
        for j in range(len(seq1)+1):
            if i==0 and j!=0:
                line.append(j*d)
                move.append("l")
            elif j!=0:
                best_three=[]
                left=line[j-1]+d
                upper=F[i-1][j]+d
                diagonal=F[i-1][j-1]+dictionary[seq2[i-1]+seq1[j-1]]
                best_three.append(left)
                best_three.append(upper)
                best_three.append(diagonal)
                line.append(max(best_three))
                if max(best_three)==left:
                    move.append("l")
                elif max(best_three)==upper:
                    move.append("u")
                elif max(best_three)==diagonal:
                    move.append("d")
        F.append(line)
        P.append(move)
    return(F)
F=matrices_F(seq1,seq2,BLOSUM52,-2)
#I'm repeating the same function but returning the other matrix
def matrices_P(seq1,seq2,dictionary,d):
    F=[]
    P=[]
    for i in range(len(seq2)+1):
        line=[]
        line.append(i*d)
        move=[]
        if i==0:
            move.append("0")
        else:
            move.append("u")
        for j in range(len(seq1)+1):
            if i==0 and j!=0:
                line.append(j*d)
                move.append("l")
            elif j!=0:
                best_three=[]
                left=line[j-1]+d
                upper=F[i-1][j]+d
                diagonal=F[i-1][j-1]+dictionary[seq2[i-1]+seq1[j-1]]
                best_three.append(left)
                best_three.append(upper)
                best_three.append(diagonal)
                line.append(max(best_three))
                if max(best_three)==left:
                    move.append("l")
                elif max(best_three)==upper:
                    move.append("u")
                elif max(best_three)==diagonal:
                    move.append("d")
        F.append(line)
        P.append(move)
    return(P)
P=matrices_P(seq1,seq2,BLOSUM52,-2)
print(F,P)
def needleman_wunsch_alignment(F,P,seq1,seq2):
    a=len(seq2)
    b=len(seq1)
    s1=""
    s2=""
    while P[a][b]=="0":
        if P[a][b]=="d":
            s1=s1+seq1[b-1]
            s2=s2+seq2[a-1]
            a=a-1
            b=b-1
        elif P[a][b]=="l":
            s1=s1+seq1[b-1]
            s2=s2+"-"
            b=b-1
        elif P[a][b]=="u":
            s1=s1+"-"
            s2=s2+seq2[a-1]
            a=a-1
    s1=s1[::-1]
    s2=s2[::-1]
    score=F[len(seq2)][len(seq1)]
    return(s1,s2,score)
print(needleman_wunsch_alignment(F,P,seq1,seq2))


