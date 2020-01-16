def dictionary(fil):
    diz={}
    l=fil.readline()
    l=l.rstrip()
    base=l.split()
    for i in range(len(base)):
        line=fil.readline()
        line=line.rstrip()
        line=line.split()
        for j in range(len(line)):
            diz[base[i]+base[j]]=int(line[j][:len(line[j])-1])
            diz[base[j]+base[i]]=int(line[j][:len(line[j])-1])
            diz[base[i]+"-"]=-2
            diz["-"+base[i]]=-2
    return(diz)

PAM250=open("./PAM250.txt","r")
PAM=dictionary(PAM250)
BLOSUM62=open("./BLOSUM62.txt","r")
BLOSUM=dictionary(BLOSUM62)

def alignments(fasta):
    alignment=[]
    list_a=[]
    for x in range(3):
        line=fasta.readline()
        line=fasta.readline()
        line=line.rstrip()
        list_a.append(line)
        line=fasta.readline()
        line=fasta.readline()
        line=line.rstrip()
        list_a.append(line)
        alignment.append(list_a)
        list_a=[]
    return(alignment)
f=open("./alignments.fasta","r")
sequence=alignments(f)

def score(lista):
    for al in lista:
        key_list=[]
        for i in range(len(al[0])):
            key_list.append(al[0][i]+al[1][i])
        score_blosum=0
        score_pam=0

        for x in key_list:
            score_blosum=score_blosum+BLOSUM[x]
            score_pam=score_pam+PAM[x]
        print("\n",al[0],"\n",al[1],"\n","score blosum62:",score_blosum,"\n","score pam250:", score_pam)

score(sequence)
        


        

