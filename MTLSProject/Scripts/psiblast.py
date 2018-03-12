#PSI-BLAST

#Al final tendremos un PSSM para cada proteina. En este caso, 42 PSSMs. 
#para todo lo que acabe en fasta do psi-blast + flag
#database es swissprot
#number of iterations is 3

#################################################################
#          SAVE EACH PROTEIN IN A SEPARATED FASTA FILE          #
#################################################################

filehandle=open("membrane-beta_3state.3line.txt","r")
text=filehandle.read().splitlines()  

for i in range(len(text)):
    seq=[]
    if text[i].startswith(">"):
        seq.append(text[i])
        seq.append(text[i+1])  
        writehandle=open("../Datasets/FASTAfiles/"+seq[0]+".fasta","w")
        writehandle.write(seq[0]+"\n"+seq[1])
        

filehandle.close()
writehandle.close()
    
    
    

if __name__ == "__main__":
    
    print("")

