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

