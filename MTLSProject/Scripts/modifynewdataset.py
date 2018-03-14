#Adjust the new protein dataset to work with the script I have coded.

filehandle=open("newdataset.txt","r")
filelines=filehandle.read().splitlines() 
writehandle=open("modifiednewdataset.txt","w")
writehandle2=open("fastadataset.fasta","w")

for i in range(len(filelines)):
    if filelines[i].startswith(">"):
        writehandle2.write(filelines[i]+"\n")
        writehandle2.write(filelines[i+1]+"\n")
        filelines[i+2]=filelines[i+2].replace("X","i")
        filelines[i+2]=filelines[i+2].replace("m","M")
        filelines[i+2]=filelines[i+2].replace("O","o")
        filelines[i+2]=filelines[i+2].replace("I","i")
        
    writehandle.write(filelines[i]+"\n")
