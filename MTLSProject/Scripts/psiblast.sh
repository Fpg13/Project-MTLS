#################################################################
#                      PSI-BLAST BASH SCRIPT                    #
#################################################################

cd ../Datasets/FASTAfiles

for file in *.fasta 

do

echo "Running psiblast on $file at $(date)..."
time psiblast -query $file -db /../../../../../../scratch/uniref90.fasta -num_iterations 3 -evalue 0.001 -num_threads 8 -out ../PSIBLASTfiles/$file.psiblast -out_ascii_pssm ../PSSMfiles/$file.pssm 
echo "Finished running psiblast on $file at $(date)."

done

echo 'PSI-BLAST run is complete'
