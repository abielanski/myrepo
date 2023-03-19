#STEP 1#

#Import Python libraries
import gzip
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

#Generate a new directory/folder to store outputfiles
os.system("mkdir PipelineProject_Anna_Bielanski")

#Change to new directory to run commands
os.system("cd ~/PipelineProject_Anna_Bielanski")

#Open log file for writing & to store output
logfile = open("PipelineProject.log","w")

#Use wget to retrieve the transcriptomes from the 2 patients donors from SRA
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045")

#Convert the transcriptomes to paired-end fastq files using fastq-dump command
os.system("fastq-dump -I --split-files SRR5660030")
os.system("fastq-dump -I --split-files SRR5660033")
os.system("fastq-dump -I --split-files SRR5660044")
os.system("fastq-dump -I --split-files SRR5660045")

#STEP 2 (TRACK 2)#

#Run this python code to retrieve HCMV genome and convert to FASTA
exec(open("Accession_to_FASTA.py").read())

#Create an index for HCMV(NCBI accession NC_006273.2) using bowtie2-build command
#This command should output a set of 6 files
os.system("bowtie2-build ref_genome_out.txt NC_006273.2")

#Map and save reads to the generated HCMV index
os.system("bowtie2 --quiet -x NC_006273.2 -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S SRR5660030map.sam --al-conc SRR5660030_mapped_%.fq")
os.system("bowtie2 --quiet -x NC_006273.2 -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S SRR5660033map.sam --al-conc SRR5660033_mapped_%.fq")
os.system("bowtie2 --quiet -x NC_006273.2 -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S SRR566044map.sam --al-conc SRR5660044_mapped_%.fq")
os.system("bowtie2 --quiet -x NC_006273.2 -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S SRR5660045map.sam --al-conc SRR5660045_mapped_%.fq")


#Calculate the number of lines for one of the paired fastq files and the mapped fastq files and divide it by 4
#in order to receive the number of read pairs in each
#Repeat step for each of the paired reads (fastq & mapped fastq)
Donor_1_before_2dpi = os.system("wc -l SRR5660030_1.fastq")
Donor_1_reads_before = (int(Donor_1_before_2dpi))/4
Donor_1_after_2dpi = os.system("wc -l SRR5660030_mapped_1.fq")
Donor_1_reads_after = (int(Donor_1_after_2dpi))/4

Donor_1_before_6dpi = os.system("wc -l SRR5660033_1.fastq")
Donor_1_reads_before_6 = (int(Donor_1_before_6dpi))/4
Donor_1_after_6dpi = os.system("wc -l SRR5660033_mapped_1.fq")
Donor_1_reads_after_6 = (int(Donor_1_after_6dpi))/4

Donor_3_before_2dpi = os.system("wc -l SRR5660044_1.fastq")
Donor_3_reads_before_2 = (int(Donor_3_before_2dpi))/4
Donor_3_after_2dpi = os.system("wc -l SRR5660044_mapped_1.fq")
Donor_3_reads_after_2 = (int(Donor_3_after_2dpi))/4


Donor_3_before_6dpi = os.system("wc -l SRR5660045_1.fastq")
Donor_3_reads_before_6 = (int(Donor_3_before_6dpi))/4
Donor_3_after_6dpi = os.system("wc -l SRR5660045_mapped_1.fq")
Donor_3_reads_after_6 = (int(Donor_3_after_6dpi))/4

#Write the number of read pairs to the log file before and after Bowtie2 mapping/filtering 

logfile.write("Donor 1 (2dpi) had " + str(Donor_1_reads_before) + " read pairs before Bowtie2 filtering and " + str(Donor_1_reads_after) + " read pairs after.")
logfile.write("Donor 1 (6dpi) had " + str(Donor_1_reads_before_6) + " read pairs before Bowtie2 filtering and " + str(Donor_1_reads_after_6) + " read pairs after.")
logfile.write("Donor 3 (2dpi) had " + str(Donor_3_reads_before_2) + " read pairs before Bowtie2 filtering and " + str(Donor_3_reads_after_2) + " read pairs after.")
logfile.write("Donor 3 (6dpi) had " + str(Donor_3_reads_before_6) + " read pairs before Bowtie2 filtering and " + str(Donor_3_reads_after_6) + " read pairs after.")

#STEP 3#
#Assemble all 4 transcriptomes to produce 1 assembly via Spades
os.system("spades.py -k 77,99,127 -t 4 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SRR5660030_mapped_2.fq --pe-1 2 SRR5660033_mapped_1.fq --pe-2 2 SRR5660033_mapped_2.fq --pe-1 3 SRR5660044_mapped_1.fq --pe-2 3 SRR5660044_mapped_2.fq --pe-1 4 SRR5660045_mapped_1.fq --pe-2 4 SRR5660045_mapped_2.fq -o HCMV_transcriptome_assembly/")

#STEP 4#
#Change directory to the HCMV transcriptome assembly
os.system("cd ~/PipelineProject_Anna_Bielanski/HCMV_transcriptome_assembly")

#Open the fasta w/ the contigs
handle = open("contigs.fasta")

contig_count = 0   #Set contig count to 0
contig_list = []   #Open up a contig list to store the contigs in 
for record in SeqIO.parse(handle, format = "fasta"):   #Set up for loop to parse the fasta file of the contigs
    len_record = len(record)   #Retrieve length of contigs w/ len(record)
    if len_record > 1000:      #If the length > 1000 bp       
        contig_list.append(int(len_record))  #Append the len_record to contig_list
        contig_count += 1               #Add 1 to counter

bp_assembly = sum(contig_list)

#Calculate number of contigs with a length > 1000bp
logfile.write("There are " + str(contig_count) + " contigs > 1000 bp in the assembly.")


#Calculate length of the assembly(the total number of bp in all of the contigs > 1000 bp in length)
logfile.write("There are " + str(bp_assembly) + " bp in the assembly.")


#STEP 5#

#Retrieve longest contig from Spades assembly
max_len_contig = 0
max_contig =""
for record in SeqIO.parse(handle, format = "fasta"):
    if len(record) > max_len_contig:
        max_len_contig = len(record)
        max_contig = record.seq

#Need to retrieve fasta files of subfamily since unable to transfer fasta file to terminal
#Use esearch to look for record sequences in the nucleotide database that belong to the subfamily
#Use Entrez.read to read the records retrieved using esearch in handle1
#Save each record's id in id_list
#Turn the id_list into a comma separated string, id_string
#Use id=id_string in efetch to retrieve the sequeces of the records present in that id_string in FASTA format
#Open a fasta file to write & store all fasta files retrieved by efetch
Entrez.email = "abielanski@luc.edu"
handle1 = Entrez.esearch(db="nucleotide", term="Betaherpesvirinae")
record = Entrez.read(handle1)
id_list = record["IdList"]
id_string = ','.join(id_list)
handle2 = Entrez.efetch(db="nucleotide", id=id_string, rettype="fasta", retmode="text")
output_fas = open("betaherpesvirinae.txt", "w")
output_fas.write(handle2)

#Make a local database of sequences present in the betaherpesvirinae subfamily
os.system("makeblastdb -in betaherpesvirinae.fasta -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl")


#input file is the max_contig from Spades assembly
#Open an output file(virus_results) to store BLAST results
#Generate a blast command to blast query sequence against the generated betaherpesvirinae database 
#Create a hits command to pull the top 10 hits & write them into the log file
output_file = "virus_results.csv"
blast_command = "blastn -query "+ max_contig +" -db betaherpesvirinae -out "+ output_file +" -outfmt 6 sacc pident length qstart qend sstart send bitscore evalue stitle -max_hsps 1"
os.system("blast_command")
hits_command = "head -n 10 virus_results.csv | -out" + logfile
print(os.system("hits_command"))


