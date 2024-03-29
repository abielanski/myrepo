Pipeline Project Track 2

The tools that need to be installed or accessible to run this pipeline are Python 3.11, Bowtie2, SPAdes, and BLAST+. The sample test data uploaded to repo are paired-end fastq files after the fastq-dump command. This sample test data can be used starting with step 2. Please use SampleTestData.

This pipeline should run in the terminal by typing the command $python PipelineProject.py

Step 2 runs the Accession_to_FASTA.py code to retrieve the HCMV genome from NCBI and converts it to FASTA format. Next, it creates an index for HCMV(NCBI accession NC_006273.2) using bowtie2-build command, which should output a set of 6 files. After that, the reads will be mapped and saved to the generated HCMV index using bowtie2 --quiet command. This command is used for all the paired end fastq files. 
Finally, the number of reads in each transcriptome (1 set of the paired-end fastq files) before and after Bowtie2 filtering and mapping will be written into the log file generated in step 1. The wc -l command is used to figure out the number of lines in each (or one of) fastq file. This .py code should take the number of lines in each fastq and divide it by 4 to retrieve the read pair count. This is done on the paired end fastq files and the mapped fastq files. 

Step 3 takes the 4 paired mapped fastq files and assembles them into 1 assembly using SPAdes. The spades.py -k 77,99,127 -t 4 --only-assembler commands should complete this task. The command also includes --pe- since there are several pairs of the mapped fastq files. The assembly and additional files generated are saved in the HCMV_transcriptome_assembly folder. 

Step 4 calculates the number of contigs with a length greater than 1000bp and also calculates the bp length of the assembly using the total number of bp in all of the contigs greater than 1000 bp in length. The number of contigs and the length of the assembly are found using the Python code embedded under this step, and the outputs are written into the the log file. 

Step 5 blasts the query sequence against the generated betaherpesvirinae database and a hits command is created to pull the top 10 hits & write them into the log file.

