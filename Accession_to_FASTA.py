#from Bio import SeqIO
from Bio import Entrez

Entrez.email="abielanski@luc.edu"

ref_genome = (Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="NC_006273.2").read())

with open("ref_genome_out.txt", "w") as outfile:
    print(ref_genome, file=outfile)
