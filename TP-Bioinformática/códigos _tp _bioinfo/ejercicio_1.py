import argparse
from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
import re
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

def write_fasta(file_name, sequences):
    with open(file_name, 'a') as fasta_file:
        for name, sequence in sequences.items():
            fasta_file.write('>' + name + '\n')
            fasta_file.write(sequence + '\n')

def main():
    parser = argparse.ArgumentParser(description='Process a GenBank file and generate reading frames.')
    parser.add_argument('genbank_path', type=str, help='Path to the GenBank file')
    parser.add_argument('output_path', type=str, help='Output path for the FASTA file')

    args = parser.parse_args()
    genbank_path = args.genbank_path
    output_path = args.output_path

    genbank = SeqIO.read(genbank_path, 'genbank')
    mARN = genbank.seq

    marco_lectura = 0
    for i in range(3):
        seq_proteina = Seq(mARN[marco_lectura:]).translate()
        ORF = {f'ORF{marco_lectura+1}': str(seq_proteina)}
        write_fasta(output_path, ORF)
        marco_lectura += 1

    mARN_invertido = mARN[::-1]

    marco_lectura = 0
    for i in range(3):
        seq_proteina_invertida = Seq(mARN_invertido[marco_lectura:]).translate()
        ORF = {f'ORF{marco_lectura+4}': str(seq_proteina_invertida)}
        write_fasta(output_path, ORF)
        marco_lectura += 1

if __name__ == "__main__":
    main()

          
