# Instalar Biopython y MUSCLE (sudo apt-get install muscle)

import subprocess
from Bio import SeqIO
import os
import argparse


def MSA(input_file_1, input_file_2, output_file):


   # Leer la secuencia query del primer archivo FASTA
   query_seq = None
   with open(input_file_1, "r") as f1:
        query_seq = list(SeqIO.parse(f1, "fasta"))

   
   # Leer las 10 primeras secuencias del segundo archivo FASTA
   blast_seqs = []
   with open(input_file_2, "r") as f2:
        for record in SeqIO.parse(f2, "fasta"):
            blast_seqs.append(record)
            if len(blast_seqs) == 10:
                break


   # Crear un archivo intermedio para las secuencias a alinear
   temp_file = "temp_msa.fasta"
   with open(temp_file, "w") as temp_f:
        SeqIO.write(query_seq, temp_f, "fasta")
        SeqIO.write(blast_seqs, temp_f, "fasta")


   # Ejecutar MUSCLE y redirigir la salida a un archivo de log
   log_file = "muscle_output.log"
   with open(log_file, "w") as log:

      try:
           subprocess.run(["muscle", "-in", temp_file, "-out", output_file], stdout=log, stderr=log, check=True)

      except subprocess.CalledProcessError as e:
           print(f"Error ejecutando MUSCLE: {e}")


   # Limpiar el archivo temporal
   os.remove(temp_file)



if __name__ == "__main__":

   
   parser = argparse.ArgumentParser(description = "Realizar un MSA y guardar el resultado en un archivo FASTA.")
   parser.add_argument("archivo_fasta_entrada_1", help = "Archivo FASTA con la secuencia query")
   parser.add_argument("archivo_fasta_entrada_2", help = "Archivo FASTA con los resultados BLAST")
   parser.add_argument("archivo_fasta_salida", help = "Archivo FASTA alineado") 


   args = parser.parse_args()


   MSA(args.archivo_fasta_entrada_1, args.archivo_fasta_entrada_2, args.archivo_fasta_salida)







