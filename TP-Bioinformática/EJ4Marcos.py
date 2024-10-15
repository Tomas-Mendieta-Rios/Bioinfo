# -*- coding: utf-8 -*-
import sys
import subprocess
from Bio import SeqIO
import os

GenBankPath = "opsina.gb"
GenBank = SeqIO.read(GenBankPath, "genbank")
nucleotide_sequence = str(GenBank.seq)
print(nucleotide_sequence[:10])

fasta_path = 'secuencia.fasta'
with open(fasta_path, 'w') as fasta_file:
    fasta_file.write(f">Opsina_sequence\n{nucleotide_sequence}\n")

prositeData_path = 'prosite.dat'
with open(prositeData_path, 'r') as file:
    prosite_db = file.read()

orf_output = "orf_proteins.fasta"
motifs_output = "motifs_results.txt"

def run_command(command):
    """Ejecuta un comando en la terminal y retorna la salida"""
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"Error ejecutando el comando: {command}\n{result.stderr}")
        return None
    return result.stdout

"""##ORF"""

if os.path.exists(orf_output):
    os.remove(orf_output)
# Ejecuta getorf para encontrar los ORFs
print("Encontrando ORFs en la secuencia de nucleótidos...")
getorf_command = f"getorf -sequence {fasta_path} -outseq {orf_output}"
getorf_result = run_command(getorf_command)

# Verifica si getorf se ejecutó correctamente
if getorf_result is None:
    print("Error al ejecutar getorf.")

# Verifica si se encontraron ORFs
if not os.path.exists(orf_output) or os.path.getsize(orf_output) == 0:
    print("No se encontraron ORFs en la secuencia de nucleótidos.")
else:
    print(f"ORFs encontrados y guardados en {orf_output}")
    orf_count = sum(1 for _ in SeqIO.parse(orf_output, "fasta"))
    print(f"Cantidad de ORFs encontrados: {orf_count}")

orf_count = sum(1 for _ in SeqIO.parse(orf_output, "fasta"))
print(f"Cantidad de ORFs encontrados: {orf_count}")

# Llamar a EMBOSS 'geecee' para calcular el contenido de GC
os.system('geecee -sequence secuencia.fasta -outfile resultado_gc.txt')

# Mostrar el contenido del archivo de resultados
with open('resultado_gc.txt', 'r') as result_file:
    print(result_file.read())

#procesamiento manual
Citocina = nucleotide_sequence.count("C")
Adenina = nucleotide_sequence.count("A")
Timina = nucleotide_sequence.count("T")
Guanina = nucleotide_sequence.count("G")
print("Cantidad de Citocina:", Citocina)
print("Cantidad de Adenina:", Adenina)
print("Cantidad de Timina:", Timina)
print("Cantidad de Guanina:", Guanina)
Total = Citocina + Adenina + Timina + Guanina
print("Total de nucleotidos:", Total)
print("Contenido de G o C:", (Citocina + Guanina) / Total * 100, "%")

# Contar la cantidad de veces que aparece la secuencia "GCA" en la secuencia de nucleótidos
ATG_count = nucleotide_sequence.count("ATG")

print(f'La secuencia "ATG" aparece {ATG_count} veces en la cadena de nucleótidos.')

# Contar la cantidad de veces que aparece la secuencia "GCA" en la secuencia de nucleótidos
stop_count = nucleotide_sequence.count("TGA")+nucleotide_sequence.count("TAA")+nucleotide_sequence.count("TAG")

print(f'La secuencia "stop" aparece {stop_count} veces en la cadena de nucleótidos.')

"""##Motif"""

# Ejecuta patmatmotifs para analizar los motivos en las secuencias de proteínas
print("Analizando motivos en las secuencias de proteínas obtenidas...")
patmatmotifs_command = f"patmatmotifs -sequence {orf_output} -db {prosite_db} -outfile {motifs_output}"
patmatmotifs_result = run_command(patmatmotifs_command)

# Verifica si patmatmotifs se ejecutó correctamente
if patmatmotifs_result is None:
    print("Error al ejecutar patmatmotifs.")
else:
    print(f"Análisis de motivos completado. Resultados guardados en {motifs_output}")

print(f"Análisis de motivos completado. Resultados guardados en {motifs_output}")