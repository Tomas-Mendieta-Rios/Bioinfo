from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Leer la secuencia en formato fasta
fasta_sequence = open("tu_archivo.fasta").read()

# Ejecutar BLAST de manera remota
result_handle = NCBIWWW.qblast("blastp", "nr", fasta_sequence)

# Guardar el resultado en un archivo XML
with open("blast_output.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

# Cerrar el resultado
result_handle.close()

# Leer el archivo XML y parsear los resultados
with open("blast_output.xml") as result_file:
    blast_record = NCBIXML.read(result_file)

# Interpretar los resultados (ejemplo)
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < 0.01:
            print(f"***Aligment***")
            print(f"sequence: {alignment.title}")
            print(f"length: {alignment.length}")
            print(f"E-value: {hsp.expect}")
