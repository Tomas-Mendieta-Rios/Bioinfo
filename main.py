import ssl
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import urllib.request

# Deshabilitar la verificación de certificados SSL (no recomendado para producción)
ssl._create_default_https_context = ssl._create_unverified_context

# Secuencia en formato FASTA
fasta_sequence = """>sp|P01308|INS_HUMAN Insulin precursor OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFFVQPCSQITTCGILICSSLSTNQVEALYLVCGERGFFYTPKT
RRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"""

# Realiza el BLAST de manera remota en la base de datos de NCBI
result_handle = NCBIWWW.qblast("blastp", "nr", fasta_sequence)

# Guarda el resultado del BLAST en un archivo XML
with open("blast_remote_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

print("Resultado del BLAST guardado en 'blast_remote_result.xml'.")

# Cierra el handle de resultados
result_handle.close()

# Lee y parsea el archivo XML con los resultados
with open("blast_remote_result.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)

# Imprimir información básica del primer resultado
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print("****Align****")
        print(f"sequence: {alignment.title}")
        print(f"length: {alignment.length}")
        print(f"e-value: {hsp.expect}")
        print(f"score: {hsp.score}")
