from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

# Lista de archivos de resultados BLAST
blast_result_files = [f"resultado_blast_{i}.xml" for i in range(1, 3)]

top_hits_sequences = []

# Iterar sobre los archivos XML de resultados del BLAST
for blast_result_file in blast_result_files:
    with open(blast_result_file, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments[:10]:
                for hsp in alignment.hsps:
                    sequence = SeqRecord(Seq(hsp.sbjct), id=alignment.title.split()[1])
                    top_hits_sequences.append(sequence)

# Ordenar las secuencias por longitud y seleccionar las 10 primeras
top_hits_sequences.sort(key=lambda x: len(x.seq), reverse=True)
top_hits_sequences = top_hits_sequences[:10]

# Guardar las secuencias en un archivo FASTA temporal
input_fasta = "temp_sequences.fasta"
SeqIO.write(top_hits_sequences, input_fasta, "fasta")

# Ejecutar Clustal Omega para alinear las secuencias
output_fasta = "alignment.fasta"
clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta, outfile=output_fasta, verbose=True, auto=True)
stdout, stderr = clustalomega_cline()

# Leer el alineamiento resultante
alignment = SeqIO.parse(output_fasta, "fasta")

# Convertir el alineamiento en un objeto MultipleSeqAlignment
alignment = MultipleSeqAlignment(list(alignment))

# Imprimir el alineamiento y guardar en archivo
print(alignment)
SeqIO.write(alignment, output_fasta, "fasta")

print(f"Alineamiento guardado en {output_fasta}")

