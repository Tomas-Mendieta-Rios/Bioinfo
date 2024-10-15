# ejercicio_2.py
from Bio import SeqIO
from Bio.Blast import NCBIWWW

def main(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Verifica que hay al menos 6 secuencias en el archivo fasta
    if len(records) < 6:
        print("El archivo FASTA debe contener al menos 6 secuencias.")
        return

    # Iterar sobre cada secuencia en el archivo y distribuir en 6 archivos diferentes
    for idx, record in enumerate(records):
        # Realizar una búsqueda BLAST en la base de datos 'nr' para cada secuencia
        result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)

        # Determinar el archivo de salida en función del índice
        output_file = f"resultado_blast_{(idx % 6) + 1}.xml"

        # Guardar los resultados en el archivo correspondiente
        with open(output_file, "a") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        print(f"¡Búsqueda BLAST completada para {record.id} en {output_file}!")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Uso: python3 ejercicio_2.py <archivo_fasta>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    main(fasta_file)

