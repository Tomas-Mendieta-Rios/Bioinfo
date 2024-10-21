import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def editar_transcripto(input_file, output_file):
    # Leer el archivo GenBank
    with open(input_file, "r") as file:
        records = list(SeqIO.parse(file, "genbank"))

    # Reemplazar T por C en la posición 689
    for record in records:
        if len(record.seq) >= 689:  # Asegurarse de que la posición existe
            seq_list = list(record.seq)  # Convertir a lista para modificar
            seq_list[688] = 'C'  # Reemplazar en la posición 689 (índice 688)
            record.seq = Seq(''.join(seq_list))  # Usar Seq para la nueva secuencia

    # Guardar los cambios en un nuevo archivo GenBank
    with open(output_file, "w") as file:
        SeqIO.write(records, file, "genbank")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Editar un transcripto en formato GenBank.")
    parser.add_argument("input_file", help="Archivo GenBank de entrada")
    parser.add_argument("output_file", help="Archivo GenBank de salida con el transcripto editado")

    args = parser.parse_args()

    editar_transcripto(args.input_file, args.output_file)

