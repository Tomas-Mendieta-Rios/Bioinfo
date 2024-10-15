#pip install primer3-py esto va en el bash
import primer3
import argparse
from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW

from Bio.Blast.Applications import NcbiblastpCommandline
import re #re es de regular expresion. Permite trabajar f cilmente con strings y caracteres
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Process a GenBank file and generate reading frames.')
parser.add_argument('genbank_path', type=str, help='Path to the GenBank file')
args = parser.parse_args()
genbank_path = args.genbank_path
genbank = SeqIO.read(genbank_path, 'genbank')
mARN=genbank.seq

sequence_args = {
    'SEQUENCE_ID': 'example',
    'SEQUENCE_TEMPLATE':genbank.seq,
    'SEQUENCE_INCLUDED_REGION': [0, len(genbank.seq)]
}

default_params = {
    'SEQUENCE_ID': 'example',
    'SEQUENCE_TEMPLATE': genbank.seq,
    'SEQUENCE_INCLUDED_REGION': [0, len(genbank.seq)],
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MAX_TM': 67.0,
    'PRIMER_MIN_GC': 50.0,
    'PRIMER_MAX_GC': 60.0,
    'PRIMER_NUM_RETURN': 5
}

# Preguntar al usuario si quiere utilizar los valores por defecto
use_default = input("¿Desea utilizar los valores por defecto? (yes/no): ").strip().lower()

# Si el usuario no quiere los valores por defecto, permitir ingresar valores personalizados
if use_default != 'yes':
	primer3_params = {
	    'SEQUENCE_ID': 'example',
	    'SEQUENCE_TEMPLATE': genbank.seq,
	    'SEQUENCE_INCLUDED_REGION': [0, len(genbank.seq)],
	    'PRIMER_OPT_SIZE': int(input("Ingrese el tamaño óptimo del primer: ")),
	    'PRIMER_PICK_LEFT_PRIMER': int(input("Seleccione si desea escoger el primer izquierdo: ")),
	    'PRIMER_PICK_INTERNAL_OLIGO': int(input("Seleccione si desea escoger el oligo interno: ")),
	    'PRIMER_PICK_RIGHT_PRIMER': int(input("Seleccione si desea escoger el primer derecho: ")),
	    'PRIMER_MIN_SIZE': int(input("Ingrese el tamaño mínimo del primer: ")),
	    'PRIMER_MAX_SIZE': int(input("Ingrese el tamaño máximo del primer: ")),
	    'PRIMER_OPT_TM': float(input("Ingrese la temperatura óptima del primer: ")),
	    'PRIMER_MIN_TM': float(input("Ingrese la temperatura mínima del primer: ")),
	    'PRIMER_MAX_TM': float(input("Ingrese la temperatura máxima del primer: ")),
	    'PRIMER_MIN_GC': float(input("Ingrese el contenido mínimo de GC del primer: ")),
	    'PRIMER_MAX_GC': float(input("Ingrese el contenido máximo de GC del primer: ")),
	    'PRIMER_NUM_RETURN': int(input("Ingrese el número de primers a retornar: "))
	}
else:
    primer3_params = default_params

primer_results = primer3.bindings.design_primers(sequence_args, primer3_params)

# Lista para almacenar los registros SeqRecord
records = []

# Recorrer los primers diseñados
for i in range(primer_results['PRIMER_PAIR_NUM_RETURNED']):
    left_primer = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
    right_primer = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']

    # Crear objetos SeqRecord para los primers
    left_record = SeqRecord(Seq(left_primer), id=f"Primer_{i+1}_left", description="Left primer")
    right_record = SeqRecord(Seq(right_primer), id=f"Primer_{i+1}_right", description="Right primer")

    # Agregar los registros a la lista
    records.append(left_record)
    records.append(right_record)
    
# Guardar los registros en un archivo FASTA
output_file = "primers.fasta"
with open(output_file, "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")

print(f"Se han guardado los primers en el archivo {output_file}.")
    

