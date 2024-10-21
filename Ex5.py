# Es necesario tener instalado primer3-py (pip install primer3-py o python3 -m pip install primer3-py, dependiendo la versión de pip instalada)

import primer3
import argparse
from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.etree.ElementTree as ET
import os


def leer_parametros_xml(ruta_xml):
    """Lee parámetros desde un archivo XML."""
    if not os.path.exists(ruta_xml):
        raise FileNotFoundError(f"El archivo {ruta_xml} no existe.")
    
    try:
        tree = ET.parse(ruta_xml)
        root = tree.getroot()
        
        params = {
            'PRIMER_OPT_SIZE': int(root.find('PRIMER_OPT_SIZE').text),
            'PRIMER_PICK_LEFT_PRIMER': int(root.find('PRIMER_PICK_LEFT_PRIMER').text),
            'PRIMER_PICK_INTERNAL_OLIGO': int(root.find('PRIMER_PICK_INTERNAL_OLIGO').text),
            'PRIMER_PICK_RIGHT_PRIMER': int(root.find('PRIMER_PICK_RIGHT_PRIMER').text),
            'PRIMER_MIN_SIZE': int(root.find('PRIMER_MIN_SIZE').text),
            'PRIMER_MAX_SIZE': int(root.find('PRIMER_MAX_SIZE').text),
            'PRIMER_OPT_TM': float(root.find('PRIMER_OPT_TM').text),
            'PRIMER_MIN_TM': float(root.find('PRIMER_MIN_TM').text),
            'PRIMER_MAX_TM': float(root.find('PRIMER_MAX_TM').text),
            'PRIMER_MIN_GC': float(root.find('PRIMER_MIN_GC').text),
            'PRIMER_MAX_GC': float(root.find('PRIMER_MAX_GC').text),
            'PRIMER_NUM_RETURN': int(root.find('PRIMER_NUM_RETURN').text)
        }
        
        return params
    
    except ET.ParseError as e:
        raise ValueError(f"Error en el formato del archivo XML: {e}")


def guardar_parametros_xml(ruta_xml, params):
    """Guarda los parámetros en un archivo XML."""
    root = ET.Element('parametros')
    
    for key, value in params.items():
        param_element = ET.SubElement(root, key)
        param_element.text = str(value)
    
    tree = ET.ElementTree(root)
    tree.write(ruta_xml)



parser = argparse.ArgumentParser(description='Process a GenBank file and generate reading frames.')
parser.add_argument('archivo_genbank', help='Archivo GenBank a procesar')
parser.add_argument('archivo_parametros', help='Archivo XML con los parámetros de diseño')
parser.add_argument('archivo_primers', help='Archivo FASTA de salida con los primers')
args = parser.parse_args()

genbank_path = args.archivo_genbank
genbank = SeqIO.read(genbank_path, 'genbank')
transcripto = genbank.seq


sequence_args = {
    'SEQUENCE_ID': 'example',
    'SEQUENCE_TEMPLATE': str(transcripto),
    'SEQUENCE_INCLUDED_REGION': [0, len(transcripto)]
}


# Cargar parámetros desde el archivo XML
try:
    default_params = leer_parametros_xml(args.archivo_parametros)
    
except (FileNotFoundError, ValueError) as e:
    print(e)
    exit(1)


# Preguntar al usuario si quiere utilizar los valores por defecto
use_default = input("¿Desea utilizar los valores por defecto? (si/no): ").strip().lower()


# Si el usuario no quiere los valores por defecto, permitir ingresar valores personalizados
if use_default != 'si':
    primer3_params = default_params.copy()  # Hacer una copia de los parámetros por defecto
    primer3_params.update({
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
    })

    # Guardar los parámetros actualizados en el archivo XML
    guardar_parametros_xml(args.archivo_parametros, primer3_params)
else:
    primer3_params = default_params


# Generar los primers
primer_results = primer3.bindings.designPrimers(sequence_args, primer3_params)

# Lista para almacenar los registros SeqRecord
records = []

# Recorrer los primers diseñados
for i in range(primer_results['PRIMER_PAIR_NUM_RETURNED']):
    left_primer = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
    right_primer = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']

    # Crear objetos SeqRecord para los primers
    left_record = SeqRecord(Seq(left_primer), id=f"Primer_{i+1}_l", description="Primer Izquierdo")
    right_record = SeqRecord(Seq(right_primer), id=f"Primer_{i+1}_r", description="Primer Derecho")

    # Agregar los registros a la lista
    records.append(left_record)
    records.append(right_record)
    
# Guardar los registros en un archivo FASTA
with open(args.archivo_primers, "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")

print(f"Se han guardado los primers en el archivo {output_file}.")


