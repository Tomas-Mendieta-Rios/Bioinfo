from Bio import SeqIO
import requests
import argparse
import json
import os

parser = argparse.ArgumentParser(description='Process a GenBank file and generate reading frames.')
parser.add_argument('ORF_path', type=str, help='Path to the O file')
args = parser.parse_args()
ORF_path = args.ORF_path


# Ruta al archivo FASTA
fasta_file = ORF_path


# Directorio de salida para guardar los resultados
output_dir = "prosite_results"

# Crear el directorio de salida si no existe
os.makedirs(output_dir, exist_ok=True)

# Endpoint de la API de PROSITE
prosite_url = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"

# Leer el archivo FASTA y procesar cada secuencia
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = str(record.seq)
        # Parámetros para la API de PROSITE
        params = {
            'seq': sequence,
            'output': 'json'
        }
        # Realizar la solicitud GET
        response = requests.get(prosite_url, params=params)
        if response.status_code == 200:
            data = response.json()
            # Construir la estructura de salida en el formato esperado
            results = {
                'n_match': data['n_match'],
                'n_seq': data['n_seq'],
                'matchset': []
            }
            for hit in data['matchset']:
                match_info = {
                    'sequence_ac': hit['sequence_ac'],
                    'start': hit['start'],
                    'stop': hit['stop'],
                    'signature_ac': hit['signature_ac'],
                    'score': hit['score'],
                    'level': hit['level']
                }
                results['matchset'].append(match_info)
            
            # Nombre del archivo de salida basado en el ID de la secuencia
            output_file = os.path.join(output_dir, f"{record.id}_prosite_results.json")
            # Guardar los resultados en un archivo JSON
            with open(output_file, "w") as output_handle:
                json.dump(results, output_handle, indent=2)
            
            print(f"Resultados guardados para la secuencia {record.id} en {output_file}")
        else:
            print(f"Error al procesar la secuencia {record.id}: {response.status_code}")
