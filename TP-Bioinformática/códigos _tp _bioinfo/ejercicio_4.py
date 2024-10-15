import sys
import subprocess
import os

def run_command(command):
    """Ejecuta un comando en la terminal y retorna la salida"""
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"Error ejecutando el comando: {command}\n{result.stderr}")
        return None
    return result.stdout

def main(sequence_file):
    # Archivos de entrada y salida
    orf_output = "orf_proteins.fasta"
    prosite_db = "prosite.dat"
    motifs_output = "motifs_results.txt"

    # Verifica que el archivo de secuencia de nucleótidos exista
    if not os.path.exists(sequence_file):
        print(f"El archivo {sequence_file} no existe.")
        return
    
    # Ejecuta getorf para encontrar los ORFs
    print("Encontrando ORFs en la secuencia de nucleótidos...")
    getorf_command = f"getorf -sequence {sequence_file} -outseq {orf_output}"
    getorf_result = run_command(getorf_command)
    
    # Verifica si getorf se ejecutó correctamente
    if getorf_result is None:
        print("Error al ejecutar getorf.")
        return

    # Verifica si se encontraron ORFs
    if not os.path.exists(orf_output) or os.path.getsize(orf_output) == 0:
        print("No se encontraron ORFs en la secuencia de nucleótidos.")
        return

    print(f"ORFs encontrados y guardados en {orf_output}")

    # Ejecuta patmatmotifs para analizar los motivos en las secuencias de proteínas
    print("Analizando motivos en las secuencias de proteínas obtenidas...")
    patmatmotifs_command = f"patmatmotifs -sequence {orf_output} -db {prosite_db} -outfile {motifs_output}"
    patmatmotifs_result = run_command(patmatmotifs_command)

    # Verifica si patmatmotifs se ejecutó correctamente
    if patmatmotifs_result is None:
        print("Error al ejecutar patmatmotifs.")
        return

    print(f"Análisis de motivos completado. Resultados guardados en {motifs_output}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python3 script.py <archivo_secuencia>")
        sys.exit(1)

    sequence_file = sys.argv[1]
    main(sequence_file)


