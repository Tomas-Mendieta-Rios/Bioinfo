import os
import sys
import shutil
import subprocess
from Bio import SeqIO

def remove_file_if_exists(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)

def check_result(command_result, step_name):
    if command_result is None:
        print(f"Error al ejecutar {step_name}.")
        sys.exit(1)

def individual_orf(archivo_orfs_path):
  if os.path.exists("orf_individual.fasta"):
    # Si el archivo existe, eliminarlo
    os.remove("orf_individual.fasta")

  ruta_salida = 'orf_individual.fasta'
  # Leer el archivo de entrada
  with open(archivo_orfs_path, 'r') as archivo_entrada:
      lineas = archivo_entrada.readlines()

  # Lista para almacenar las líneas del primer ORF
  lineas_a_borrar = []
  lineas_restantes = []

  # Variables de control
  esta_guardando_orf = False

  # Recorrer las líneas y clasificar entre las que se guardarán y las que se conservarán
  for linea in lineas:
      # Si la línea comienza con ">", indica el inicio de un ORF
      if linea.startswith('>'):
          # Si ya estábamos guardando el primer ORF, dejamos de guardar y pasamos al siguiente
          if esta_guardando_orf:
              esta_guardando_orf = False
          # Si no estábamos guardando, es el primer ORF que encontramos, así que empezamos a guardar
          elif not lineas_a_borrar:
              esta_guardando_orf = True

      # Guardar las líneas del primer ORF o del resto
      if esta_guardando_orf:
          lineas_a_borrar.append(linea)
      else:
          lineas_restantes.append(linea)

  # Guardar las líneas del primer ORF en el archivo de salida
  with open(ruta_salida, 'w') as archivo_salida:
      archivo_salida.writelines(lineas_a_borrar)

  # Sobrescribir el archivo original con las líneas restantes
  with open(archivo_orfs_path, 'w') as archivo_entrada:
      archivo_entrada.writelines(lineas_restantes)

def run_command(command):
    """Ejecuta un comando en la terminal y retorna la salida"""
    try:
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        if result.returncode != 0:
            print(f"Error ejecutando el comando: {command}\n{result.stderr}")
            return None
        return result.stdout
    except Exception as e:
        print(f"Excepción al ejecutar el comando: {command}\n{e}")
        return None


# Remover archivos antiguos si existen
remove_file_if_exists("motifs_final_results.txt")
remove_file_if_exists("orf_proteins.fasta")

# Leer el archivo GenBank
GenBankPath = "opsina.gb"
try:
    GenBank = SeqIO.read(GenBankPath, "genbank")
    nucleotide_sequence = str(GenBank.seq)
except FileNotFoundError:
    print(f"Archivo {GenBankPath} no encontrado.")
    sys.exit(1)
except ValueError:
    print("Error al leer el archivo GenBank. Verifique el formato.")
    sys.exit(1)

# Guardar la secuencia en formato FASTA
fasta_path = 'secuencia.fasta'
with open(fasta_path, 'w') as fasta_file:
    fasta_file.write(f">Opsina_sequence\n{nucleotide_sequence}\n")

# Ejecutar getorf
orf_output = "orf_proteins.fasta"
remove_file_if_exists(orf_output)
print("Encontrando ORFs en la secuencia de nucleótidos...")
getorf_command = f"getorf -sequence {fasta_path} -outseq {orf_output} -find {1}"
getorf_result = run_command(getorf_command)

# Verificar resultado de getorf
if not os.path.exists(orf_output) or os.path.getsize(orf_output) == 0:
    print("No se encontraron ORFs en la secuencia de nucleótidos.")
    sys.exit(1)

print(f"ORFs encontrados y guardados en {orf_output}")
orf_count = sum(1 for _ in SeqIO.parse(orf_output, "fasta"))
print(f"Cantidad de ORFs encontrados: {orf_count}")

# Procesar motivos y guardar resultados
"""
#Estas lineas solo correr cuando se pone a punto el patmatmotifs
prosite_path = "/PROSITE"
prosextract_command = f"prosextract -prositedir {prosite_path}"
check_result(run_command(prosextract_command), "prosextract")
"""
shutil.copy(orf_output, 'orf_copia.fasta')

for i in range(1, orf_count + 1):
    individual_orf('orf_copia.fasta')

    # Ejecutar patmatmotifs
    patmatmotifs_command = f"patmatmotifs -sequence orf_individual.fasta -outfile motifs_individual_results.txt"
    check_result(run_command(patmatmotifs_command), "patmatmotifs")

    # Procesar resultados de motivos
    with open("motifs_individual_results.txt", 'r') as archivo:
        for linea in archivo:
            if 'HitCount:' in linea:
                hit_count = linea.split(':')[-1].strip()
                if hit_count != "0":
                    print(f"Motifs encontrado en la secuencia {i}\n")
                    print(f'HitCount: {hit_count}\n')
                    shutil.copyfileobj(archivo, open("motifs_final_results.txt", 'a'))

remove_file_if_exists("motifs_individual_results.txt")
remove_file_if_exists("orf_individual.fasta")
remove_file_if_exists("orf_copia.fasta")
remove_file_if_exists("secuencia.fasta")
