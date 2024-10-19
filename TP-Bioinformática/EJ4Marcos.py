import os
import sys
import shutil
import subprocess
from Bio import SeqIO

if os.path.exists("motifs_final_results.txt"):
    os.remove("motifs_final_results.txt")
if os.path.exists("orf_proteins.fasta"):
    os.remove("orf_proteins.fasta")
    
GenBankPath = "opsina.gb"
GenBank = SeqIO.read(GenBankPath, "genbank")
nucleotide_sequence = str(GenBank.seq)

fasta_path = 'secuencia.fasta'
with open(fasta_path, 'w') as fasta_file:
    fasta_file.write(f">Opsina_sequence\n{nucleotide_sequence}\n")

def run_command(command):
    """Ejecuta un comando en la terminal y retorna la salida"""
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"Error ejecutando el comando: {command}\n{result.stderr}")
        return None
    return result.stdout

"""##ORF"""
orf_output = "orf_proteins.fasta"

if os.path.exists(orf_output):
    os.remove(orf_output)
# Ejecuta getorf para encontrar los ORFs
print("Encontrando ORFs en la secuencia de nucleótidos...")
getorf_command = f"getorf -sequence {fasta_path} -outseq {orf_output} -find {1}"
getorf_result = run_command(getorf_command)

# Verifica si getorf se ejecutó correctamente
if getorf_result is None:
    print("Error al ejecutar getorf.")

# Verifica si se encontraron ORFs
if not os.path.exists(orf_output) or os.path.getsize(orf_output) == 0:
    print("No se encontraron ORFs en la secuencia de nucleótidos.")
else:
    print(f"ORFs encontrados y guardados en {orf_output}")
    orf_count = sum(1 for _ in SeqIO.parse(orf_output, "fasta"))
    print(f"Cantidad de ORFs encontrados: {orf_count}")

# Contar la cantidad de veces que aparece la secuencia "GCA" en la secuencia de nucleótidos
ATG_count = nucleotide_sequence.count("ATG")
stop_count = nucleotide_sequence.count("TGA")+nucleotide_sequence.count("TAA")+nucleotide_sequence.count("TAG")

print(f'La secuencia "start" aparece {ATG_count} veces en la cadena de nucleótidos.')
print(f'La secuencia "stop" aparece {stop_count} veces en la cadena de nucleótidos.')

"""##Motif"""

motifs_output = "motifs_individual_results.txt"
motifs_final_output = "motifs_final_results.txt"
prosite_path = "/PROSITE"
ruta_copia = 'orf_copia.fasta'

shutil.copy(orf_output, ruta_copia)

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

prosextract_command = f"prosextract -prositedir {prosite_path}"
prosextract_result = run_command(prosextract_command)

i=1
while i<=orf_count:
  print("/n", i, "/n" )
  individual_orf('orf_copia.fasta')


  # Ejecuta patmatmotifs para analizar los motivos en las secuencias de proteínas
  print(f"Analizando motivos de la secuencia {i}")
  orf_output = "orf_individual.fasta"
  patmatmotifs_command = f"patmatmotifs -sequence {orf_output} -outfile {motifs_output} "
  patmatmotifs_result = run_command(patmatmotifs_command)
  # Verifica si patmatmotifs se ejecutó correctamente
  if patmatmotifs_result is None:
      print("Error al ejecutar patmatmotifs.")

  with open(motifs_output, 'r') as archivo:
      for linea in archivo:
          if 'HitCount:' in linea:
              # Extraer el valor del HitCount
              hit_count = linea.split(':')[-1].strip()
              print(f'HitCount: {hit_count}')
  if hit_count == "0":
    print("No se encontraron coincidencias.")
  else:

    # Leer el contenido del archivo fuente
    with open(motifs_output, 'r') as archivo_fuente:
        contenido = archivo_fuente.read()

    # Abrir el archivo destino en modo 'append' y agregar el contenido
    with open(motifs_final_output, 'a') as archivo_destino:
        archivo_destino.write(contenido)

    print(f'El contenido de {motifs_output} se ha agregado al final de {motifs_final_output}.')

  i=i+1

os.remove("motifs_individual_results.txt")
os.remove("orf_individual.fasta")
os.remove("orf_copia.fasta")
os.remove("secuencia.fasta")
