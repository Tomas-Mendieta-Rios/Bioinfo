import os
import sys
import shutil
import subprocess
from Bio import SeqIO
import argparse



def redirectLog(archivo_log, string):

    print(string)

    with open(archivo_log, 'a') as logueo:
        logueo.write(string + '\n')



def remove_file_if_exists(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)



def check_result(command_result, step_name):
    if command_result is None:
        redirectLog(error_log_file, f"Error al ejecutar {step_name}.")
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
        result = subprocess.run(command, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            redirectLog(error_log_file, f"Error ejecutando el comando: {command}\n{result.stderr}")
            return None
        return result.stdout
    except Exception as e:
        redirectLog(error_log_file, f"Excepción al ejecutar el comando: {command}\n{e}")
        return None





if __name__ == "__main__":

   
   log_file = os.environ.get('LOGFILE')
   error_log_file = os.environ.get('ERRORLOGFILE')

   parser = argparse.ArgumentParser(description = "Encuentra ORFs y motivos.")
   parser.add_argument("archivo_Genbank_entrada", help = "La opsina mutada en formato Genbank")
   parser.add_argument("archivo_fasta_salida", help = "Archivo FASTA con los ORFs")
   parser.add_argument("archivo_txt_salida", help = "Archivo TXT con los motivos") 


   args = parser.parse_args()




   # Remover archivos antiguos si existen
   remove_file_if_exists("motifs_final_results.txt")
   remove_file_if_exists("orf_proteins.fasta")


   # Leer el archivo GenBank
   GenBankPath = args.archivo_Genbank_entrada
   try:
       GenBank = SeqIO.read(GenBankPath, "genbank")
       nucleotide_sequence = str(GenBank.seq)
   except FileNotFoundError:
       redirectLog(error_log_file, f"Archivo {GenBankPath} no encontrado.")
       sys.exit(1)
   except ValueError:
       redirectLog(error_log_file, "Error al leer el archivo GenBank. Verifique el formato.")
       sys.exit(1)


   # Guardar la secuencia en formato FASTA
   fasta_path = 'secuencia.fasta'
   with open(fasta_path, 'w') as fasta_file:
       fasta_file.write(f">Opsina_sequence\n{nucleotide_sequence}\n")


   # Ejecutar getorf
   orf_output = args.archivo_fasta_salida
   remove_file_if_exists(orf_output)
   redirectLog(log_file, "Buscando ORFs en la secuencia de nucleótidos...")
   getorf_command = f"getorf -sequence {fasta_path} -outseq {orf_output} -find {1}"
   getorf_result = run_command(getorf_command)

   # Verificar resultado de getorf
   if not os.path.exists(orf_output) or os.path.getsize(orf_output) == 0:
       redirectLog(error_log_file, "No se encontraron ORFs en la secuencia de nucleótidos.")
       sys.exit(1)

   redirectLog(log_file, f"ORFs encontrados y guardados en {orf_output}")
   orf_count = sum(1 for _ in SeqIO.parse(orf_output, "fasta"))
   redirectLog(log_file, f"Cantidad de ORFs encontrados: {orf_count}")


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
                       redirectLog(log_file, f"Motifs encontrado en la secuencia {i}")
                       redirectLog(log_file, f"HitCount: {hit_count}")
                       

       if hit_count != "0":

           # Leer el contenido del archivo fuente
           with open("motifs_individual_results.txt", 'r') as archivo_fuente:
               contenido = archivo_fuente.read()

           # Abrir el archivo destino en modo 'append' y agregar el contenido
           with open(args.archivo_txt_salida, 'a') as archivo_destino:
               archivo_destino.write(contenido)
   
   


   remove_file_if_exists("motifs_individual_results.txt")
   remove_file_if_exists("orf_individual.fasta")
   remove_file_if_exists("orf_copia.fasta")
   remove_file_if_exists("secuencia.fasta")
