from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import argparse
import os


def redirectLog(archivo_log, string):

    print(string)

    with open(archivo_log, 'a') as logueo:
        logueo.write(string + '\n')



def run_blastp(query_file, blast_db, evalue=0.001, outfmt=5, output_file="blast_results.xml"):
    """
    Ejecuta una búsqueda BLASTP en el archivo de consulta contra la base de datos especificada y analiza los resultados.

    Parametros:
    - query_file (str): Ubicación del archivo FASTA que contiene la secuencia query.
    - blast_db (str): Ubicación de la base de datos BLAST (e.g., Swiss-Prot).
    - evalue (float): Umbral para reportar alineamientos.
    - outfmt (int): Formato de salida (5 es XML).
    - output_file (str): Archivo para guardar los resultados de BLASTP.
    
    Returns:
    - None: Prints parsed BLAST results.
    """
    
    # Configurar el comando BLASTP
    blastp_cline = NcbiblastpCommandline(query=query_file, db=blast_db, evalue=evalue, outfmt=outfmt, out=output_file)
    
    # Ejecutar el comando BLASTP
    if log_file:
      redirectLog(log_file, "Ejecutando búsqueda BLASTP...")

    stdout, stderr = blastp_cline()

    if log_file:
      redirectLog(log_file, f"Búsqueda BLAST completada. Resultados guardados en {output_file}")
    

def blast_fasta(blast_xml_file, output_fasta_file="ordered_blast_output.fasta"):
    """
    Analiza los resultados XML de BLAST, los ordena por mejor coincidencia (según el valor e) 
    y crea un archivo FASTA a partir de las secuencias alineadas.
    
    Parametros:
    - blast_xml_file (str): La ruta a los resultados de BLAST en formato XML.
    - output_fasta_file (str): El nombre del archivo FASTA de salida.
    
    Returns:
    - None: Writes the aligned sequences to a FASTA file ordered by best match.
    """
    results = []
    
    # Analizar los resultados de BLAST y recopilar los alineamientos junto a los e-values.
    with open(blast_xml_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        for blast_record in blast_records:
            query_id = blast_record.query
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    results.append({
                        'title': alignment.title,
                        'evalue': hsp.expect,  # E-value
                        'sequence': hsp.sbjct  # Subject sequence
                    })
    
    # Ordenar los resultados por e-value (de menor a mayor).
    sorted_results = sorted(results, key=lambda x: x['evalue'])
    
    # Incorporar los resultados ordenados en el archivo FASTA de salida.
    with open(output_fasta_file, "w") as fasta_file:
        for result in sorted_results:
            fasta_file.write(f">{result['title']} | e-value: {result['evalue']}\n")
            fasta_file.write(f"{result['sequence']}\n")
            fasta_file.write("\n")
    



if __name__ == "__main__":

   
   log_file = os.environ.get('LOGFILE')
   
   parser = argparse.ArgumentParser(description = "Realizar un BLAST y guardar el resultado en un archivo FASTA.")
   parser.add_argument("archivo_fasta_entrada_1", help = "Archivo FASTA con la secuencia query.")
   parser.add_argument("path_database", help = "Ubicación de la base de datos necesaria para realizar el BLAST.")
   parser.add_argument("archivo_fasta_salida", help = "Archivo FASTA con los resultados de BLAST.") 


   args = parser.parse_args()

   run_blastp(args.archivo_fasta_entrada_1, args.path_database)
   blast_fasta("blast_results.xml", output_fasta_file = args.archivo_fasta_salida)
   os.remove("blast_results.xml")
