from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable


def generar_cadena_negativa(arnm_seq_positiva):
	# Lista inicial vacía
	largo = len(arnm_seq_positiva)
	arnm_seq_negativa = [''] * largo
	
	# Generación de cadena ARN complementaria
	for i in range(largo):
		nucleotido = arnm_seq_positiva[i]
		complementario = ''

		if nucleotido == 'A':
			complementario += 'T'
		if nucleotido == 'T':
			complementario += 'A'
		if nucleotido == 'C':
			complementario += 'G'
		if nucleotido == 'G':
			complementario += 'C'
		
		# La lista se encuentra invertida, para que la traducción 5' 3' se de desde el comienzo al final de la lista
		arnm_seq_negativa[largo - 1 - i] += complementario

	# Conversión de la lista de caracteres a cadena
	secuencia_str = ''.join(arnm_seq_negativa)

	# Crear un objeto Seq
	secuencia_final = Seq(secuencia_str)

	return secuencia_final


def traducir_en_tres_marcos(arnm_seq):
	# Traduccion en los tres marcos de lectura
	traducciones = []
	
	# Diccionario de traduccion. Cada key es un codon y los codones de stop son keys nulas
	dic_codon = CodonTable.standard_dna_table.forward_table
	


	for marco in range(3):
		traduccion = ""
		seq = arnm_seq[marco:]  # Desplazar el marco de lectura
		for i in range(0, len(seq) -2, 3):
			codon = seq[i:i+3]
			if len(codon) < 3:
				break

			amino_acido = dic_codon.get(codon) # Aqui se da la traduccion
			if len(traduccion) == 0:
				if amino_acido == "M": # Asegurar que el primer aminoacido sea metionina
					traduccion += amino_acido
			
			else:
				if amino_acido != None: # Detener la traduccion al llegar al codon de Stop
					traduccion += amino_acido
				else:
					break	
				
		traducciones.append(traduccion)


	return traducciones


def seleccionar_traduccion_mas_larga(traducciones):
      	# Encontrar la traduccion con la mayor longitud (criterio heuristico)
      	traduccion_mas_larga = max(traducciones, key = len)
	posicion_traduccion = traducciones.index(traduccion_mas_larga)
      	return posicion_traduccion, traduccion_mas_larga



def traducir_arnm_a_aminoacidos(archivo_genbank, archivo_fasta_salida):
      	# Leer la secuencia de ARN mensajero del archivo GenBank
      	record = SeqIO.read(archivo_genbank, "genbank")
	cadena_positiva = record.seq


	# Generar la cadena de ARN mensajero complementaria a la del genbank
	cadena_negativa = generar_cadena_negativa(cadena_positiva)
  

      	# Traducir la secuencia positiva en los tres marcos de lectura
	traduccion_positiva = traducir_en_tres_marcos(cadena_positiva)


	# Traducir la secuencia negativa en los tres marcos de lectura
	traduccion_negativa = traducir_en_tres_marcos(cadena_negativa)


	# Combinar todas las traducciones en una sola lista
	traducciones = traduccion_positiva + traduccion_negativa


      	# Seleccionar la traduccion mas larga
	posicion_traduccion, traduccion_mas_larga = seleccionar_traduccion_mas_larga(traducciones)


      	# Guardar la secuencia de aminoacidos en formato FASTA
	record_fasta = SeqRecord(Seq(traduccion_mas_larga), id = f"NP_000504." + str(posicion_traduccion), description = "medium-wave-sensitive opsin 1 [Homo sapiens]")


      	with open(archivo_fasta_salida, 'w') as archivo_fasta:
          	SeqIO.write(record_fasta, archivo_fasta, "fasta")


if __name__ == "__main__":

   import argparse


   parser = argparse.ArgumentParser(description = "Traducir una secuencia de ARNm y guardar la cadena de mayor longitud en formatoFASTA.")
   parser.add_argument("archivo_genbank", help = "Archivo GenBank a traducir")
   parser.add_argument("archivo_fasta_salida", help = "Archivo FASTA donde se guardara la secuencia de aminoacidos") 


   args = parser.parse_args()


   traducir_arnm_a_aminoacidos(args.archivo_genbank, args.archivo_fasta_salida)








