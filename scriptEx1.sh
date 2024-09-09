# Variables de entrada y salida
ARCHIVO_GENBANK="opsina.gb"
ARCHIVO_FASTA_SALIDA="traduccion.fasta"

# Comprobar si el archivo GenBank existe
if [[ ! -f "$ARCHIVO_GENBANK" ]]; then
	echo "Error: El archivo GenBank '$ARCHIVO_GENBANK' no existe."
	exit 1
fi


# Ejecutar el script Python
python3 Ex1.py "$ARCHIVO_GENBANK" "$ARCHIVO_FASTA_SALIDA"


# Verificar si el script Python se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script Python fallo."
	exit 1
fi


echo "La traduccion se completo y se guardo en '$ARCHIVO_FASTA_SALIDA'."
