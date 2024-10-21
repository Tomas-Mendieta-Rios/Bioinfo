# Variables 
ARCHIVO_GENBANK_1="opsina.gb"
ARCHIVO_GENBANK_2="opsinamutada.gb"
ARCHIVO_FASTA_1="traduccion.fasta"
ARCHIVO_FASTA_2="secuencias_BLAST.fasta"
ARCHIVO_FASTA_3="MSA.fasta"
ARCHIVO_FASTA_5="primers.fasta"
SCRIPT_MUTADOR="mutador.py"
PARAMETROS_PRIMERS="param.xml"
SCRIPT_EJ_1="Ex1.py"
SCRIPT_EJ_2="Ex2.py"
SCRIPT_EJ_3="Ex3.py"
SCRIPT_EJ_4="Ex4.py"
SCRIPT_EJ_5="Ex5.py"



# Comprobar si el archivo GenBank existe
if [[ ! -f "$ARCHIVO_GENBANK_1" ]]; then
	echo "Error: El archivo GenBank '$ARCHIVO_GENBANK_1' no existe."
	exit 1
fi



# Ejectura el script que incorpora la mutación al transcripto
python3 "$SCRIPT_MUTADOR" "$ARCHIVO_GENBANK_1" "$ARCHIVO_GENBANK_2"


# Verificar si el script mutador se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_MUTADOR' fallo."
	exit 1
fi


echo "Se incorporó la mutación y se guardo en '$ARCHIVO_GENBANK_2'."




# Ejecutar el script Python para la traduccion
python3 "$SCRIPT_EJ_1" "$ARCHIVO_GENBANK_2" "$ARCHIVO_FASTA_1"


# Verificar si el script Python de la traduccion se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_1' fallo."
	exit 1
fi


echo "La traduccion se completo y se guardo en '$ARCHIVO_FASTA_1'."




# Ejecutar el script Python del alineamiento BLAST
python3 "$SCRIPT_EJ_1"


# Verificar si el script Python del BLAST se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_2' fallo."
	exit 1
fi




# Ejecutar el script Python del alineamiento múltiple
python3 "$SCRIPT_EJ_3" "$ARCHIVO_FASTA_1" "$ARCHIVO_FASTA_2" "$ARCHIVO_FASTA_3"


# Verificar si el script Python del BLAST se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_3' fallo."
	exit 1
fi


echo "El alineamiento múltiple se completo y se guardo en '$ARCHIVO_FASTA_3'."




# Ejecutar el script Python del generador de primers
python3 "$SCRIPT_EJ_5" "$ARCHIVO_GENBANK_2" "$PARAMETROS_PRIMERS" "$ARCHIVO_FASTA_5"


# Verificar si el script Python del BLAST se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_5' fallo."
	exit 1
fi


echo "El alineamiento múltiple se completo y se guardo en '$ARCHIVO_FASTA_5'."
