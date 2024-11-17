# Entradas y salidas 
ARCHIVO_GENBANK_1="opsina.gb"
ARCHIVO_GENBANK_2="opsinamutada.gb"
ARCHIVO_FASTA_1="traduccion.fasta"
ARCHIVO_FASTA_2="secuencias_BLAST.fasta"
DATABASE_PATH_BLAST="$BLASTDB_PATH"
ARCHIVO_FASTA_3="MSA.fasta"
ARCHIVO_FASTA_4="orf_proteins.fasta"
ARCHIVO_TXT_SALIDA="motifs_final_results.txt"
ARCHIVO_FASTA_5="primers.fasta"


# Scripts
SCRIPT_MUTADOR="mutador.py"
PARAMETROS_PRIMERS="param.xml"
SCRIPT_EJ_1="Ex1.py"
SCRIPT_EJ_2="Ex2.py"
SCRIPT_EJ_3="Ex3.py"
SCRIPT_EJ_4="Ex4.py"
SCRIPT_EJ_5="Ex5.py"


# Logs
LOGFILE="Logs$(date +%Y-%m-%d).log"
ERRORLOGFILE="Errorlogs$(date +%Y-%m-%d).log"
LOGMUSCLE="muscle_output.log"


# Carpetas
LOGFOLDER="TP Logs"
RESULTSFOLDER="TP Resultados"


touch "$LOGFILE"
touch "$ERRORLOGFILE"


# Exportar ambos logs para que los scripts puedan acceder a los mismos
export LOGFILE
export ERRORLOGFILE




# Verificar si la carpeta TPLogs existe, si no, crearla
if [ ! -d "$LOGFOLDER" ]; then
	mkdir "$LOGFOLDER"
	echo "Carpeta '$LOGFOLDER' creada." | tee -a "$LOGFILE"
fi




# Verificar si la carpeta TPResultados existe, si no, crearla
if [ ! -d "$RESULTSFOLDER" ]; then
	mkdir "$RESULTSFOLDER"
	echo "Carpeta '$RESULTSFOLDER' creada." | tee -a "$LOGFILE"
fi




# Comprobar si el archivo GenBank existe
if [[ ! -f "$ARCHIVO_GENBANK_1" ]]; then
	echo "Error: El archivo GenBank '$ARCHIVO_GENBANK_1' no existe." | tee -a "$ERRORLOGFILE"
	exit 1
fi




# Ejectura el script que incorpora la mutación al transcripto
python3 "$SCRIPT_MUTADOR" "$ARCHIVO_GENBANK_1" "$ARCHIVO_GENBANK_2"


# Verificar si el script mutador se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_MUTADOR' fallo." | tee -a "$ERRORLOGFILE"
	exit 1
fi


echo "Se incorporó la mutación y se guardo en '$ARCHIVO_GENBANK_2'." | tee -a "$LOGFILE"





# Ejecutar el script Python para la traduccion
python3 "$SCRIPT_EJ_1" "$ARCHIVO_GENBANK_2" "$ARCHIVO_FASTA_1"


# Verificar si el script Python de la traduccion se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_1' fallo." | tee -a "$ERRORLOGFILE"
	exit 1
fi


echo "La traduccion se completo y se guardo en '$ARCHIVO_FASTA_1'." | tee -a "$LOGFILE"





# Ejecutar el script Python del alineamiento BLAST
python3 "$SCRIPT_EJ_2" "$ARCHIVO_FASTA_1" "$DATABASE_PATH_BLAST" "$ARCHIVO_FASTA_2"


# Verificar si el script Python del BLAST se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_2' fallo." | tee -a "$ERRORLOGFILE"
	exit 1
fi


echo "El archivo FASTA '$ARCHIVO_FASTA_2' se creó correctamente, ordenado por e-value creciente." | tee -a "$LOGFILE"





# Ejecutar el script Python del alineamiento múltiple
python3 "$SCRIPT_EJ_3" "$ARCHIVO_FASTA_1" "$ARCHIVO_FASTA_2" "$ARCHIVO_FASTA_3"


# Verificar si el script Python del BLAST se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_3' fallo." | tee -a "$ERRORLOGFILE"
	exit 1
fi


echo "El alineamiento múltiple se completo y se guardo en '$ARCHIVO_FASTA_3'." | tee -a "$LOGFILE"





# Ejecutar el script Python que utiliza EMBOSS
python3 "$SCRIPT_EJ_4" "$ARCHIVO_GENBANK_2" "$ARCHIVO_FASTA_4" "$ARCHIVO_TXT_SALIDA"


# Verificar si el script Python de EMBOSS se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_4' fallo." | tee -a "$ERRORLOGFILE"
	exit 1
fi


echo "Se obtuvieron los ORFs y Motifs, y se guardaron en '$ARCHIVO_FASTA_4' y '$ARCHIVO_TXT_SALIDA', respectivamente." | tee -a "$LOGFILE"





# Ejecutar el script Python del generador de primers
python3 "$SCRIPT_EJ_5" "$ARCHIVO_GENBANK_2" "$PARAMETROS_PRIMERS" "$ARCHIVO_FASTA_5"


# Verificar si el script Python del BLAST se ejecuto correctamente
if [[ $? -ne 0 ]]; then
	echo "Error: El script '$SCRIPT_EJ_5' fallo." | tee -a "$ERRORLOGFILE"
	exit 1
fi


echo "El alineamiento múltiple se completo y se guardo en '$ARCHIVO_FASTA_5'." | tee -a "$LOGFILE"





# Mover los archivos especificos a la carpeta TPResultados
mv "$ARCHIVO_GENBANK_2" "$RESULTSFOLDER"/ 2>/dev/null
mv "$ARCHIVO_FASTA_1" "$RESULTSFOLDER"/ 2>/dev/null
mv "$ARCHIVO_FASTA_2" "$RESULTSFOLDER"/ 2>/dev/null
mv "$ARCHIVO_FASTA_3" "$RESULTSFOLDER"/ 2>/dev/null
mv "$ARCHIVO_FASTA_4" "$RESULTSFOLDER"/ 2>/dev/null
mv "$ARCHIVO_TXT_SALIDA" "$RESULTSFOLDER"/ 2>/dev/null
mv "$ARCHIVO_FASTA_5" "$RESULTSFOLDER"/ 2>/dev/null

echo "Archivos FASTA guardados en la carpeta TPResultados." | tee -a "$LOGFILE"

mv "$LOGFILE" "$LOGFOLDER"/ 2>/dev/null
mv "$ERRORLOGFILE" "$LOGFOLDER"/ 2>/dev/null
mv "$LOGMUSCLE" "$LOGFOLDER"/ 2>/dev/null


# Actualizo la ruta del archivo LOGFILE después de moverlo
LOGFILE="$LOGFOLDER/$(basename "$LOGFILE")"


echo "Archivos log guardados en la carpeta TPLogs." | tee -a "$LOGFILE"




