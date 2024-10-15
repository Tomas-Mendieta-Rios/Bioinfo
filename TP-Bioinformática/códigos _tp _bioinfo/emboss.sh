#!/bin/bash

# Verificar que se ha proporcionado un archivo de secuencia
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 archivo_secuencia.gb"
    exit 1
fi

# Archivo de entrada
archivo_secuencia="$1"
archivo_fasta="${archivo_secuencia%.gb}.fasta"

# Definir directorios y archivos
prosite_dir=$(pwd)
log_file="$prosite_dir/prosite.log"
resultados_dir="$prosite_dir/resultados"
orf_file="$resultados_dir/orf_secuencias.fasta"


# Crear los directorios necesarios
mkdir -p "$resultados_dir"

# Convertir GenBank a FASTA
echo "Convirtiendo archivo GenBank a FASTA..."
seqret -sequence "$archivo_secuencia" -outseq "$archivo_fasta" -auto 2>>"$log_file"

# Verificar si el archivo FASTA fue creado correctamente
if [ -f "$archivo_fasta" ]; then
    echo -e "\033[1;32mConversión a FASTA completada. Archivo FASTA guardado en $archivo_fasta.\033[0m"
else
    echo -e "\033[1;31mFallo al convertir el archivo GenBank a FASTA.\033[0m" | tee -a "$log_file"
    exit 1
fi

# Calcular los ORF y obtener las secuencias de proteínas
echo "Calculando ORF y obteniendo secuencias de proteínas..."
getorf -sequence "$archivo_fasta" -outseq "$orf_file" -minsize 100 2>>"$log_file"

# Verificar si el archivo ORF fue creado correctamente
if [ -f "$orf_file" ]; then
    echo -e "\033[1;32mORF calculados y secuencias de proteínas guardadas en $orf_file.\033[0m"
else
    echo -e "\033[1;31mFallo al calcular ORF o guardar secuencias de proteínas en $orf_file.\033[0m" | tee -a "$log_file"
    exit 1
fi

#obtener dominios
chmod 777 dominios.py
python3 dominios.py "$orf_file"

