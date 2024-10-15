# Verifica que se ha pasado el argumento correcto
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 <path_del_archivo_genbank>"
    exit 1
fi

# Obtén el path del archivo GenBank
path_gb="$1"

# Define la ruta del archivo de salida fasta
archivo_salida_fasta="ORFs.fasta"

# Elimina el archivo de salida fasta si ya existe
if [ -f "$archivo_salida_fasta" ]; then
    rm "$archivo_salida_fasta"
fi

#Ejecuta el primer script de Python para generar el archivo fasta
chmod 777 ejercicio_1.py
python3 ejercicio_1.py "$path_gb" "$archivo_salida_fasta"

#Ejecuta el segundo script de Python para  realizar las búsquedas BLAST
chmod 777 ejercicio_2.py
python3 ejercicio_2.py "$archivo_salida_fasta"

#Ejecuta el tercer script de Python para realizar el MSA
chmod 777 ejercicio_3.py
conda install -c bioconda clustalo
python3 ejercicio_3.py

#Ejecuta el quinto script de Python para generar los primers
pip install primer3-py
chmod 777 ejercicio_5.py
#Ejecutar la creación de primers
python3 ejercicio_5.py "$path_gb" 



