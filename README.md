# Bioinfo

El trabajo comprende los cinco scripts de python (uno por cada ejercicio), el script mutador, el param.xml con los parámetros para el diseño
de primers, el arhivo genbank de la Opsina y el script orquestador de bash. Es necesario que los archivos mencionados se encuentren todos en una misma carpeta, en una pc con linux.


Antes de ejecutar el orquestador, es necesario realizar los siguientes pasos para que todos los scripts funcionen correctamente:

Al momento de ejecutar los scripts, tener en cuenta que se utilizó la versión 3.6.9 de Python para el armado de los mismos.

## Ex1 

   Instalar Biopython  (Desde terminal: python3 -m pip install biopython)

## Ex2 

   Se requiere tener Blast+ y una copia de la base de datos de swissprot instaladas. En las siguientes líneas se dará un paso a paso:
   
   - Crear 2 directorios, uno destinado al aplicativo de Blast+ y otro destinado a la base de datos 
     (Ejemplo: mkdir Apps/Blast+ y mkdir Apps/Blastdb).

   - Ingresar al directorio de Blast+ (cd Apps/Blast+) y ejecutar la siguiente línea para descargar el aplicativo dentro del directorio:
     (wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz).

   - Descomprimir el aplicativo de Blast+ (tar -xzvf ncbi-blast-2.16.0+-x64-linux.tar.gz).

   - El siguiente paso implica agregar el path del aplicativo, a la configuración inicial del entorno de ejecución. Ésto permitirá ejecutar
     cualquiera de las aplicaciones de Blast+ desde cualquier lugar de la terminal, sin necesidad de escribir la ruta completa de cada
     aplicación. Lo primero es obtener el path de la carpeta bin del aplicativo (ej: /home/fran/Apps/Blast+/ncbi-blast-2.16.0+/bin), lo
     siguiente es abrir el script bashrc con un editor (vi ~/.bashrc) Al final de dicho script debe incorporarse la siguiente línea:
     (export PATH="$PATH:/home/fran/Apps/Blast+/ncbi-blast-2.16.0+/bin"), siendo PATH: el correspondiente a la ubicación de la carpeta bin
     del aplicativo de Blast+. Una vez incorporada dicha línea de código, guardar los cambios y cerrar el script.

   - Ingresar al directorio de Blastdb (cd Apps/Blastdb) y ejecutar la siguiente línea para descargar la base de datos dentro del directorio:
     (wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz).

   - Descomprimir el archivo que contiene la base de datos (gunzip swissprot.gz).

   - La base de datos es un .FASTA, debe cambiarse su formato utilizando la aplicación makeblastdb de Blast+. Desde el directorio de la 
     base de datos, aplicar el siguiente comando en la terminal: (makeblastdb -in swissprot -dbtype prot -out swissprot_db -parse_seqids).

   - Por último, incorporar el path de la base de datos al archivo bash. Abrir el script (vi ~/.bashrc) e incorporar al final la siguiente
     línea: (export BLASTDB_PATH="/home/fran/Apps/Blastdb/swissprot_db") guardar cambios y cerrar.
     

## Ex3

   Instalar MUSCLE para poder realizar el alineamiento múltiple (sudo apt-get install muscle).

## Ex4
   - Primero se requiere instalar EMBOSS, para ello desde la consola descargar el zip que contiene el programa 
     (wget ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz).

   - Descomprimir el archivo descargado (tar -zxvf emboss-latest.tar.gz).

   - El siguiente paso es compilar EMBOSS. Ingresar al directorio que se ha descomprimido. Una vez dentro ejecutar el siguiente comando en la
     terminal: (./configure). Éste último prepara el entorno de compilación, detectando las características del sistema y asegurándose de que
     estén presentes todas las dependencias necesarias. En caso de faltar alguna dependencia, se notificará en la terminal y será necesario
     instalarla previo a realizar el siguiente paso (En mi caso se notificó que no se encontraba instalado en mi sistema las bibliotecas de
     desarrollo de X11, por lo que tuve que instalar las mismas: sudo apt-get install libx11-dev).

   - Una vez instaladas las dependencias que requiera el archivo configure, ejecutar nuevamente (./configure) dentro del directorio 
     descomprimido de EMBOSS. Seguidamente compilar ejecutando (make), y por último (sudo make install).

   - Es necesario incorporar al script bashrc, las rutas de la carpeta bin donde se encuentran los ejecutables de EMBOSS y la carpeta
     lib con las bibliotecas compartidas de EMBOSS. Abrir el script de bash (vi ~/.bashrc) e incorporar al final las siguientes dos
     líneas (export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib) siendo /usr/local/lib la ubicación de las biblioteca de EMBOSS, y
     (export PATH="$PATH:/home/fran/Apps/Blast+/ncbi-blast-2.16.0+/bin:usr/local/bin"), siendo :usr/local/bin la ubicación de los ejecutables
     de EMBOSS. En éste último caso ya tenía definido por el ejercicio 2 un export PATH. Es por ello que a éste último le incorporo al final
     el :usr/local/bin, en lugar de crear un (export PATH="$PATH:usr/local/bin")

   - Si la configuración se dió de manera correcta, al escribir (embossversion) en la terminal se debería imprimir la versión de EMBOSS 
     instalada.

   - El siguiente paso es crear un directorio en la cual se guardarán las bases de datos de PROSITE. Para ello crear el directorio
     (mkdir Apps/PROSITE), ingresar al directorio creado y descargar los siguientes dos archivos:
     (wget https://ftp.expasy.org/databases/prosite/prosite.doc) y (wget https://ftp.expasy.org/databases/prosite/prosite.dat)

   - Por último se requiere configurar a EMBOSS para que, al momento de determinar motifs, utilice la base de datos descargada de PROSITE.
     Para ello es necesario ejecutar el siguiente comando desde la terminal (sudo prosextract -prositedir /home/fran/Apps/PROSITE), donde
     /home/fran/Apps/PROSITE corresponde a la ubicación del directorio creado en el inciso anterior.
     


## Ex5 
   
   Instalar primer3-py, la versión 0.6.1 (Desde terminal: python3 -m pip install primer3-py==0.6.1).



Al ejecutar el orquestador, se crearán dos directorios, uno en el cual se alojarán todos los archivos .log (TP Logs) y otro en el
cual se podrán encontrar todos los archivos generados por los scripts (Tp Resultados).



