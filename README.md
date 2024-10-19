# Bioinfo
EL ejercicio 4 esta listo para usar en linux. Si queiren la verssión colab aca el link: https://colab.research.google.com/drive/1zYaHwGXvc6yfroYE7PSQGOJOBWSu5b5f?usp=sharing.
Si tenes todo ok solo tenes que hacer una cosa. Correr el programa en conjunto con el archivo del genbank, en mi caso opsina.gb.
Si no tenes todo okk tenes que instalar varios paquetes de python, entre ellos: emboss, biophyton y lo que te pida. para instalar biopython use: sudo apt install python3-biopython
Ademas, para poner a punto el patmatmotifs primero tenes que instalar el prosextract. Para este lo que tenes que hacer es hacerte una carpeta llamada PROSITE y adentro meterle los archivos de la base de datos de prosite "prosite.dat" y "prosite.doc". 
Una vez que tenes todo ok el scrip te genera dos archivos, el archivo con los orfs y el archivo con los motivos de los orf que tuvieron al meos un hit. En la opsina son 2 no más.
para correrlo en linux asegurate de tener todo en una carpeta sola. Entonces en esta carpeta tendrías: EJ4Marcos.py opsina.gb PROSITE/    siendo PROSITE/ la carpeta con los archivos prosite.dat y prosite.doc
se corre con esta linea: python3 EJ4Marcos.py opsina.gb
errores posibles: se queda en las lineas de libreria? hay que instalar un paquete
                  patmatmotifs no labura? chequeate lo de prosite y prosextract
basicamente eso es todo!!


PD: Si ya hiciste una vez lo de prosite no hace falta que lo hagas más. Es más evitalo que tira un error boludo. Basicamente elimina las lineas donde aparezca PROSITE y prosextract....SOLO SI YA LO INSTALASTE
