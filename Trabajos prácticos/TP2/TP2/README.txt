Comenzando por Python, abrimos PyCharm, Visual Studio Code, o el intérprete que usemos
y cargamos toda la carpeta para que todos los archivos
se carguen juntos.

Una vez abierto el proyecto, nos vamos a "main.py"
Ejecutando únicamente la función "covarianceMatrixPCA()" vamos a generar
la matriz de covarianza de PCA.
Lo mismo sucede para 2DPCA si ejecutamos también "imageCovarianceMatrix()"

Luego tenemos para el ejercicio 3 c) unas 3 funciones más: "covarianceMatrixPCA_reduced()", "covarianceMatrixPCA_reduced_sinPrimero()" y "imageCovarianceMatrix_sinPrimero()"


Todas estas funciones nos generarán 5 archivos de texto:

--> "PCATest.txt"
--> "2DPCA.txt"
--> "PCATest_reduced.txt"
--> "PCATest_reduced_sinPrimero.txt"
--> "2DPCA_sinPrimero.txt"

Con ellos, podemos ejecutar "main.cpp", para poner en funcionamiento el método
de potencia-deflación y conseguir los autovalores y autovectores de cada matriz.
Adentro de "main.cpp" se explica el funcionamiento del código, y como se puede editar, según se desee.

Luego de correr "main.cpp", obtendremos 10 archivos en total

- Para PCA:
 --> "autovalores_PCA.txt"
 --> "autovectores_PCA.txt"
 --> "autovalores_PCA_reduced.txt"
 --> "autovectores_PCA_reduced.txt"
 --> "autovalores_PCA_reduced_sinPrimero.txt"
 --> "autovectores_PCA_reduced_sinPrimero.txt"

- Para 2DPCA:
 --> "autovalores_2DPCA.txt"
 --> "autovectores_2DPCA.txt"
 --> "autovalores_2DPCA_sinPrimero.txt"
 --> "autovectores_2DPCA_sinPrimero.txt"

Ahora podemos desplazarnos libremente a lo largo de los archivos:
"EJ_2A.py", "EJ_2B.py", "EJ_2C.py", "EJ_2D.py"
y tan solo ejecutarlos, ya que nos darán los resultados que esperamos.
(De igual forma, en el informe se incluyeron los debidos gráficos, fotos y explicaciones)
Todas las funciones llamadas en estos archivos pueden verse implementadas en
- "PCA.py"
- "TwoDPCA.py"

Podemos hacer lo mismo si nos dirigimos a los archivos:
"EJ_3A.py", "EJ_3B.py", "EJ_3C.py", "EJ_3D.py"
y observar los gráficos esperados.
En el archivo "Matriz_Similaridad.py" se encuentran las implementaciones usadas a lo largo
de estos archivos.

ACLARACIONES:
- Para correr el código de C++ en CLion, seleccionar la opción llamada TP2 al lado del botón de play
- Tuvimos unos problemas al principio a la hora de incluir Eigen en nuestros archivos de C++, y a cada uno nos funcionaba de manera distinta.
Si no funciona el include tal como está en el código, probar con #include <Eigen/Dense> o directamente con #include <Dense>.
- A veces, al tener el PyCharm y el CLion abiertos al mismo tiempo, se genera un pop-up en el PyCharm indicando que no se pudo
cargar correctamente el módulo TP2. Clickeando en el link del pop-up, nos dice que se debe a un archivo .iml y preguntas si deseamos borrarlo.
Cuando aparezca dicha opción, borrar ese archivo.
- En "EJ_3B.py" se usan muchos Ks para generar distintas matrices de similaridad.
- Si se desea, podemos "acortar" la lista, o usar un único K en el caso de no querer crear
tantos gráficos.
- Otra opción es dirigirse al archivo "Matriz_Similaridad.py" y comentar en las funciones correspondientes la líneas que imprimen los gráficos




