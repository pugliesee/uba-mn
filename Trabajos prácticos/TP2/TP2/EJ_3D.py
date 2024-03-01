from Auxiliares import *

# ITEM D

# GRÁFICO DE TIEMPO PARA PCA -- UTILIZANDO EL MÉTODO DE POTENCIA-DEFLACIÓN

x = [1278.62, 5922.43, 10423.3, 19975.9, 31585.6]	 # Tiempos calculados desde C++
y = [10, 100, 200, 500, 1000]						 # Cantidad de iteraciones utilizadas
plt.plot(x, y, marker="o", label="PCA", color="orchid")
plt.title("Medición de tiempos del método potencia-deflación con PCA")
plt.xlabel("Segundos")
plt.ylabel("Iteraciones")
plt.legend()
plt.show()


# GRÁFICO DE TIEMPO PARA 2DPCA -- UTILIZANDO EL MÉTODO DE POTENCIA-DEFLACIÓN

x = [0.130564, 0.1617879, 0.16941, 0.180686, 0.2601]	 # Tiempos calculados desde C++
y = [10, 100, 200, 500, 1000]  							 # Cantidad de iteraciones utilizadas
plt.plot(x, y, marker="o", label="2DPCA", color="royalblue")
plt.title("Medición de tiempos del método potencia-deflación con 2DPCA")
plt.xlabel("Segundos")
plt.ylabel("Iteraciones")
plt.legend()
plt.show()

