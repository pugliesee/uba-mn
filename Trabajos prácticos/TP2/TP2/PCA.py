import matplotlib.pyplot as plt
from Auxiliares import *

# TODO: Ver como linkear la carpeta "caras" para que pathlib las extraiga bien.


def covarianceMatrixPCA():
	X = imgs_apiladas()
	X_flat = aplanar_imgs(X)

	n = len(X_flat[0])
	X_c = centrar_matriz(X_flat)
	# Nos falta crear la matriz de covarianza C
	C = (X_c.T @ X_c) / (n - 1)

	m = len(C)
	f = open("PCATest.txt", "w")
	f.write(str(m))
	f.write(" ")
	f.write(str(m))
	f.write("\n")
	for i in range(m):
		for j in range(m):
			f.write(str(C[i][j]))		# Queremos los autovectores de la matriz de covarianza
			f.write(" ")
		f.write("\n")
	f.close()

	return 		# Generamos un archivo de texto que contenga las imágenes. Luego, ese archivo lo llamaremos desde main.cpp


def covarianceMatrixPCA_reduced():
	X = imgs_apiladas_reduced()
	X_flat = aplanar_imgs(X)

	n = len(X_flat[0])
	X_c = centrar_matriz(X_flat)
	# Nos falta crear la matriz de covarianza C
	C = (X_c.T @ X_c) / (n - 1)

	m = len(C)
	f = open("PCATest_reduced.txt", "w")
	f.write(str(m))
	f.write (" ")
	f.write(str(m))
	f.write("\n")
	for i in range(m):
		for j in range(m):
			f.write(str(C[i][j]))		# Queremos los autovectores de la matriz de covarianza
			f.write (" ")
		f.write("\n")
	f.close()

	return 		# Generamos un archivo de texto que contenga las imágenes. Luego, ese archivo lo llamaremos desde main.cpp


def covarianceMatrixPCA_reduced_sinPrimero():
	X = imgs_apiladas_reduced_sinPrimero()
	X_flat = aplanar_imgs(X)

	n = len(X_flat[0])
	X_c = centrar_matriz(X_flat)
	# Nos falta crear la matriz de covarianza C
	C = (X_c.T @ X_c) / (n - 1)

	m = len(C)
	f = open("PCATest_reduced_sinPrimero.txt", "w")
	f.write(str(m))
	f.write (" ")
	f.write(str(m))
	f.write("\n")
	for i in range(m):
		for j in range(m):
			f.write(str(C[i][j]))		# Queremos los autovectores de la matriz de covarianza
			f.write (" ")
		f.write("\n")
	f.close()

	return 		# Generamos un archivo de texto que contenga las imágenes. Luego, ese archivo lo llamaremos desde main.cpp


def comprimirimagenPCA(X, V_k, i, k):
	# X es la matriz con todas las imagenes apiladas (no aplanadas)
	# V_k es la matriz con los primeros K autovectores de la matriz de covarianza
	# i --> Para seleccionar la i-ésima imagen
	# k --> Para saber sobre cuántos autovectores proyectamos la imágen
	X_flat = aplanar_imgs(X)

	img = []
	j = 0
	while j < k:
		z_i = X_flat[i] @ V_k[:, j]
		img.append(z_i) # Vamos armando el vector que surge de hacer X_i * V_k (es decir, reducimos la dimensionalidad de la imagen
		j += 1			# o simplemente le hacemos un cambio de base en el caso de que k = n)

	ret_var = (X[i], img)

	return ret_var  # Retornamos una tupla cuya primera componente es la imagen original y cuya segunda componente es la versión comprimida de la imagen (vector)


def descomprimir_PCA(z):
	# z es la imagen comprimida
	k = len(z)
	# queremos tomar los primeros k autovectores de la matriz de covarianza
	V_k = mat_eigenvectors_k("autovectores_PCA.txt", k)
	V_k = V_k.T
	res = z @ V_k
	res = np.reshape(res, (112, 92))

	plt.imshow(res, cmap='gray')  # COMENTAR ESTAS LÍNEAS EN CASO DE SER NECESARIO PARA NO GENERARTANTOS GRÁFICOS DURANTE EJECUIÓN
	plt.show()

	return res


def descomprimir_PCA_reduced(z, V_k):
	# z es la imagen comprimida
	res = z @ V_k.T
	res = np.reshape(res, (56, 46))

	plt.imshow(res, cmap='gray')  # COMENTAR ESTAS LÍNEAS EN CASO DE SER NECESARIO PARA NO GENERARTANTOS GRÁFICOS DURANTE EJECUIÓN
	plt.show()

	return res


def descomprimir_PCA_reduced_sinPrimero(z, V_k):
	# z es la imagen comprimida
	res = z @ V_k.T
	res = np.reshape(res, (56, 46))

	plt.imshow(res, cmap='gray')  # COMENTAR ESTAS LÍNEAS EN CASO DE SER NECESARIO PARA NO GENERARTANTOS GRÁFICOS DURANTE EJECUIÓN
	plt.show()

	return res


def imprimirEigenFacePCA(i):
	# X contiene a todas las imágenes sin aplanar
	X = imgs_apiladas()
	a = len(X[0]) # Cantidad de filas de cada imagen
	b = len(X[0][0]) # Cantidad de columnas de cada imagen
	# Se cumple que a*b = n
	V = mat_eigenvectors("autovectores_PCA.txt")
	eface = V[:, i] # Tomamos el i-ésimo autovector de la matriz de covarianza (i-ésima eigenface)
	eface = np.reshape(eface, (a, b))

	plt.imshow(eface, cmap='gray')
	plt.show()

	return


def print_eigenvalues_PCA():
	e_values = []
	f = open("autovalores_PCA.txt", "r")
	for line in f.readlines():
		for nro in line.split():
			e_values.append(float(nro))
	plt.plot(e_values, color="green")
	plt.yscale("log")
	plt.title("Gráfico primeros 410 autovalores PCA")
	plt.show()

	return
