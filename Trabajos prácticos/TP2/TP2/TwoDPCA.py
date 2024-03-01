import matplotlib.pyplot as plt

from Auxiliares import *


def imageCovarianceMatrix():
	A = imgs_apiladas()  # A queda como un arreglo tridimensional de m x (a x b)
	i = 0
	n = len(A)  # en n queda guardada la cantidad de imágenes
	a = len(A[0])  # en a queda guardada la cantidad de filas de cada imagen
	b = len(A[0][0])  # en b queda guardada la cantidad de columnas de cada imagen
	AVG = np.zeros((a, b))
	while i < n:
		AVG = (AVG + A[i])
		i += 1
	AVG = AVG / n

	# Ya tenemos "A-barra"
	j = 0
	G = np.zeros(()) # la matriz G queda de b x b
	while j < n:		# Vamos a calcular la matriz G
		temp = (A[j] - AVG).T @ (A[j] - AVG)
		G = G + temp
		j += 1
	G = G / n		# Ahora queremos los autovectores de G

	m = len(G) # con una sola variable alcanza pues G es cuadrada
	f = open("2DPCA.txt", "w")		# Creamos el .txt y luego llamamos a main.cpp para calcularlos
	f.write(str(m))
	f.write(" ")
	f.write(str(m))
	f.write("\n")
	for i in range(m):
		for j in range(m):
			f.write(str(G[i][j]))
			f.write(" ")
		f.write("\n")
	f.close()

	return A


def imageCovarianceMatrix_sinPrimero():
	A = imgs_apiladas_sinPrimero()  # A queda como un arreglo tridimensional de m x (a x b)
	i = 0
	n = len(A)  # en n queda guardada la cantidad de imágenes
	a = len(A[0])  # en a queda guardada la cantidad de filas de cada imagen
	b = len(A[0][0])  # en b queda guardada la cantidad de columnas de cada imagen
	AVG = np.zeros((a, b))
	while i < n:
		AVG = (AVG + A[i])
		i += 1
	AVG = AVG / n

	# Ya tenemos "A-barra"
	j = 0
	G = np.zeros(()) # la matriz G queda de b x b
	while j < n:		# Vamos a calcular la matriz G
		temp = (A[j] - AVG).T @ (A[j] - AVG)
		G = G + temp
		j += 1
	G = G / n		# Ahora queremos los autovectores de G

	m = len(G) # con una sola variable alcanza pues G es cuadrada
	f = open("2DPCA_sinPrimero.txt", "w")		# Creamos el .txt y luego llamamos a main.cpp para calcularlos
	f.write(str(m))
	f.write(" ")
	f.write(str(m))
	f.write("\n")
	for i in range(m):
		for j in range(m):
			f.write(str(G[i][j]))
			f.write(" ")
		f.write("\n")
	f.close()

	return A


def comprimir_img_2DPCA(X, U_k, nro_imagen):
	A = X[nro_imagen]  # imagen de a x b
	V = A @ U_k   # V tiene en sus columnas a los primeros k feature vectors de la imagen

	return V


def print_subimage_twoDPCA(A, i, k):
	# i --> Para la i-ésima imágen. k --> Para el k-ésimo feature
	# Acá vamos a calcular V = [Y1,...,Yb] y U = [X1,...,Xb]
	# Y_i siendo el i-ésimo feature vector y X_i el i-ésimo autovector ortogonal de G
	# REMINDER: Y = AX
	a = len(A[0])
	b = len(A[0][0])
	U = mat_eigenvectors("autovectores_2DPCA.txt")  # en U me quedan todos los autovectores ortonormales de G

	V = A[i] @ U 			# V es la matriz cuyas columnas son los feature vectors
	R = np.zeros((a, b))
	k_feature = V[:, k]
	k_feature = np.reshape(k_feature, (len(k_feature), 1))
	k_autovec = U[:, k]
	k_autovec = np.reshape(k_autovec, (1, len(k_autovec)))
	e_face = k_feature @ k_autovec
	R = R + e_face

	plt.imshow(R, cmap='gray')
	plt.show()

	return


def descomprimir_2DPCA(z_2DPCA, U_k):
	# En esta función vamos a sumar k subimágenes para reconstruir la imagen comprimida
	# k es el nro de feature vectors
	# si multiplicamos el i-ésimo feature vector por el i-ésimo autovector de la matriz de covarianza traspuesto,
	# obtendremos la i-ésima subimagen
	a = len(z_2DPCA)
	b = len(U_k)
	k = len(z_2DPCA[0])
	R = np.zeros((a, b))
	j = 0
	while j < k:
		j_feature = z_2DPCA[:, j]
		j_feature = np.reshape(j_feature, (len(j_feature), 1))
		j_autovec = U_k[:, j]
		j_autovec = np.reshape(j_autovec, (1, len(j_autovec)))
		e_face = j_feature @ j_autovec
		R = R + e_face
		j += 1

	# plt.imshow(R, cmap='gray')
	# plt.show()

	return R


def print_eigenvalues_twoDPCA():
	e_values = []
	f = open("autovalores_2DPCA.txt", "r")
	for line in f.readlines():
		for nro in line.split():
			e_values.append(float(nro))
	plt.plot(e_values)
	plt.title("Gráfico autovalores 2DPCA")
	plt.yscale("log")
	plt.show()

	return
