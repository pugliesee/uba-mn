from Matriz_Similaridad import *
from Auxiliares import *

# ITEM A
X = imgs_apiladas()

# PCA
# covarianceMatrixPCA()  # luego de esta función habría que correr el método de la potencia con delfación en c++
k_PCA = 1  # ir cambiándolos para probar con distintos parámetros e imágenes.
nro_imagen_PCA = 260
V_k = mat_eigenvectors_k("autovectores_PCA.txt", k_PCA)
z_PCA = comprimirimagenPCA(X, V_k, nro_imagen_PCA, k_PCA)  # z_PCA[1] es el vector con la imagen comprimida
print(z_PCA[1])

# 2DPCA
# imageCovarianceMatrix()  # luego de esta función habría que correr el método de la potencia con delfación en c++
k_2DPCA = 1 # ir cambiándolo para probar con distintos parámetros e imágenes
nro_imagen_2DPCA = 260
U_k = mat_eigenvectors_k("autovectores_2DPCA.txt", k_2DPCA)
z_2DPCA = comprimir_img_2DPCA(X, U_k, nro_imagen_2DPCA)  # z_2DPCA tiene como columnas a los primeros k feature vectors de la imagen y es de a x k
print(z_2DPCA)
