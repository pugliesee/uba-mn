from TwoDPCA import *
from Matriz_Similaridad import *
from Auxiliares import *

# ITEM D
X = imgs_apiladas()
# PCA
k_PCA = 420
nro_imagen_PCA = 181
V_k = mat_eigenvectors_k("autovectores_PCA.txt", k_PCA)
z_PCA = comprimirimagenPCA(X, V_k, nro_imagen_PCA, k_PCA)  # z_PCA[1] es el vector con la imagen comprimida
z_descomp_PCA = descomprimir_PCA(z_PCA[1])

# 2DPCA
k_2DPCA = 90
nro_imagen_2DPCA = 260
U_k = mat_eigenvectors_k("autovectores_2DPCA.txt", k_2DPCA)
z_2DPCA = comprimir_img_2DPCA(X, U_k, nro_imagen_2DPCA)  # z_2DPCA tiene como columnas a los primeros k feature vectors de la imagen y es de a x k
z_descomp_2DPCA = descomprimir_2DPCA(z_2DPCA, U_k)
