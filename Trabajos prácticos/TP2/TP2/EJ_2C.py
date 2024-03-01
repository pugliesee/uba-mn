from Matriz_Similaridad import *
from Auxiliares import *

# ITEM C

# PCA
nro_eigenface_PCA = 0  # para ir probando con distintos valores
imprimirEigenFacePCA(nro_eigenface_PCA)

# 2DPCA
A = imgs_apiladas()
nro_imagen = 0
nro_subimagen_2DPCA = 0
print_subimage_twoDPCA(A, nro_imagen, nro_subimagen_2DPCA)
