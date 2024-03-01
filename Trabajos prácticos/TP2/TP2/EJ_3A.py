from Matriz_Similaridad import *
from Auxiliares import *

# ITEM A
X = imgs_apiladas()
# Matriz de correlación para las imágenes en formato original (centradas)
X_flat = aplanar_imgs(X)
similarityMatrix(X_flat)

# Matriz de correlación para PCA
# k_PCA = 400
# V_k_PCA = mat_eigenvectors_k("autovectores_PCA.txt", k_PCA)
# similarityMatrixPCA(X, V_k_PCA, k_PCA)

# Matriz de correlación para 2DPCA
# k_2DPCA = 30
# U_k_2DPCA = mat_eigenvectors_k("autovectores_2DPCA.txt", k_2DPCA)
# similarityMatrix2DPCA(X, U_k_2DPCA, k_2DPCA)

