from Matriz_Similaridad import *
import matplotlib.pyplot as plt
from Auxiliares import *

# ITEM B
X = imgs_apiladas()
# PCA
k_list_PCA = [10, 50, 100, 200, 410, 600, 800, 1000, 2500]
metrica_mismo_PCA = []
metrica_distinto_PCA = []
# Obtenemos las métricas para los distintos valores de k
for k in k_list_PCA:
    V_k_PCA = mat_eigenvectors_k("autovectores_PCA.txt", k)
    R = similarityMatrixPCA(X, V_k_PCA, k)
    metrica_mismo_PCA.append(metrica_mismo(R))
    metrica_distinto_PCA.append(metrica_distinto(R))
# Genero los gráficos
plt.plot(k_list_PCA, metrica_mismo_PCA, marker="o", label="Métrica mismo", color="red")
plt.plot(k_list_PCA, metrica_distinto_PCA, marker="o", label="Métrica distinto", color="green")
plt.title("Métricas para el análisis de similaridad en PCA")
plt.xlabel("k")
plt.ylabel("Métricas")
plt.yscale("log")
plt.legend()
plt.show()

print(metrica_mismo_PCA)
print(metrica_distinto_PCA)

# 2DPCA
k_list_2DPCA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 60, 70, 80, 92]
metrica_mismo_2DPCA = []
metrica_distinto_2DPCA = []
# Obtenemos las métricas para los distintos valores de k
for k in k_list_2DPCA:
    U_k_2DPCA = mat_eigenvectors_k("autovectores_2DPCA.txt", k)
    R = similarityMatrix2DPCA(X, U_k_2DPCA, k)
    metrica_mismo_2DPCA.append(metrica_mismo(R))
    metrica_distinto_2DPCA.append(metrica_distinto(R))

plt.plot(k_list_2DPCA, metrica_mismo_2DPCA, marker="o", label="Métrica mismo", color="blue")
plt.plot(k_list_2DPCA, metrica_distinto_2DPCA, marker="o", label="Métrica distinto", color="orange")
plt.title("Métricas para el análisis de similaridad en 2DPCA")
plt.xlabel("k")
plt.ylabel("Métricas")
plt.yscale("log")
plt.legend()
plt.show()

print(metrica_mismo_2DPCA)
print(metrica_distinto_2DPCA)
