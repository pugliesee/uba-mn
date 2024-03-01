from Matriz_Similaridad import *
from ECM import *
from Auxiliares import *

# ITEM C
# VERSIÓN CON TODAS LAS IMÁGENES
X = imgs_apiladas()  # para 2DPCA
X_reduced = imgs_apiladas_reduced()  # para PCA
# covarianceMatrixPCA_reduced()  # lo ejecutamos una única vez para poder generar la matriz de autovectores con C++

k_list_2DPCA = [1, 5, 10, 15, 20, 30, 45, 60, 75, 92]
k_list_PCA = [10, 50, 100, 200, 410, 600, 800, 1000, 1500, 2000, 2576]

# PCA
m = len(X_reduced)
valores_ecm_PCA = []
for k in k_list_PCA:
    acum_err = 0
    V_k = mat_eigenvectors_k("autovectores_PCA_reduced.txt", k)
    i = 0
    while i < m:
        z = comprimirimagenPCA(X_reduced, V_k, i, k)
        z_descomp = descomprimir_PCA_reduced(z[1], V_k)
        err = ECM(X_reduced[i], z_descomp) # calculamos el ecm para cada imagen del dataset
        acum_err += err  # y lo vamos acumulando
        i += 1

    acum_err = acum_err / m # finalmente calculamos el promedio del ecm para todas las imágenes
    valores_ecm_PCA.append(acum_err)
    print(f"Listo con k = {k}")

plt.plot(k_list_PCA, valores_ecm_PCA, marker="o")
plt.title("Error de compresión con ECM para PCA con dataset completo")
plt.xlabel("k")
plt.ylabel("err")
plt.show()

print(valores_ecm_PCA)

# 2DPCA
m = len(X)
valores_ECM_2DPCA = []
for k in k_list_2DPCA:
    acum_err = 0
    U_k = mat_eigenvectors_k("autovectores_2DPCA.txt", k)
    i = 0
    while i < m:
        z = comprimir_img_2DPCA(X, U_k, i)
        z_descomp = descomprimir_2DPCA(z, U_k)
        err = ECM(X[i], z_descomp)
        acum_err += err
        i += 1

    acum_err = acum_err / m
    valores_ECM_2DPCA.append(acum_err)
    print(f"Listo con k = {k}")

plt.plot(k_list_2DPCA, valores_ECM_2DPCA, marker="o", color="purple")
plt.title("Error de compresión con ECM para 2DPCA con dataset completo")
plt.xlabel("k")
plt.ylabel("err")
plt.show()

print(valores_ECM_2DPCA)

# VERSIÓN CON EL DATASET SIN LAS IMÁGENES DE LA PERSONA 1
X_soloPrimero = imgs_apiladas_soloPrimero()
X_reduced_soloPrimero = imgs_apiladas_reduced_soloPrimero()
X_sinPrimero= imgs_apiladas_sinPrimero()
X_reduced_sinPrimero = imgs_apiladas_reduced_sinPrimero()
# covarianceMatrixPCA_reduced_sinPrimero()
# imageCovarianceMatrix_sinPrimero()

# PCA
m = len(X_reduced_soloPrimero)
valores_ecm_PCA_sinPrimero = []
for k in k_list_PCA:
    acum_err = 0
    V_k = mat_eigenvectors_k("autovectores_PCA_reduced_sinPrimero.txt", k)
    i = 0
    while i < m:
        z = comprimirimagenPCA(X_reduced_soloPrimero, V_k, i, k)
        z_descomp = descomprimir_PCA_reduced_sinPrimero(z[1], V_k)
        err = ECM(X_reduced_soloPrimero[i], z_descomp) # calculamos el ecm para cada imagen de la primera persona
        acum_err += err  # y lo vamos acumulando
        i += 1

    acum_err = acum_err / m # finalmente calculamos el promedio del ecm para todas las imágenes
    valores_ecm_PCA_sinPrimero.append(acum_err)
    print(f"Listo con k = {k}")

plt.plot(k_list_PCA, valores_ecm_PCA_sinPrimero, marker="o")
plt.title("Error de compresión con ECM para PCA con dataset sin persona 1 \n para imágenes de persona 1")
plt.xlabel("k")
plt.ylabel("err")
plt.show()

print(valores_ecm_PCA_sinPrimero)

# 2DPCA
m = len(X_soloPrimero)
valores_ECM_2DPCA_sinPrimero = []
for k in k_list_2DPCA:
    acum_err = 0
    U_k = mat_eigenvectors_k("autovectores_2DPCA_sinPrimero.txt", k)
    i = 0
    while i < m:
        z = comprimir_img_2DPCA(X_soloPrimero, U_k, i)
        z_descomp = descomprimir_2DPCA(z, U_k)
        err = ECM(X_soloPrimero[i], z_descomp)
        acum_err += err
        i += 1

    acum_err = acum_err / m
    valores_ECM_2DPCA_sinPrimero.append(acum_err)
    print(f"Listo con k = {k}")

plt.plot(k_list_2DPCA, valores_ECM_2DPCA_sinPrimero, marker="o", color="purple")
plt.title("Error de compresión con ECM para 2DPCA con dataset sin persona 1 \n para imágenes de persona 1")
plt.xlabel("k")
plt.ylabel("err")
plt.show()

print(valores_ECM_2DPCA_sinPrimero)
