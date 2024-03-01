from PCA import *
import math
from TwoDPCA import *
from Auxiliares import *


def similarityMatrix(X_flat):
    X_c = centrar_matriz(X_flat)
    C_matrix = (X_c @ X_c.T) / (len(X_flat[0]) - 1)

    m = len(X_flat)
    R = np.zeros((m, m))

    i = 0
    while i < m:
        j = 0
        while j < m:
            R[i][j] = C_matrix[i][j] / math.sqrt(C_matrix[i][i] * C_matrix[j][j])
            j += 1
        i += 1

    plt.pcolor(R, cmap='GnBu')
    plt.title("Matriz de similaridad para el conjunto original de imágenes (centrado)")
    plt.show()

    return


def similarityMatrixPCA(X, V_k, k):
    # X --> imágenes apiladas
    # V_k --> matriz con primeros k autovectores de la matriz de covarianza para PCA
    m = len(X)
    a = len(X[0])
    b = len(X[0][0])
    i = 0
    X_c_flat = centrar_matriz(aplanar_imgs(X))
    X_c = []
    for img_flat in X_c_flat:
        X_c.append(np.reshape(img_flat, (a, b)))
    X_c = np.stack(X_c)
    matrizDeComprimidas = []
    while i < m:
        z_i = comprimirimagenPCA(X_c, V_k, i, k)
        matrizDeComprimidas.append(z_i[1])
        i += 1
    matrizDeComprimidas = np.stack(matrizDeComprimidas)

    C = (matrizDeComprimidas @ matrizDeComprimidas.T) / (k-1)

    R = np.zeros((m, m))

    i = 0
    while i < m:
        j = 0
        while j < m:
            R[i][j] = C[i][j] / math.sqrt(C[i][i] * C[j][j])
            j += 1
        i += 1

    plt.pcolor(R, cmap="GnBu")  # COMENTAR ESTAS LÍNEAS EN CASO DE SER NECESARIO PARA NO GENERARTANTOS GRÁFICOS DURANTE EJECUIÓN
    plt.title(f"Matriz de similaridad para PCA con k = {k}")
    plt.colorbar()
    plt.show()

    return R  # la devuelvo porque nos sirve para el ej 3 b)


def similarityMatrix2DPCA(A, U_k, k):
    feature_vecs_mats_flat = []
    m = len(A)
    a = len(A[0])  # en a queda guardada la cantidad de filas de cada imagen
    b = len(A[0][0])  # en b queda guardada la cantidad de columnas de cada imagen
    AVG = np.zeros((a, b))
    i = 0
    while i < m:
        AVG = (AVG + A[i])
        i += 1
    AVG = AVG / m
    A_c = [] # Aquí guardaremos las imágenes centradas
    for img in A:
        A_c.append(img - AVG)
    A_c = np.stack(A_c)
    i = 0
    while i < m:
        i_esima_feat_vecs_mat = comprimir_img_2DPCA(A_c, U_k, i)
        feature_vecs_mats_flat.append(i_esima_feat_vecs_mat.flatten())
        i += 1
    feature_vecs_mats_flat = np.stack(feature_vecs_mats_flat)

    n = len(feature_vecs_mats_flat[0])

    C_mat = (feature_vecs_mats_flat @ feature_vecs_mats_flat.T) / (n - 1)

    R = np.zeros((m, m))

    i = 0
    while i < m:
        j = 0
        while j < m:
            R[i][j] = C_mat[i][j] / math.sqrt(C_mat[i][i] * C_mat[j][j])
            j += 1
        i += 1

    plt.pcolor(R, cmap="GnBu")  # COMENTAR ESTAS LÍNEAS EN CASO DE SER NECESARIO PARA NO GENERARTANTOS GRÁFICOS DURANTE EJECUIÓN
    plt.title(f"Matriz de similaridad para 2DPCA con k = {k}")
    plt.colorbar()
    plt.show()

    return R  # la devuelvo porque nos sirve para el ej 3 b)


def metrica_mismo(R):
    # asumimos que R es la matriz de similaridad de 410 x 410
    # son 10 imágenes por persona
    mat_sum = np.zeros((10, 10))
    n = 10 * 10 * 41  # son 10 x 10 coordenadas para 41 personas
    acum = 0
    i = 0
    j = 0
    while i < 410:  # vamos a ir tomando las matrices de 10 x 10 de la diagonal
        mat_sum += R[i:i+10, j:j+10]
        i += 10
        j += 10

    # ya tenemos la matriz de 10 x 10 con los valores acumulados
    # ahora hay que sumar todas las coordenadas de esa matriz para devolver un escalar
    for fila in mat_sum:
        for nro in fila:
            acum += nro

    acum = acum / n

    return acum


def metrica_distinto(R):
    # asumimos que R es la matriz de similaridad de 410 x 410
    # son 10 imágenes por persona
    # IDEA: la matriz R es simétrica. Aprovechando esto, vamos a sumar todas las coordenadas del triángulo superior
    # (sin la diagonal) y luego vamos a multiplicar por 2.
    mat_sum = np.zeros((10, 10))
    n = 410 * 410 - (10 * 10 * 41)  # vamos a sumar todas las coordenadas de la matriz
                                    # salvo las de los bloques de 10 x 10 en la diagonal
    acum = 0
    i = 0
    while i < 410:
        j = i + 10
        while j < 410:
            mat_sum += R[i:i+10, j:j+10]
            j += 10
        i += 10

    # ya sumamos las coordenadas del triángulo superior
    # ahora sumamos a la matriz con sí misma para multiplicar los valores por 2
    # y así tener la suma de ambos triángulos de la matriz
    mat_sum += mat_sum

    # ya tenemos la matriz de 10 x 10 con los valores acumulados
    # ahora hay que sumar todas las coordenadas de esa matriz para devolver un escalar
    for fila in mat_sum:
        for nro in fila:
            acum += nro

    acum = acum / n

    return acum
