import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def imgs_apiladas():
    paths = []  # Vamos almacenando las imágenes una por una
    imgs = []
    for path in sorted(list(Path('./caras').rglob('*/*.pgm'))):
        paths.append(path)
        imgs.append(plt.imread(path)[::2, ::2] / 255)
    X = np.stack(imgs)  # Apilamos las imágenes, una arriba de la otra
    # Conseguimos una matriz X de tamaño M x N.

    return X


def imgs_apiladas_soloPrimero():
    paths = []  # Vamos almacenando las imágenes una por una
    imgs = []
    for path in sorted(list(Path('./caras/s1').rglob('*.pgm'))):
        paths.append(path)
        imgs.append(plt.imread(path)[::2, ::2] / 255)
    X = np.stack(imgs)  # Apilamos las imágenes, una arriba de la otra
    # Conseguimos una matriz X de tamaño 10 x N.

    return X


def imgs_apiladas_sinPrimero():
    paths = []  # Vamos almacenando las imágenes una por una
    imgs = []
    for path in sorted(list(Path('./caras').rglob('*/*.pgm'))):
        paths.append(path)
        imgs.append(plt.imread(path)[::2, ::2] / 255)
    imgs = imgs[10:]  # saco las 10 fotos de la primera persona
    X = np.stack(imgs)  # Apilamos las imágenes, una arriba de la otra
    # Conseguimos una matriz X de tamaño (M-1) x N.

    return X


def imgs_apiladas_reduced():
    paths = []  # Vamos almacenando las imágenes una por una
    imgs = []
    for path in sorted(list(Path('./caras').rglob('*/*.pgm'))):
        paths.append(path)
        imgs.append(plt.imread(path)[::2, ::2] / 255)
    X = np.stack(imgs)  # Apilamos las imágenes, una arriba de la otra
    # Conseguimos una matriz X de tamaño M x N.

    return X


def imgs_apiladas_reduced_soloPrimero():
    paths = []  # Vamos almacenando las imágenes una por una
    imgs = []
    for path in sorted(list(Path('./caras/s1').rglob('*.pgm'))):
        paths.append(path)
        imgs.append(plt.imread(path)[::2, ::2] / 255)
    X = np.stack(imgs)  # Apilamos las imágenes, una arriba de la otra
    # Conseguimos una matriz X de tamaño 10 x N.

    return X


def imgs_apiladas_reduced_sinPrimero():
    paths = []  # Vamos almacenando las imágenes una por una
    imgs = []
    for path in sorted(list(Path('./caras').rglob('*/*.pgm'))):
        paths.append(path)
        imgs.append(plt.imread(path)[::2, ::2] / 255)
    imgs = imgs[10:]  # saco las 10 fotos de la primera persona
    X = np.stack(imgs)  # Apilamos las imágenes, una arriba de la otra
    # Conseguimos una matriz X de tamaño (M-1) x N.

    return X


def aplanar_imgs(X):
    X_flat = []
    for i in range(len(X)):
        X_flat.append(X[i].flatten())  # Ahora sí cada imagen nos queda como un solo vector horizontal
    X_flat = np.stack(X_flat)

    return X_flat


def centrar_matriz(X):
    # debemos calcular el promedio (μ) de cada columna, y restarselo.
    j = 0
    n = len(X[0])
    res = []
    while j < n:
        Mu_j = np.average(X[:, j])  # Tomamos el j-ésimo píxel de c/ imagen y hacemos el promedio
        res.append(Mu_j)
        j += 1
    j = 0
    mylist = []  # Acá vamos a almacenar los (x_i - μ), x_i i-ésima imágen (en filas, luego tendremos que transponerlos para que nos queden como columnas)
    while j < n:
        Mu = np.full(len(X), res[j])  # Vector tamaño M con el μ_j repetido M veces
        x_j = X[:, j] - Mu
        mylist.append(x_j)
        j += 1

    X_c = np.stack(mylist)  # Transformamos la lista en un array de NumPy (nos queda de n x m)
    X_c = X_c.T  # Como guardamos las columnas 1 x 1, se almacenaban como filas. Las queremos como columnas.

    return X_c


def mat_eigenvectors(filename):
    V = []
    f = open(filename, "r")
    for line in f.readlines():
        G = []
        for nro in line.split():
            G.append(float(nro))
        V.append(G)
    V = np.stack(V)  # Los autovectores vienen ordenados de mayor a menor, dependiendo del módulo de sus autovalores asociados.


    return V


def mat_eigenvectors_k(filename, k):
    f = open(filename, "r")
    V_k = []
    # queremos tomar los primeros k autovectores de la matriz de covarianza
    for linea in f.readlines():
        i = 0
        row = []
        for nro in linea.split():
            if i < k:
                row.append(float(nro))
                i += 1
        V_k.append(row)
    V_k = np.stack(V_k)

    return V_k
