import numpy as np
import matplotlib.pyplot as plt


def graficar_particulares():
    # El archivo es de la siguiente forma:
    # FILA 1: EDD
    # FILA 2: IDENTIDAD
    # FILA 3: SDP
    # FILA 4: TRIANGULAR
    errores = []
    metodos = ["J Mat", "J Sum", "G-S Mat", "G-S Sum", "LU"]
    x = np.arange(len(metodos))
    width = 0.15
    f = open("Errores_Matrices_Particulares.txt", "r")
    for line in f.readlines():
        e = []
        for nro in line.split():
            e.append(float(nro))
        errores.append(e)

    ax = plt.subplot()
    ax.bar(x - 1.5 * width, errores[0], width, label="EDD")
    ax.bar(x - 0.5 * width, errores[1], width, label="Identidad")
    ax.bar(x + 0.5 * width, errores[2], width, label="SDP")
    ax.bar(x + 1.5 * width, errores[3], width, label="Triangular")
    ax.set_title("Error de ejecución de Métodos Iterativos y LU sin iteraciones fijas \n para matrices particulares")
    ax.set_xticks(x)
    ax.set_xticklabels(metodos)
    ax.set_xlabel("Método")
    ax.set_ylabel("Error")
    plt.yscale("log")
    plt.legend()
    plt.show()
