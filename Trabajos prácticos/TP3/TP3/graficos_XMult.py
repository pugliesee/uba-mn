import numpy as np
import matplotlib.pyplot as plt

iteraciones = [10, 13, 15, 17, 20, 30, 40, 50, 100]


def grafico_J_Mat_XMult():
    errores = []
    x = np.arange(len(iteraciones))
    width = 0.15
    f = open("Error_XMult_J_Mat.txt", "r")
    for line in f.readlines():
        e = []
        for nro in line.split():
            e.append(float(nro))
        errores.append(e)

    ax = plt.subplot()
    ax.bar(x - 1.5 * width, errores[0], width, label="x_init * 2")
    ax.bar(x - 0.5 * width, errores[1], width, label="x_init * 10")
    ax.bar(x + 0.5 * width, errores[2], width, label="x_init * 100")
    ax.bar(x + 1.5 * width, errores[3], width, label="x_init * 1000")
    ax.set_title("Comparación de error para Jacobi Matricial para \n distintos valores iniciales de x")
    ax.set_xticks(x)
    ax.set_xticklabels(iteraciones)
    ax.set_xlabel("Cantidad de Iteraciones")
    ax.set_ylabel("Error")
    plt.yscale("log")
    ax.legend()
    plt.show()


def grafico_J_Sum_XMult():
    errores = []
    x = np.arange(len(iteraciones))
    width = 0.15
    f = open("Error_XMult_J_Sum.txt", "r")
    for line in f.readlines():
        e = []
        for nro in line.split():
            e.append(float(nro))
        errores.append(e)

    ax = plt.subplot()
    ax.bar(x - 1.5 * width, errores[0], width, label="x_init * 2")
    ax.bar(x - 0.5 * width, errores[1], width, label="x_init * 10")
    ax.bar(x + 0.5 * width, errores[2], width, label="x_init * 100")
    ax.bar(x + 1.5 * width, errores[3], width, label="x_init * 1000")
    ax.set_title("Comparación de error para Jacobi Sumatoria para \n distintos valores iniciales de x")
    ax.set_xticks(x)
    ax.set_xticklabels(iteraciones)
    ax.set_xlabel("Cantidad de Iteraciones")
    ax.set_ylabel("Error")
    plt.yscale("log")
    ax.legend()
    plt.show()


def grafico_GS_Mat_XMult():
    errores = []
    x = np.arange(len(iteraciones))
    width = 0.15
    f = open("Error_XMult_GS_Mat.txt", "r")
    for line in f.readlines():
        e = []
        for nro in line.split():
            e.append(float(nro))
        errores.append(e)

    ax = plt.subplot()
    ax.bar(x - 1.5 * width, errores[0], width, label="x_init * 2")
    ax.bar(x - 0.5 * width, errores[1], width, label="x_init * 10")
    ax.bar(x + 0.5 * width, errores[2], width, label="x_init * 100")
    ax.bar(x + 1.5 * width, errores[3], width, label="x_init * 1000")
    ax.set_title("Comparación de error para Gauss-Seidel Matricial para \n distintos valores iniciales de x")
    ax.set_xticks(x)
    ax.set_xticklabels(iteraciones)
    ax.set_xlabel("Cantidad de Iteraciones")
    ax.set_ylabel("Error")
    plt.yscale("log")
    ax.legend()
    plt.show()


def grafico_GS_Sum_XMult():
    errores = []
    x = np.arange(len(iteraciones))
    width = 0.15
    f = open("Error_XMult_GS_Sum.txt", "r")
    for line in f.readlines():
        e = []
        for nro in line.split():
            e.append(float(nro))
        errores.append(e)

    ax = plt.subplot()
    ax.bar(x - 1.5 * width, errores[0], width, label="x_init * 2")
    ax.bar(x - 0.5 * width, errores[1], width, label="x_init * 10")
    ax.bar(x + 0.5 * width, errores[2], width, label="x_init * 100")
    ax.bar(x + 1.5 * width, errores[3], width, label="x_init * 1000")
    ax.set_title("Comparación de error para Gauss-Seidel Sumatoria para \n distintos valores iniciales de x")
    ax.set_xticks(x)
    ax.set_xticklabels(iteraciones)
    ax.set_xlabel("Cantidad de Iteraciones")
    ax.set_ylabel("Error")
    plt.yscale("log")
    ax.legend()
    plt.show()


def generarGraficosXMult():
    grafico_J_Mat_XMult()
    grafico_J_Sum_XMult()
    grafico_GS_Mat_XMult()
    grafico_GS_Sum_XMult()
