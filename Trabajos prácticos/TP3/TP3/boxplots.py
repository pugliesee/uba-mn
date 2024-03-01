import numpy as np
import matplotlib.pyplot as plt


def generarBoxPlotsTiempos():
    sizes = [128, 256, 512]  # son los más interesantes de observar
    metodos = ["J Mat", "J Sum", "G-S Mat", "G-S Sum", "LU"]

    i = 0
    for tamano in sizes:
        tiempos = []
        f = open("TiempoPorItTamaño" + str(tamano) + ".txt")
        for line in f.readlines():
            linea = []
            for nro in line.split():
                linea.append(float(nro))
            tiempos.append(linea)

        plt.boxplot(tiempos, labels=metodos)
        plt.title("Tiempos de 10 ejecuciones con matrices EDD de " + str(tamano) + "x" + str(tamano))
        plt.ylabel("Tiempo (milisegundos)")
        plt.xlabel("Métodos")
        plt.yscale("log")
        plt.show()


