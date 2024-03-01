import numpy as np
import matplotlib.pyplot as plt


def generarGraficosErrores():
    sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    iteraciones = [10, 50, 100, 300, 500, 1000, 2500]
    x = np.arange(len(iteraciones))
    width = 0.20

    # los errores los voy a leer como una matriz de 5 x len(iteraciones)
    # los leo de su archivo correspondiente para cada tamaño de matriz
    # la distribucion de los datos en cada matriz va a ser de la siguiente forma:
    # FILA 1: Jacobi Matricial
    # FILA 2: Jacobi Sumatoria
    # FILA 3: Gauss-Seidel Matricial
    # FILA 4: Gauss-Seidel Sumatoria
    # FILA 5 = LU (PERO ESTA NO LA USAMOS EN ESTE GRÁFICO)

    for tamano in sizes:
        errores = []
        # este gráfico va a ser de barras
        f = open("ErroresTamaño" + str(tamano) + ".txt", "r")
        for line in f.readlines():
            e = []
            for nro in line.split():
                e.append(float(nro))
            errores.append(e)

        ax = plt.subplot()
        ax.bar(x - 1.5*width, errores[0], width, label="Jacobi Matricial")
        ax.bar(x - 0.5*width, errores[1], width, label="Jacobi Sumatoria")
        ax.bar(x + 0.5*width, errores[2], width, label="G-S Matricial")
        ax.bar(x + 1.5*width, errores[3], width, label="G-S Sumatoria")
        ax.set_title("Error promedio de métodos iterativos para \nmatrices EDD de " + str(tamano) + "x" + str(tamano))
        ax.set_xticks(x)
        ax.set_xticklabels(iteraciones)
        ax.set_xlabel("Cantidad de Iteraciones")
        ax.set_ylabel("Error")
        plt.yscale("log")
        ax.legend()
        plt.show()


def generarGraficosErroresLineas():
    sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    iteraciones = [10, 50, 100, 300, 500, 1000, 2500]

    # los errores los voy a leer como una matriz de 5 x len(iteraciones)
    # los leo de su archivo correspondiente para cada tamaño de matriz
    # la distribucion de los datos en cada matriz va a ser de la siguiente forma:
    # FILA 1: Jacobi Matricial
    # FILA 2: Jacobi Sumatoria
    # FILA 3: Gauss-Seidel Matricial
    # FILA 4: Gauss-Seidel Sumatoria
    # FILA 5 = LU (lo mostramos constante como una linea discontinua pues no depende de una cantidad de iteraciones)

    for tamano in sizes:
        errores = []
        # este gráfico va a ser de lineas
        f = open("ErroresTamaño" + str(tamano) + ".txt", "r")
        for line in f.readlines():
            e = []
            for nro in line.split():
                e.append(float(nro))
            errores.append(e)

        promedio_LU = 0
        i = 0
        while i < len(errores[4]):
            promedio_LU += errores[4][i] / len(errores[4])
            i += 1

        i = 0
        while i < len(errores[4]):
            errores[4][i] = promedio_LU
            i += 1

        plt.plot(iteraciones, errores[0], marker="o", label="Jacobi Matricial")
        plt.plot(iteraciones, errores[1], marker="o", label="Jacobi Sumatoria")
        plt.plot(iteraciones, errores[2], marker="o", label="G-S Matricial")
        plt.plot(iteraciones, errores[3], marker="o", label="G-S Sumatoria")
        plt.plot(iteraciones, errores[4], "--", label="LU")
        plt.yscale("log")
        plt.xscale("log")
        plt.ylabel("Error")
        plt.xlabel("Cantidad de Iteraciones")
        plt.title("Error promedio de métodos iterativos para \nmatrices EDD de " + str(tamano) + "x" + str(tamano))
        plt.legend(loc="best")
        plt.show()


def generarGraficosTiemposXSize():
    sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    tiempos = []
    # este gráfico es de lineas
    f = open("TimesMetodosXTamaño.txt", "r")
    for line in f.readlines():
        t = []
        for nro in line.split():
            t.append(float(nro))
        tiempos.append(t)

    plt.plot(sizes, tiempos[0], marker="o", label="Jacobi Matricial")
    plt.plot(sizes, tiempos[1], marker="o", label="Jacobi Sumatoria")
    plt.plot(sizes, tiempos[2], marker="o", label="G-S Matricial")
    plt.plot(sizes, tiempos[3], marker="o", label="G-S Sumatoria")
    plt.plot(sizes, tiempos[4], marker="o", label="LU")
    plt.title("Tiempo de convergencia promedio de métodos iterativos y LU para \nmatrices EDD según su tamaño")
    plt.xlabel("Tamaño de la matriz")
    plt.ylabel("Tiempo")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.show()


def generarGraficosTiemposSegunIteraciones():

    # los errores los voy a leer como una matriz de 5 x len(iteraciones)
    # los leo de su archivo correspondiente para cada tamaño de matriz
    # la distribucion de los datos en cada matriz va a ser de la siguiente forma:
    # FILA 1: Jacobi Matricial
    # FILA 2: Jacobi Sumatoria
    # FILA 3: Gauss-Seidel Matricial
    # FILA 4: Gauss-Seidel Sumatoria
    # FILA 5 = LU (esta la vamos a mostrar constante, ya que no depende de una cantidad de iteraciones)
    sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    iteraciones = [10, 50, 100, 300, 500, 1000, 2500]

    for tamano in sizes:
        tiempos = []
        f = open("TimesTamaño" + str(tamano) + ".txt", "r")
        for linea in f.readlines():
            t = []
            for nro in linea.split():
                t.append(float(nro))
            tiempos.append(t)

        promedio_LU = 0
        i = 0
        while i < len(tiempos[4]):
            promedio_LU += tiempos[4][i] / len(tiempos[4])
            i += 1

        i = 0
        while i < len(tiempos[4]):
            tiempos[4][i] = promedio_LU
            i += 1

        plt.plot(iteraciones, tiempos[0], marker="o", label="Jacobi Matricial")
        plt.plot(iteraciones, tiempos[1], marker="o", label="Jacobi Sumatoria")
        plt.plot(iteraciones, tiempos[2], marker="o", label="G-S Matricial")
        plt.plot(iteraciones, tiempos[3], marker="o", label="G-S Sumatoria")
        plt.plot(iteraciones, tiempos[4], "--", label="LU")
        plt.yscale("log")
        plt.ylabel("Tiempo")
        plt.xlabel("Cantidad de Iteraciones")
        plt.legend(loc="best")
        plt.title("Tiempo de ejecución promedio de métodos iterativos y LU \npara matrices de " + str(tamano) + "x" + str(tamano))
        plt.show()
