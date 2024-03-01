import numpy as np
import matplotlib.pyplot as plt


def ECM(vec_eigen, vec_propio):
    n = len(vec_propio)
    ecm = 0
    i = 0
    while i < n:
        ecm += (vec_eigen[i] - vec_propio[i])**2
        i += 1
    ecm = ecm / n

    return ecm


def imprimir_graficos_tiempos():
    nro_autovalor = ["λ1", "λ2", "λ3", "λ4", "λ5", "λ6", "λ7", "λ8", "λ9", "λ10"]
    # Vamos primero con el gráfico de los tests 1 a 4
    tiempos_t1 = []
    tiempos_t2 = []
    tiempos_t3 = []
    tiempos_t4 = []
    tiempos_promedio = []
    f = open("TiemposTest1.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t1.append(float(nro))

    f = open("TiemposTest2.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t2.append(float(nro))

    f = open("TiemposTest3.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t3.append(float(nro))

    f = open("TiemposTest4.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t4.append(float(nro))

    f = open("TiemposPromedio.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_promedio.append(float(nro))

    plt.plot(nro_autovalor, tiempos_t1, marker="o", label="Identidad")
    plt.plot(nro_autovalor, tiempos_t2, marker="o", label="Simétrica")
    plt.plot(nro_autovalor, tiempos_t3, marker="o", label="Triangular")
    plt.plot(nro_autovalor, tiempos_t4, marker="o", label="Aleatoria")
    plt.plot(nro_autovalor, tiempos_promedio, marker="o", label="Promedio")
    plt.title("Tiempos de convergencia del método de la potencia con deflación \n para autovectores de cada autovalor de matrices particulares en ns")
    plt.xlabel("Autovalor asociado al autovector")
    plt.ylabel("Tiempo")
    plt.legend(loc="best")
    plt.yscale("log")
    plt.show()

    # Vamos ahora con el gráfico de los test 5 a 8
    tiempos_t5 = []
    tiempos_t6 = []
    tiempos_t7 = []
    tiempos_t8 = []
    f = open("TiemposTest5.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t5.append(float(nro))

    f = open("TiemposTest6.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t6.append(float(nro))

    f = open("TiemposTest7.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t7.append(float(nro))

    f = open("TiemposTest8.txt", "r")
    for line in f.readlines():
        for nro in line.split():
            tiempos_t8.append(float(nro))

    plt.plot(nro_autovalor, tiempos_t5, marker="o", label="Cercanos Grande")
    plt.plot(nro_autovalor, tiempos_t6, marker="o", label="Cercanos Chico")
    plt.plot(nro_autovalor, tiempos_t7, marker="o", label="5 Grandes, 5 Chicos")
    plt.plot(nro_autovalor, tiempos_t8, marker="o", label="Saltos de a 1000")
    plt.title("Tiempos de convergencia del método de la potencia con deflación \n en ns para matrices según relación entre módulos de autovalores")
    plt.xlabel("Autovalor asociado al autovector")
    plt.ylabel("Tiempo")
    plt.legend(loc="best")
    plt.yscale("log")
    plt.show()


def graficar_errores_epsilon():
    nro_autovalor = ["λ1", "λ2", "λ3", "λ4", "λ5", "λ6", "λ7", "λ8", "λ9", "λ10"]
    errores_e12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    errores_e6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    t = 5
    while t < 9:
        mat_eigen = []
        mat_e12 = []
        mat_e6 = []
        f = open("EigenvectorsTest" + str(t) + ".txt", "r")
        for line in f.readlines():
            fila = []
            for nro in line.split():
                fila.append(float(nro))
            mat_eigen.append(fila)
        mat_eigen = np.stack(mat_eigen)

        f = open("AutovectoresTest" + str(t) + "E12.txt", "r")
        for line in f.readlines():
            fila = []
            for nro in line.split():
                fila.append(float(nro))
            mat_e12.append(fila)
        mat_e12 = np.stack(mat_e12)

        f = open("AutovectoresTest" + str(t) + "E6.txt", "r")
        for line in f.readlines():
            fila = []
            for nro in line.split():
                fila.append(float(nro))
            mat_e6.append(fila)
        mat_e6 = np.stack(mat_e6)

        j = 0
        while j < 10:
            eigenvec = mat_eigen[:, j]
            autovec_e12 = mat_e12[:, j]
            autovec_e6 = mat_e6[:, j]

            errores_e12[j] += ECM(eigenvec, autovec_e12) / 4  # dividimos por 4 porque es el promedio entre los 4 tests
            errores_e6[j] += ECM(eigenvec, autovec_e6) / 4
            j += 1

        t += 1

    plt.plot(nro_autovalor, errores_e12, marker="o", label="ε = 10⁻¹²")
    plt.plot(nro_autovalor, errores_e6, marker="o", label="ε = 10⁻⁶")
    plt.title("Error promedio en cálculo de autovectores de método de la potencia \n con deflación comparado con Eigen")
    plt.xlabel("Autovalor asociado al autovector")
    plt.ylabel("Error")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.show()

    print(errores_e12)
    print(errores_e6)


def graficar_error_iteraciones():
    iteraciones = [10, 20, 30, 40, 50, 80, 100, 150, 200, 300, 500]
    errores = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    mat_eigen = []

    f = open("EigenvectorsTest5.txt", "r")
    for line in f.readlines():
        fila = []
        for nro in line.split():
            fila.append(float(nro))
        mat_eigen.append(fila)
    mat_eigen = np.stack(mat_eigen)

    i = 0
    while i < len(iteraciones):
        mat_iters = []
        f = open("AutovectoresTest5" + str(iteraciones[i]) + "Iters.txt", "r")
        for line in f.readlines():
            fila = []
            for nro in line.split():
                fila.append(float(nro))
            mat_iters.append(fila)
        mat_iters = np.stack(mat_iters)

        j = 0
        while j < 10:
            eigenvec = mat_eigen[:, j]
            autovec_iters = mat_iters[:, j]
            errores[i] += ECM(eigenvec, autovec_iters) / 10  # divido por 10 pues en cada posicion ponemos el error promedio entre todos los autovectores para la cantidad de iteraciones dada

            j += 1

        i += 1

    plt.plot(iteraciones, errores, marker="o", color="orange")
    plt.title("Error de convergencia promedio para matriz diagonal según \n cantidad de iteraciones en método de la potencia con deflación")
    plt.xlabel("N° de iteraciones")
    plt.ylabel("Error")
    plt.show()
