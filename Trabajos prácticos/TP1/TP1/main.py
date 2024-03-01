# This is a sample Python script.

# Press Mayús+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lng
import math
from time import time


# EJERCICIO 1
def elim_gauss_sin_pivoteo(mat_a, b):
    n = len(mat_a)
    if n != len(b):
        print("Error, el sistema no posee solución")
        return
    # armamos la matriz extendida con los elementos de b
    j = 0
    while j < n:
        mat_a[j].append(b[j])
        j = j + 1
    # comienza la eliminación Gaussiana
    i = 0
    while i < n-1:
        j = i + 1
        while j < n:
            # como este algoritmo es sin pivoteo, si algún elemento de la diagonal es 0, se rompe.
            if mat_a[i][i] == 0:
                print("Error, el algoritmo no puede resolver este sistema")
                return
            m = -mat_a[j][i]/mat_a[i][i]
            k = i
            # acá no usamos n ya que también hay que modificar el valor de la columna con la que extendimos a A
            while k < len(mat_a[j]):
                mat_a[j][k] = mat_a[j][k] + m*mat_a[i][k]
                k = k + 1
            j = j + 1
        i = i + 1
    print("Así queda la matriz triangulada: ")
    for fila in mat_a:
        print(fila)
    res = []
    # con la matriz ya triangulada (si es que se pudo) primero analizamos el caso en el que hayan quedado todos ceros
    # en la última fila para ver si no existe solución o si existen infinitas
    if mat_a[n-1][n-1] == 0 and mat_a[n-1][n] == 0:
        print("El sistema posee infinitas soluciones")
        return
    elif mat_a[n-1][n-1] == 0 and mat_a[n-1][n] != 0:
        print("Error, el sistema no posee solución")
        return
    else:
        xn = mat_a[n-1][n]/mat_a[n-1][n-1]
        res.append(xn)
    # Luego continuamos con el proceso de backward substitution
    # Aquí no hace falta chequear si los elementos de la diagonal son 0, ya que si lo fueran,
    # el algoritmo ya se hubiera roto antes de llegar a esta instancia
    nro_fila = n-2
    while nro_fila >= 0:
        xi = mat_a[nro_fila][n]
        t = 0
        while t < len(res):
            xi = xi - res[t]*mat_a[nro_fila][n-t-1]
            t = t + 1
        xi = xi / mat_a[nro_fila][nro_fila]
        res.append(xi)
        nro_fila = nro_fila - 1

    # hace falta invertir el array ya que estábamos guardando el resultado como Xn,.....,X1
    res = res[::-1]
    print("El resultado es: ", res)
    return


# EJERCICIO 2
def elim_gauss_piv_parcial(mat_a, b):
    f = 0
    n = len(mat_a)
    if n != len(b):
        print("Error, el sistema no posee solución")
        return
    # armamos la matriz extendida con los elementos de b
    j = 0
    while j < n:
        mat_a[j] = np.append(mat_a[j], [b[j]])
        j = j + 1
    # comienza la eliminación Gaussiana
    i = 0
    while i < n - 1:
        # si el elemento en la diagonal es 0
        # (o capaz no es exactamente 0 por error numérico pero está muy cerca), recurrimos al pivoteo
        if abs(mat_a[i][i]) < 0.0000000001:
            # fila con la que voy a intercambiar
            p = i
            # auxiliar para calcular p
            v = i
            while v < n:
                if abs(mat_a[v][i]) > abs(mat_a[p][i]):
                    p = v
                v = v + 1
            t = i
            # aquí no usamos n ya que al hacer el pivoteo, también debemos intercambiar elementos de la columna b
            while t < len(mat_a[0]):
                mat_a[i][t], mat_a[p][t] = mat_a[p][t], mat_a[i][t]
                t = t + 1
        # una vez hecho el pivoteo, continuamos normalmente
        # si el elemento de la diagonal se mantuvo en 0, quiere decir que todos los elementos debajo de él en esa
        # columna también eran 0. Luego, pasamos al siguiente paso de la eliminación Gaussiana
        if mat_a[i][i] != 0:
            j = i + 1
            while j < n:
                m = -mat_a[j][i] / mat_a[i][i]
                k = i
                if abs(mat_a[i][i]) < 0.0000000001:
                    f = 1
                # acá no usamos n ya que también hay que modificar el valor de la columna con la que extendimos a A
                while k < len(mat_a[j]):
                    mat_a[j][k] = mat_a[j][k] + m * mat_a[i][k]
                    k = k + 1
                j = j + 1
        i = i + 1
    print("Así queda la matriz triangulada: ")
    for fila in mat_a:
        print(fila)
    res = []
    # analizamos primero la última fila para ver si ya podemos decir que el sistema no tiene solución,
    # o tiene infinitas, o si tenemos que continuar con el proceso de backward substitution
    if mat_a[n-1][n-1] == 0 and mat_a[n-1][n] == 0:
        print("El sistema posee infinitas soluciones")
        return
    elif mat_a[n-1][n-1] == 0 and mat_a[n-1][n] != 0:
        print("Error, el sistema no posee solución")
        return
    else:
        xn = mat_a[n-1][n]/mat_a[n-1][n-1]
        if abs(mat_a[n-1][n-1]) < 0.0001:
                    f = 1
        res.append(xn)
    # En este caso SÍ hace falta chequear para el resto de las filas si el elemento de la diagonal es 0.
    # Si así fuera, el sistema tendría infinitas soluciones
    nro_fila = n - 2
    while nro_fila >= 0:
        if mat_a[nro_fila][nro_fila] == 0:
            print("El sistema tiene infinitas soluciones")
            return
        # inicializamos xi en el valor que se encuentra en la fila i y columna n+1
        xi = mat_a[nro_fila][n]
        t = 0
        while t < len(res):
            xi = xi - res[t] * mat_a[nro_fila][n - t - 1]
            t = t + 1
        xi = xi / mat_a[nro_fila][nro_fila]
        if abs(mat_a[nro_fila][nro_fila]) < 0.0001:
            f = 1
        res.append(xi)
        nro_fila = nro_fila - 1
    # hace falta invertir el array ya que estábamos guardando el resultado como Xn,.....,X1
    res = res[::-1]
    print("El resultado es: ", res)
    if f == 1:
        print("El resultado puede contener errores de cálculo numérico")
    return

# EJERCICIO 2 Versión sin matriz extendida
def elim_gauss_piv_parcial2(mat_a, b):
    f = 0
    n = len(mat_a)
    if n != len(b):
        print("Error, el sistema no posee solución")
        return
    # comienza la eliminación Gaussiana
    i = 0
    while i < n - 1:
        # si el elemento en la diagonal es 0, recurrimos al pivoteo
        if mat_a[i][i] == 0:
            # fila con la que voy a intercambiar
            p = i
            # auxiliar para calcular p
            v = i
            while v < n:
                if abs(mat_a[v][i]) > abs(mat_a[p][i]):
                    p = v
                v = v + 1
            mat_a[i], mat_a[p] =  mat_a[p],  mat_a[i]
            b[i], b[p] = b[p], b[i]
        # una vez hecho el pivoteo, continuamos normalmente
        # si el elemento de la diagonal se mantuvo en 0, quiere decir que todos los elementos debajo de él en esa
        # columna también eran 0. Luego, pasamos al siguiente paso de la eliminación Gaussiana
        if mat_a[i][i] != 0:
            j = i + 1
            while j < n:
                m = -mat_a[j][i] / mat_a[i][i]
                k = i
                if mat_a[i][i] < 0.0001:
                    f = 1
                # acá no usamos n ya que también hay que modificar el valor de la columna con la que extendimos a A
                while k < n:
                    mat_a[j][k] = mat_a[j][k] + m * mat_a[i][k]
                    k = k + 1
                b[j] = b[j] + m * b[i]
                j = j + 1
        i = i + 1

    res = []
    # analizamos primero la última fila para ver si ya podemos decir que el sistema no tiene solución,
    # o tiene infinitas, o si tenemos que continuar con el proceso de backward substitution
    if mat_a[n-1][n-1] == 0 and b[n-1] == 0:
        print("El sistema posee infinitas soluciones")
        return
    elif mat_a[n-1][n-1] == 0 and b[n] != 0:
        print("Error, el sistema no posee solución")
        return
    else:
        xn = b[n-1]/mat_a[n-1][n-1]
        if mat_a[n-1][n-1] < 0.0001:
                    f = 1
        res.append(xn)
    # En este caso SÍ hace falta chequear para el resto de las filas si el elemento de la diagonal es 0.
    # Si así fuera, el sistema tendría infinitas soluciones
    nro_fila = n - 2
    while nro_fila >= 0:
        if mat_a[nro_fila][nro_fila] == 0:
            print("El sistema tiene infinitas soluciones")
            return
        # inicializamos xi en el valor que se encuentra en la fila i y columna n+1
        xi = b[nro_fila]
        t = 0
        while t < len(res):
            xi = xi - res[t] * mat_a[nro_fila][n - t - 1]
            t = t + 1
        xi = xi / mat_a[nro_fila][nro_fila]
        if mat_a[nro_fila][nro_fila] < 0.0001:
            f = 1
        res.append(xi)
        nro_fila = nro_fila - 1
    # hace falta invertir el array ya que estábamos guardando el resultado como Xn,.....,X1
    res = res[::-1]
    # print("El resultado es: ", res)
    # if f == 1:
        # print("El resultado puede contener errores de cálculo numérico")
    return


# EJERCICIO 3A
def EC_Tridiagonal(mat_a, b):
    n = len(mat_a)
    # armamos la matriz extendida con los elementos de b
    i = 0
    while i < n - 1:
        if mat_a[i][i] == 0 and mat_a[i+1][i] != 0:     # SI TENGO QUE PIVOTEAR
            mat_a[i], mat_a[i+1] = mat_a[i+1], mat_a[i]
        k = -(mat_a[i+1][i]/mat_a[i][i])
        j = i
        while j < i + 2:        # CICLA 3 VECES, pues j comienza siendo igual a i. Entonces este ciclo es O(3) = O(1)
            mat_a[i+1][j] = mat_a[i+1][j] + (mat_a[i][j] * k)
            j = j + 1
        b[i+1] = b[i+1] + (b[i] * k)
        i = i + 1
    j = 0
    while j < n:
        mat_a[j].append(b[j])
        j = j + 1

    print("Así queda la matriz triangulada: ")
    for fila in mat_a:
        print(fila)

    res = []
    # analizamos primero la última fila para ver si ya podemos decir que el sistema no tiene solución,
    # o tiene infinitas, o si tenemos que continuar con el proceso de backward substitution
    if mat_a[n - 1][n - 1] == 0 and mat_a[n - 1][n] == 0:
        print("El sistema posee infinitas soluciones")
        return
    elif mat_a[n - 1][n - 1] == 0 and mat_a[n - 1][n] != 0:
        print("Error, el sistema no posee solución")
        return
    else:
        xn = mat_a[n - 1][n] / mat_a[n - 1][n - 1]
        res.append(xn)
    # En este caso SÍ hace falta chequear para el resto de las filas si el elemento de la diagonal es 0.
    # Si así fuera, el sistema tendría infinitas soluciones
    nro_fila = n - 2
    while nro_fila >= 0:
        xi = (mat_a[nro_fila][n] - mat_a[nro_fila][nro_fila+1]*res[n-2-nro_fila])/mat_a[nro_fila][nro_fila]
        res.append(xi)
        nro_fila = nro_fila - 1
    # hace falta invertir el array ya que estábamos guardando el resultado como Xn,.....,X1
    res = res[::-1]
    print("El resultado es: ", res)
    return

# EJERCICIO 3A Versión 2 sin extendet
def EC_Tridiagonal2(mat_a, b):
    n = len(mat_a)
    # armamos la matriz extendida con los elementos de b
    i = 0
    while i < n - 1:
        if mat_a[i][i] == 0 and mat_a[i+1][i] != 0:     # SI TENGO QUE PIVOTEAR
            mat_a[i], mat_a[i+1] = mat_a[i+1], mat_a[i]
        k = -(mat_a[i+1][i]/mat_a[i][i])
        j = i
        while j < i + 2:        # CICLA 3 VECES, pues j comienza siendo igual a i. Entonces este ciclo es O(3) = O(1)
            mat_a[i+1][j] = mat_a[i+1][j] + (mat_a[i][j] * k)
            j = j + 1
        b[i+1] = b[i+1] + (b[i] * k)
        i = i + 1
    j = 0

    res = []
    # analizamos primero la última fila para ver si ya podemos decir que el sistema no tiene solución,
    # o tiene infinitas, o si tenemos que continuar con el proceso de backward substitution
    if mat_a[n - 1][n - 1] == 0 and b[n - 1] == 0:
        print("El sistema posee infinitas soluciones")
        return
    elif mat_a[n - 1][n - 1] == 0 and b[n - 1] != 0:
        print("Error, el sistema no posee solución")
        return
    else:
        xn = b[n - 1] / mat_a[n - 1][n - 1]
        res.append(xn)
    # En este caso SÍ hace falta chequear para el resto de las filas si el elemento de la diagonal es 0.
    # Si así fuera, el sistema tendría infinitas soluciones
    nro_fila = n - 2
    while nro_fila >= 0:
        xi = (b[nro_fila] - mat_a[nro_fila][nro_fila+1]*res[n-2-nro_fila])/mat_a[nro_fila][nro_fila]
        res.append(xi)
        nro_fila = nro_fila - 1
    # hace falta invertir el array ya que estábamos guardando el resultado como Xn,.....,X1
    res = res[::-1]
    # print("El resultado es: ", res)
    return


# EJERCICIO 3B
def elim_gauss_sist_tridiagonal_todojunto(a, b, c, d):
    n = len(d)  # asumimos que a, b, c y d tienen la misma longitud
    i = 0
    e = []
    # complejidad del ciclo: O(n)
    while i < n-1:
        if b[i] != 0:  # no hacemos pivoteo
            b[i + 1] = b[i + 1] - (a[i + 1] / b[i]) * c[i]
            d[i+1] = d[i+1] - (a[i+1]/b[i])*d[i]
            a[i + 1] = 0
            e.append(0)  # si hago pivoteo, en algún lugar tengo que guardar el c[i+1] que hago subir de fila
                         # para luego poder usarlo cuando tenga que calcular la solución del sistema
        elif b[i] == 0 and a[i+1] != 0:  # hacemos pivoteo y con eso ya es suficiente
            b[i], a[i+1] = a[i+1], b[i]
            c[i], b[i+1] = b[i+1], c[i]
            d[i], d[i+1] = d[i+1], d[i]
            e.append(c[i+1])
            c[i+1] = 0
        # si no se entró a ninguno de los casos del if, quiere decir que toda la columna ya está en 0
        # luego, podemos pssar al siguiente paso
        i = i + 1
    res = []
    if b[n-1] == 0 and d[n-1] == 0:
        print("El sistem tiene infinitas soluciones")
        return
    elif b[n-1] == 0 and d[n-1] != 0:
        print("El sistema no tiene solución")
        return
    else:
        xn = d[n-1]/b[n-1]
        res.append(xn)
    t = n-2
    while t >= 0:
        if b[t] == 0:
            print("El sistema tiene infinitas soluciones")
            return
        x = (d[t] - c[t]*res[n-2-t] - e[t]*res[n-2-t-1])/b[t]
        res.append(x)
        t = t - 1
    res = res[::-1]
    print("El resultado es: ", res)
    return res


# EJERCICIO 3C
# Idea: modificar los vectores a, b y c de modo que nos quede la factorización LU de A
def precomp_sist_tridiagonal(a, b, c):
    n = len(a)
    i = 0
    # triangulamos en complejidad O(n)
    # luego de este while,  quedaría como el vector L que después utilizamos para modificar los distintos d_i
    while i < n - 1:
        b[i + 1] = b[i + 1] - (a[i + 1] / b[i]) * c[i]
        a[i + 1] = a[i+1]/b[i]
        i = i + 1
    return


# con esta función resolvemos el sistema con el precómputo ya hecho para un d determinado
def elim_gauss_sist_tridiagonal_precomp(a, b, c, d):
    n = len(d)
    i = 1
    # Modificamos d en base a nuestro vector L, que quedó guardado en a
    while i < n:
        d[i] = d[i] - a[i]*d[i-1]
        i = i + 1
    # Ahora resolvemos
    res=[]
    if abs(b[n - 1]) < 0.0001 and abs(d[n - 1]) < 0:
        print("El sistem tiene infinitas soluciones")
        return
    elif abs(b[n - 1]) < 0.0001 and abs(d[n - 1]) > 0.0001:
        print("El sistema no tiene solución")
        return
    else:
        xn = d[n - 1] / b[n - 1]
        res.append(xn)
    # Hacemos Backward Substitution
    t = n - 2
    while t >= 0:
        if abs(b[t]) < 0.0001:
            print("El sistema tiene infinitas soluciones")
            return
        x = (d[t] - c[t] * res[n - 2 - t]) / b[t]
        res.append(x)
        t = t - 1
    res = res[::-1]
    # print("El resultado es: ", res)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    bEjemplo1 = [1, 2, 3]
    bEjemplo2 = [3, 7, 6]
    bEjemplo3 = [4, 9, 2]
    bEjemplo4 = [4, 7, 4]

    testDiagonal = [[1, 2, 0, 0], [1, 6, 3, 0], [0, 5, 5, 1], [0, 0, 1, 1]]
    bDiagonal = [1, 1, 5, 10]

    testDiagonalExt = [[1, 2, 0, 0, 0, 0], [1, 6, 3, 0, 0, 0], [0, 5, 5, 1, 0, 0], [0, 0, 1, 1, 3, 0], [0, 0, 0, 7, 3, 2], [0, 0, 0, 0, 5, 9]]
    bDiagonalExt = [1, 1, 5, 10, 4, 7]

    a = [0, 1, 5, 1]
    b = [1, 6, 5, 1]
    c = [2, 3, 1, 0]

    aExt = [0, 1, 5, 1, 7, 5]
    bExt = [1, 6, 5, 1, 3, 9]
    cExt = [2, 3, 1, 3, 2, 0]

    # TESTS EJERCICIO 1
    matEj1_con_solucion = [[3, 7, 5], [5, 3, 8], [7, 2, 9]]
    # para este ejemplo debería fallar ya que sería necesario realizar pivoteo
    matEj1_requiere_piv = [[0, 1, 1], [3, 7, 2], [1, 5, 9]]
    matEj1_inf_soluciones = [[1, 4, 7], [5, 3, 2], [2, 8, 14]]
    print("EJERCICIO 1 - EJEMPLO 1, debería funcionar sin problema")
    elim_gauss_sin_pivoteo(matEj1_con_solucion, bEjemplo1)
    print("----------")
    print("EJERCICIO 1 - EJEMPLO 2, debería fallar ya que requiere pivoteo")
    elim_gauss_sin_pivoteo(matEj1_requiere_piv, bEjemplo1)
    print("----------")
    print("EJERCICIO 1- EJEMPLO 3, se puede hacer la eliminación pero hay infinitas soluciones")
    elim_gauss_sin_pivoteo(matEj1_inf_soluciones, bEjemplo2)
    print("----------")

    # TESTS EJERCICIO 2
    matEj2_con_solucion = [[0, 1, 1], [3, 7, 2], [1, 5, 9]]
    matEj2_inf_soluciones = [[6, 2, 8], [7, 1, 5], [3, 1, 4]]
    matEj2_sin_solucion = [[4, 2, 6], [5, 6, 1], [8, 4, 12]]
    print("EJERCICIO 2 - EJEMPLO 1, debería funcionar sin problema")
    elim_gauss_piv_parcial(matEj2_con_solucion, bEjemplo1)
    print("----------")
    print("EJERCICIO 2 - EJEMPLO 2, debería tener infinitas soluciones")
    elim_gauss_piv_parcial(matEj2_inf_soluciones, bEjemplo3)
    print("----------")
    print("EJERCICIO 2 - EJEMPLO 3, no debería tener solución")
    elim_gauss_piv_parcial(matEj2_sin_solucion, bEjemplo4)
    print("----------")
    print("EJERCICIO 2 - EJEMPLO 4, debería contemplar el error numérico")
    matEj4_con_solucion = [[0.00000000002, 30, 1], [0.00000000004, 7, 2], [0.00000000005, 5, 9]]
    elim_gauss_piv_parcial(matEj4_con_solucion, bEjemplo1)
    print("----------")
    # TESTS EJERCICIO 3
    print("EJERCICIO 3A - EJEMPLO 1, Tiene solución")
    EC_Tridiagonal(testDiagonal, bDiagonal)
    print("----------")
    print("EJERCICIO 3A - EJEMPLO 2, Tiene solución")
    EC_Tridiagonal(testDiagonalExt, bDiagonalExt)
    print("----------")
    print("EJERCICIO 3B - EJEMPLO 1, Tiene solución")
    d = [1, 1, 5, 10]
    elim_gauss_sist_tridiagonal_todojunto(a, b, c, d)
    print("----------")
    print("EJERCICIO 3B - EJEMPLO 2, Tiene solución")

    dExt= [1, 1, 5, 10, 4, 7]
    elim_gauss_sist_tridiagonal_todojunto(aExt, bExt, cExt, dExt)
    print("----------")

    a = [0, 1, 5, 1]
    b = [1, 6, 5, 1]
    c = [2, 3, 1, 0]
    d = [1, 1, 5, 10]
    print("EJERCICIO 3C - EJEMPLO 1, Tiene solución")
    precomp_sist_tridiagonal(a, b, c)
    elim_gauss_sist_tridiagonal_precomp(a, b, c, d)
    print("----------")

    a = [0, 2, 12, 1]
    b = [1, 6, 5, 1]
    c = [2, 3, 1, 0]
    d = [1, 1, 5, 10]
    print("EJERCICIO 3C - EJEMPLO 2, ?????????")
    precomp_sist_tridiagonal(a, b, c)
    elim_gauss_sist_tridiagonal_precomp(a, b, c, d)
    print("----------")

    # GRÁFICOS EJERCICIO 5 :)
    # Gráfico 5A
    i = 0
    a = []
    b = []
    c = []
    d = []
    while i < 101:
        a.append(1)
        b.append(-2)
        c.append(1)
        if i == 51:
            d.append(4 / 101)
        else:
            d.append(0)
        i = i + 1

    resa = elim_gauss_sist_tridiagonal_todojunto(a, b, c, d)

    # Gráfico 5B
    i = 0
    a = []
    b = []
    c = []
    d = []
    while i < 101:
        a.append(1)
        b.append(-2)
        c.append(1)
        d.append(4 / (101 * 101))
        i = i + 1

    resb = elim_gauss_sist_tridiagonal_todojunto(a, b, c, d)

    # Gráfico 5C
    i = 0
    a = []
    b = []
    c = []
    d = []
    while i < 101:
        a.append(1)
        b.append(-2)
        c.append(1)
        d.append((-1 + 2 * i / 100) * 12 / (101 * 101))
        i = i + 1

    resc = elim_gauss_sist_tridiagonal_todojunto(a, b, c, d)

    x = np.arange(101)
    plt.plot(x, resa, color="blue", label='(a)')
    plt.plot(x, resb, color="orange", label='(b)')
    plt.plot(x, resc, color="green", label='(c)')
    plt.legend(loc='lower left')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.show()  # función especial para que muestre el gráfico en la consola de Pycharm :)

    # EJERCICIO 6B
    alpha = 1
    n = 101
    r = 10
    k = math.trunc(n / 2)
    A = np.diag(np.full(n - 1, -1), 1) + np.diag(np.full(n, 3)) + np.diag(np.full(n - 1, -1), -1)
    u = np.concatenate((np.full(k - r + 1, 0), np.full(k + r - (k - r + 1), 1), np.full(k - r + 1, 0)))
    m = 0
    arr_resultados = []
    while m < 1000:
        arr_resultados.append(u)
        u = lng.solve(A, u)
        m = m + 1

    arr_resultados = np.transpose(arr_resultados)

    plt.pcolor(arr_resultados)
    plt.colorbar()
    plt.xlabel('k')
    plt.ylabel('X')
    plt.show()

# EJERCICIO 4A
print("EJERCICIO 4A - Cálculo de tiempos")
print("-------------")
print("Para n = 5")
minimos_piv = []
minimos_tri = []
m = 0
s1 = 0
s2 = 0
tiempo_min_piv = 0
tiempo_min_tri = 0
set_t_min_piv = 0
set_t_min_tri = 0
while m < 100:  
    n = 5
    a = np.random.rand(n-1)
    b = np.random.rand(n)
    c = np.random.rand(n-1)
    d = np.random.rand(n)
    A_1 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    A_11 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    tiempo_inicial = time()
    elim_gauss_piv_parcial2(A_1, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s1 = s1 + tiempo_total
    if set_t_min_piv == 0:
        tiempo_min_piv = tiempo_total
        set_t_min_piv = 1
    elif tiempo_total < tiempo_min_piv:
        tiempo_min_piv = tiempo_total
    tiempo_inicial = time()
    EC_Tridiagonal2(A_11, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s2 = s2 + tiempo_total
    if set_t_min_tri == 0:
        tiempo_min_tri = tiempo_total
        set_t_min_tri = 1
    elif tiempo_total < tiempo_min_tri:
        tiempo_min_tri = tiempo_total
    m = m+1
minimos_piv.append(tiempo_min_piv)
minimos_tri.append(tiempo_min_tri)
print("Tiempo para Eliminacion con pivoteo: ", s1/100)   
print("Tiempo para Eliminación con tridiagonal: ", s2/100)
print("Tiempo mínimo para Eliminación con pivoteo: ", tiempo_min_piv)
print("Tiempo mínimo para Eliminación con tridiagonal: ", tiempo_min_tri)


print("-------------")
print("Para n = 10")
m = 0
s1 = 0
s2 = 0
n = 10
tiempo_min_piv = 0
tiempo_min_tri = 0
set_t_min_piv = 0
set_t_min_tri = 0
while m < 100:  
    a = np.random.rand(n-1)
    b = np.random.rand(n)
    c = np.random.rand(n-1)
    d = np.random.rand(n)
    A_1 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    A_11 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    tiempo_inicial = time()
    elim_gauss_piv_parcial2(A_1, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s1 = s1 + tiempo_total
    if set_t_min_piv == 0:
        tiempo_min_piv = tiempo_total
        set_t_min_piv = 1
    elif tiempo_total < tiempo_min_piv:
        tiempo_min_piv = tiempo_total
    tiempo_inicial = time()
    EC_Tridiagonal2(A_11, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s2 = s2 + tiempo_total
    if set_t_min_tri == 0:
        tiempo_min_tri = tiempo_total
        set_t_min_tri = 1
    elif tiempo_total < tiempo_min_tri:
        tiempo_min_tri = tiempo_total
    m = m+1
minimos_piv.append(tiempo_min_piv)
minimos_tri.append(tiempo_min_tri)
print("Tiempo para Eliminacion con pivoteo: ", s1/100)   
print("Tiempo para Eliminación con tridiagonal: ", s2/100)
print("Tiempo mínimo para Eliminación con pivoteo: ", tiempo_min_piv)
print("Tiempo mínimo para Eliminación con tridiagonal: ", tiempo_min_tri)


print("-------------")
print("Para n = 100")
m = 0
s1 = 0
s2 = 0
n = 100
tiempo_min_piv = 0
tiempo_min_tri = 0
set_t_min_piv = 0
set_t_min_tri = 0
while m < 0:
    a = np.random.rand(n-1)
    b = np.random.rand(n)
    c = np.random.rand(n-1)
    d = np.random.rand(n)
    A_1 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    A_11 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    tiempo_inicial = time()
    elim_gauss_piv_parcial2(A_1, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s1 = s1 + tiempo_total
    if set_t_min_piv == 0:
        tiempo_min_piv = tiempo_total
        set_t_min_piv = 1
    elif tiempo_total < tiempo_min_piv:
        tiempo_min_piv = tiempo_total
    tiempo_inicial = time()
    EC_Tridiagonal2(A_11, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s2 = s2 + tiempo_total
    if set_t_min_tri == 0:
        tiempo_min_tri = tiempo_total
        set_t_min_tri = 1
    elif tiempo_total < tiempo_min_tri:
        tiempo_min_tri = tiempo_total
    m = m+1
minimos_piv.append(tiempo_min_piv)
minimos_tri.append(tiempo_min_tri)
print("Tiempo para Eliminacion con pivoteo: ", s1/100)   
print("Tiempo para Eliminación con tridiagonal: ", s2/100)
print("Tiempo mínimo para Eliminación con pivoteo: ", tiempo_min_piv)
print("Tiempo mínimo para Eliminación con tridiagonal: ", tiempo_min_tri)

print("-------------")
print("Para n = 500")
m = 0
s1 = 0
s2 = 0
n = 500
tiempo_min_piv = 0
tiempo_min_tri = 0
set_t_min_piv = 0
set_t_min_tri = 0
while m < 0:
    a = np.random.rand(n-1)
    b = np.random.rand(n)
    c = np.random.rand(n-1)
    d = np.random.rand(n)
    A_1 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    A_11 = np.diag(a, k=-1) + np.diag(b, k=0) + np.diag(c, k=1)
    tiempo_inicial = time()
    elim_gauss_piv_parcial2(A_1, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s1 = s1 + tiempo_total
    if set_t_min_piv == 0:
        tiempo_min_piv = tiempo_total
        set_t_min_piv = 1
    elif tiempo_total < tiempo_min_piv:
        tiempo_min_piv = tiempo_total
    tiempo_inicial = time()
    EC_Tridiagonal2(A_11, d)
    tiempo_final = time()
    tiempo_total = tiempo_final - tiempo_inicial
    s2 = s2 + tiempo_total
    if set_t_min_tri == 0:
        tiempo_min_tri = tiempo_total
        set_t_min_tri = 1
    elif tiempo_total < tiempo_min_tri:
        tiempo_min_tri = tiempo_total
    print(m)
    m = m+1
minimos_piv.append(tiempo_min_piv)
minimos_tri.append(tiempo_min_tri)
print("Tiempo para Eliminacion con pivoteo: ", s1/100)   
print("Tiempo para Eliminación con tridiagonal: ", s2/100)
print("Tiempo mínimo para Eliminación con pivoteo: ", tiempo_min_piv)
print("Tiempo mínimo para Eliminación con tridiagonal: ", tiempo_min_tri)

ns = ['5', '10', '100', '500']
tiempos_eg_pivoteo = [0.0002026510238647461, 0.0012599968910217286, 0.4737288689613342, 18.37085790157318]
tiempos_eg_tridiagonal = [6.180524826049805e-05, 0.00013384819030761718, 0.0007484269142150879, 0.0036716341972351076]

x = np.arange(len(ns))
width = 0.35

ax = plt.subplot()
ax.bar(x - width/2, tiempos_eg_pivoteo, width, label='Elim. Gaussiana normal con pivoteo', color="green")
ax.bar(x + width/2, tiempos_eg_tridiagonal, width, label='Elim. Gaussiana para tridiagonal', color="red")
ax.set_title("Tiempos promedio de ejecución para algoritmos de resolución de \nsistemas tridiagonales con EG normal con pivoteo y EG para tridiagonales")
ax.set_xticks(x)
ax.set_xticklabels(ns)
ax.set_xlabel("Tamaño de la matriz")
ax.set_ylabel("Tiempo de ejecución")
ax.legend()
plt.savefig("ej_4a_green_red.png")
plt.show()

minimos_tri = [2.9802322387695312e-05, 5.9604644775390625e-05, 0.0006043910980224609, 0.0031137466430664062]
minimos_piv = [9.465217590332031e-05, 0.0005238056182861328, 0.3801295757293701, 14.462109327316284]
ax = plt.subplot()
ax.bar(x - width/2, minimos_piv, width, label='Elim. Gaussiana normal con pivoteo', color="pink")
ax.bar(x + width/2, minimos_tri, width, label='Elim. Gaussiana para tridiagonal', color="turquoise")
ax.set_title("Tiempos mínimos de ejecución para algoritmos de resolución de \nsistemas tridiagonales con EG normal con pivoteo y EG para tridiagonales")
ax.set_xticks(x)
ax.set_xticklabels(ns)
ax.set_xlabel("Tamaño de la matriz")
ax.set_ylabel("Tiempo de ejecución")
ax.legend()
plt.savefig("ej_4a_minimos.png")
plt.show()


# EJERCICIO 4B
print("EJERCICIO 4B - Cálculo de tiempos")
print("-------------")
print("Para n = 5")
# vamos a calcular el tiempo promedio para 10 matrices distintas
# para cada matriz, vamos a resolver el sistema (con y sin precómputo) unas 100 veces, y nos vamos a quedar con el tiempo promedio
# la matriz se va a mantener igual, lo que vamos a ir cambiando es el d para poder mostrar cómo se aprovecha el precómputo
minimos_precomp = []
minimos_sin_precomp = []
n_matrices = 0
s1_todas_las_matrices = 0
s2_todas_las_matrices = 0
tiempo_min_sin_precomp = 0
tiempo_min_con_precomp = 0
set_tmin_sin_precomp = 0
set_tmin_con_precomp = 0
while n_matrices < 10:
    m = 0
    s1 = 0
    s2 = 0
    n = 5
    a_sin_precomp = np.insert(np.random.rand(n-1), 0, 0)
    b_sin_precomp = np.random.rand(n)
    c_sin_precomp = np.append(np.random.rand(n-1), 0)
    a_con_precomp = a_sin_precomp.copy()
    b_con_precomp = b_sin_precomp.copy()
    c_con_precomp = c_sin_precomp.copy()
    t_init_precomp = time()
    precomp_sist_tridiagonal(a_con_precomp, b_con_precomp, c_con_precomp)
    t_end_precomp = time()
    tiempo_min_con_precomp = t_end_precomp - t_init_precomp
    while m < 100:
        d = np.random.rand(n)
        A_1 = np.diag(a_sin_precomp[1:n].copy(), k=-1) + np.diag(b_sin_precomp.copy(), k=0) + np.diag(c_sin_precomp[0:(n-1)].copy(), k=1)
        tiempo_inicial = time()
        elim_gauss_sist_tridiagonal_precomp(a_con_precomp, b_con_precomp, c_con_precomp, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s1 = s1 + tiempo_total
        if set_tmin_con_precomp == 0:
            tiempo_min_con_precomp = tiempo_min_con_precomp + tiempo_total
            set_tmin_con_precomp = 1
        elif tiempo_total < tiempo_min_con_precomp:
            tiempo_min_con_precomp = tiempo_total
        tiempo_inicial = time()
        EC_Tridiagonal2(A_1, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s2 = s2 + tiempo_total
        if set_tmin_sin_precomp == 0:
            tiempo_min_sin_precomp = tiempo_total
            set_tmin_sin_precomp = 1
        elif tiempo_total < tiempo_min_sin_precomp:
            tiempo_min_sin_precomp = tiempo_total
        m = m+1
    s1 = s1 + (t_end_precomp - t_init_precomp)
    s1_todas_las_matrices = s1_todas_las_matrices + s1/100
    s2_todas_las_matrices = s2_todas_las_matrices + s2/100
    n_matrices = n_matrices + 1
minimos_precomp.append(tiempo_min_con_precomp)
minimos_sin_precomp.append(tiempo_min_sin_precomp)
print("Tiempo promedio para Eliminacion con precómputo: ", s1_todas_las_matrices/10)
print("Tiempo promedio para Eliminación sin precómputo: ", s2_todas_las_matrices/10)
print("Tiempo mínimo para Eliminación con precómputo: ", tiempo_min_con_precomp)
print("Tiempo mínimo para Eliminación sin precómputo: ", tiempo_min_sin_precomp)


print("-------------")
print("Para n = 10")
# vamos a calcular el tiempo promedio para 10 matrices distintas
# para cada matriz, vamos a resolver el sistema (con y sin precómputo) unas 100 veces, y nos vamos a quedar con el tiempo promedio
# la matriz se va a mantener igual, lo que vamos a ir cambiando es el d para poder mostrar cómo se aprovecha el precómputo
n_matrices = 0
s1_todas_las_matrices = 0
s2_todas_las_matrices = 0
tiempo_min_sin_precomp = 0
tiempo_min_con_precomp = 0
set_tmin_sin_precomp = 0
set_tmin_con_precomp = 0
while n_matrices < 10:
    m = 0
    s1 = 0
    s2 = 0
    n = 10
    a_sin_precomp = np.insert(np.random.rand(n-1), 0, 0)
    b_sin_precomp = np.random.rand(n)
    c_sin_precomp = np.append(np.random.rand(n-1), 0)
    a_con_precomp = a_sin_precomp.copy()
    b_con_precomp = b_sin_precomp.copy()
    c_con_precomp = c_sin_precomp.copy()
    t_init_precomp = time()
    precomp_sist_tridiagonal(a_con_precomp, b_con_precomp, c_con_precomp)
    t_end_precomp = time()
    while m < 100:
        d = np.random.rand(n)
        A_1 = np.diag(a_sin_precomp[1:n].copy(), k=-1) + np.diag(b_sin_precomp.copy(), k=0) + np.diag(c_sin_precomp[0:(n-1)].copy(), k=1)
        tiempo_inicial = time()
        elim_gauss_sist_tridiagonal_precomp(a_con_precomp, b_con_precomp, c_con_precomp, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s1 = s1 + tiempo_total
        if set_tmin_con_precomp == 0:
            tiempo_min_con_precomp = tiempo_min_con_precomp + tiempo_total
            set_tmin_con_precomp = 1
        elif tiempo_total < tiempo_min_con_precomp:
            tiempo_min_con_precomp = tiempo_total
        tiempo_inicial = time()
        EC_Tridiagonal2(A_1, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s2 = s2 + tiempo_total
        if set_tmin_sin_precomp == 0:
            tiempo_min_sin_precomp = tiempo_total
            set_tmin_sin_precomp = 1
        elif tiempo_total < tiempo_min_sin_precomp:
            tiempo_min_sin_precomp = tiempo_total
        m = m+1
    s1 = s1 + (t_end_precomp - t_init_precomp)
    s1_todas_las_matrices = s1_todas_las_matrices + s1/100
    s2_todas_las_matrices = s2_todas_las_matrices + s2/100
    n_matrices = n_matrices + 1
minimos_precomp.append(tiempo_min_con_precomp)
minimos_sin_precomp.append(tiempo_min_sin_precomp)
print("Tiempo promedio para Eliminacion con precómputo: ", s1_todas_las_matrices/10)
print("Tiempo promedio para Eliminación sin precómputo: ", s2_todas_las_matrices/10)
print("Tiempo mínimo para Eliminación con precómputo: ", tiempo_min_con_precomp)
print("Tiempo mínimo para Eliminación sin precómputo: ", tiempo_min_sin_precomp)


print("-------------")
print("Para n = 100")
# vamos a calcular el tiempo promedio para 10 matrices distintas
# para cada matriz, vamos a resolver el sistema (con y sin precómputo) unas 100 veces, y nos vamos a quedar con el tiempo promedio
# la matriz se va a mantener igual, lo que vamos a ir cambiando es el d para poder mostrar cómo se aprovecha el precómputo
n_matrices = 0
s1_todas_las_matrices = 0
s2_todas_las_matrices = 0
tiempo_min_sin_precomp = 0
tiempo_min_con_precomp = 0
set_tmin_sin_precomp = 0
set_tmin_con_precomp = 0
while n_matrices < 10:
    m = 0
    s1 = 0
    s2 = 0
    n = 100
    a_sin_precomp = np.insert(np.random.rand(n-1), 0, 0)
    b_sin_precomp = np.random.rand(n)
    c_sin_precomp = np.append(np.random.rand(n-1), 0)
    a_con_precomp = a_sin_precomp.copy()
    b_con_precomp = b_sin_precomp.copy()
    c_con_precomp = c_sin_precomp.copy()
    t_init_precomp = time()
    precomp_sist_tridiagonal(a_con_precomp, b_con_precomp, c_con_precomp)
    t_end_precomp = time()
    while m < 100:
        d = np.random.rand(n)
        A_1 = np.diag(a_sin_precomp[1:n].copy(), k=-1) + np.diag(b_sin_precomp.copy(), k=0) + np.diag(c_sin_precomp[0:(n-1)].copy(), k=1)
        tiempo_inicial = time()
        elim_gauss_sist_tridiagonal_precomp(a_con_precomp, b_con_precomp, c_con_precomp, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s1 = s1 + tiempo_total
        if set_tmin_con_precomp == 0:
            tiempo_min_con_precomp = tiempo_min_con_precomp + tiempo_total
            set_tmin_con_precomp = 1
        elif tiempo_total < tiempo_min_con_precomp:
            tiempo_min_con_precomp = tiempo_total
        tiempo_inicial = time()
        EC_Tridiagonal2(A_1, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s2 = s2 + tiempo_total
        if set_tmin_sin_precomp == 0:
            tiempo_min_sin_precomp = tiempo_total
            set_tmin_sin_precomp = 1
        elif tiempo_total < tiempo_min_sin_precomp:
            tiempo_min_sin_precomp = tiempo_total
        m = m+1
    s1 = s1 + (t_end_precomp - t_init_precomp)
    s1_todas_las_matrices = s1_todas_las_matrices + s1/100
    s2_todas_las_matrices = s2_todas_las_matrices + s2/100
    n_matrices = n_matrices + 1
minimos_precomp.append(tiempo_min_con_precomp)
minimos_sin_precomp.append(tiempo_min_sin_precomp)
print("Tiempo promedio para Eliminacion con precómputo: ", s1_todas_las_matrices/10)
print("Tiempo promedio para Eliminación sin precómputo: ", s2_todas_las_matrices/10)
print("Tiempo mínimo para Eliminación con precómputo: ", tiempo_min_con_precomp)
print("Tiempo mínimo para Eliminación sin precómputo: ", tiempo_min_sin_precomp)

print("-------------")
print("Para n = 500")
# vamos a calcular el tiempo promedio para 10 matrices distintas
# para cada matriz, vamos a resolver el sistema (con y sin precómputo) unas 100 veces, y nos vamos a quedar con el tiempo promedio
# la matriz se va a mantener igual, lo que vamos a ir cambiando es el d para poder mostrar cómo se aprovecha el precómputo
n_matrices = 0
s1_todas_las_matrices = 0
s2_todas_las_matrices = 0
tiempo_min_sin_precomp = 0
tiempo_min_con_precomp = 0
set_tmin_sin_precomp = 0
set_tmin_con_precomp = 0
while n_matrices < 10:
    m = 0
    s1 = 0
    s2 = 0
    n = 500
    a_sin_precomp = np.insert(np.random.rand(n-1), 0, 0)
    b_sin_precomp = np.random.rand(n)
    c_sin_precomp = np.append(np.random.rand(n-1), 0)
    a_con_precomp = a_sin_precomp.copy()
    b_con_precomp = b_sin_precomp.copy()
    c_con_precomp = c_sin_precomp.copy()
    t_init_precomp = time()
    precomp_sist_tridiagonal(a_con_precomp, b_con_precomp, c_con_precomp)
    t_end_precomp = time()
    while m < 100:
        d = np.random.rand(n)
        A_1 = np.diag(a_sin_precomp[1:n].copy(), k=-1) + np.diag(b_sin_precomp.copy(), k=0) + np.diag(c_sin_precomp[0:(n-1)].copy(), k=1)
        tiempo_inicial = time()
        elim_gauss_sist_tridiagonal_precomp(a_con_precomp, b_con_precomp, c_con_precomp, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s1 = s1 + tiempo_total
        if set_tmin_con_precomp == 0:
            tiempo_min_con_precomp = tiempo_min_con_precomp + tiempo_total
            set_tmin_con_precomp = 1
        elif tiempo_total < tiempo_min_con_precomp:
            tiempo_min_con_precomp = tiempo_total
        tiempo_inicial = time()
        EC_Tridiagonal2(A_1, d.copy())
        tiempo_final = time()
        tiempo_total = tiempo_final - tiempo_inicial
        s2 = s2 + tiempo_total
        if set_tmin_sin_precomp == 0:
            tiempo_min_sin_precomp = tiempo_total
            set_tmin_sin_precomp = 1
        elif tiempo_total < tiempo_min_sin_precomp:
            tiempo_min_sin_precomp = tiempo_total
        m = m+1
    s1 = s1 + (t_end_precomp - t_init_precomp)
    s1_todas_las_matrices = s1_todas_las_matrices + s1/100
    s2_todas_las_matrices = s2_todas_las_matrices + s2/100
    n_matrices = n_matrices + 1
minimos_precomp.append(tiempo_min_con_precomp)
minimos_sin_precomp.append(tiempo_min_sin_precomp)
print("Tiempo promedio para Eliminacion con precómputo: ", s1_todas_las_matrices/10)
print("Tiempo promedio para Eliminación sin precómputo: ", s2_todas_las_matrices/10)
print("Tiempo mínimo para Eliminación con precómputo: ", tiempo_min_con_precomp)
print("Tiempo mínimo para Eliminación sin precómputo: ", tiempo_min_sin_precomp)

ns = ['5', '10', '100', '500']
tiempos_precomp = [1.4802455902099609e-05, 2.4304151535034178e-05, 0.00020071601867675785, 0.0010811941623687746]
tiempos_sin_precomp = [3.131628036499023e-05, 6.056380271911621e-05, 0.000592694044113159, 0.0031750121116638176]

x = np.arange(len(ns))
width = 0.35
#plt.bar(ns, tiempos_precomp)
#plt.title("Tiempos de ejecución para algoritmos de resolución de \nsistemas tridiagonales con y sin precómputo")
ax = plt.subplot()
ax.bar(x - width/2, tiempos_precomp, width, label='tridiagonal con precómputo')
ax.bar(x + width/2, tiempos_sin_precomp, width, label='tridiagonal sin precómputo')
ax.set_title("Tiempos promedio de ejecución para algoritmos de resolución de \nsistemas tridiagonales con y sin precómputo")
ax.set_xticks(x)
ax.set_xticklabels(ns)
ax.set_xlabel("Tamaño de la matriz")
ax.set_ylabel("Tiempo de ejecución")
ax.legend()
plt.savefig("ej_4b.png")
plt.show()

ax = plt.subplot()
ax.bar(x - width/2, minimos_precomp, width, label="trdiagonal con precómputo", color="purple")
ax.bar(x + width/2, minimos_sin_precomp, width, label='tridiagonal sin precómputo', color="yellow")
ax.set_title("Tiempos mínimos de ejecución para algoritmos de resolución de \nsistemas tridiagonales con y sin precómputo")
ax.set_xticks(x)
ax.set_xticklabels(ns)
ax.set_xlabel("Tamaño de la matriz")
ax.set_ylabel("Tiempo de ejecución")
ax.legend()
plt.savefig("ej_4b_minimos.png")
plt.show()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
