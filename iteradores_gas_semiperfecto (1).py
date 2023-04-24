import numpy as np

T0 = 288.15 # Temperatura bloque a nivel del mar

#### Iteradores:
## Compresor
def itera_compresor(rend_c, pi_c, Tet, aire):
    # Variables del iterador
    sigue_iterando = True
    n_iter = 0
    tolerancia = 1e-4
    iter_max = 100
    Tst_ant = Tet
    Tst = Tet

    while sigue_iterando:
        # Iteraciones
        n_iter = n_iter + 1

        # Cálculo temperatura 
        Tcomp = (Tet + Tst) / 2 # Modelo con temperatura promedio
        Tst_Tet = (pi_c**((aire.gamma(Tcomp)-1)/aire.gamma(Tcomp)) - 1)/rend_c + 1
        Tst = Tst_Tet * Tet # Temp final

        # Residuo entre iteración actual y anterior:
        residuo = np.abs(Tst - Tst_ant) / T0

        print("{:} --  Tst={:6.1f}K;  T_c={:6.1f}K;  deltaT={:8.6f};  Converge: {:}".format(n_iter, Tst, Tcomp, residuo, residuo < tolerancia))

        # Evaluar condición de salida
        if residuo < tolerancia:
            # Salida por convergencia
            sigue_iterando = False

        elif n_iter > iter_max:
            # Salida por iteraciones máximas
            sigue_iterando = False

        else:
            # Sigue iterando
            Tst_ant = Tst
    return Tst

## Cámara de combustión
def itera_combustor(Q43, Gs, Tet, aire):
    sigue_iterando = True
    n_iter = 0
    iter_max = 100
    tolerancia = 1e-5
    Tst_ant = Tet

    while sigue_iterando:
        # Iteraciones
        n_iter = n_iter + 1
        
        # Cámara de combustión
        Tcomb = (Tst_ant + Tet)/2
        Tst = Q43/Gs/aire.cp(Tcomb) + Tet

        # Residuo entre iteración actual y anterior:
        residuo = np.abs(Tst - Tst_ant) / T0

        print("{:} --  Tst={:6.1f}K;  T_c={:6.1f}K;  deltaT={:8.6f};  Converge: {:}".format(n_iter, Tst, Tcomb, residuo, residuo < tolerancia))

        # Evaluar condición de salida
        if residuo < tolerancia:
            # Salida por convergencia
            sigue_iterando = False

        elif n_iter > iter_max:
            # Salida por iteraciones máximas
            sigue_iterando = False

        else:
            # Sigue iterando
            Tst_ant = Tst
    return Tst

## Turbina
def itera_turbina(Wt, Gs, Tet, aire):
    sigue_iterando = True
    n_iter = 0
    iter_max = 100
    tolerancia = 1e-5
    Tst_ant = Tet

    while sigue_iterando:
        # Iteraciones
        n_iter = n_iter + 1
        
        # Cámara de combustión
        Tmed = (Tst_ant + Tet)/2
        Tst = Wt/Gs/aire.cp(Tmed) + Tet

        # Residuo entre iteración actual y anterior:
        residuo = np.abs(Tst - Tst_ant) / T0

        print("{:} --  Tst={:6.1f}K;  T_c={:6.1f}K;  deltaT={:8.6f};  Converge: {:}".format(n_iter, Tst, Tmed, residuo, residuo < tolerancia))

        # Evaluar condición de salida
        if residuo < tolerancia:
            # Salida por convergencia
            sigue_iterando = False

        elif n_iter > iter_max:
            # Salida por iteraciones máximas
            sigue_iterando = False

        else:
            # Sigue iterando
            Tst_ant = Tst
    return Tst