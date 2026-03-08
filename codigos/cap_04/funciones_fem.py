"""
Este archivo coontiene funciones adaptadas del método de los elementos finitos 
para el cálculo de los esfuerzos y las deformaciones al interior de vigas.

Adaptado de: https://github.com/diegoandresalvarez/elementosfinitos/blob/master/codigo/2D/ejemplo_Q6/python/funciones.py
Apoyo en IA: Claude CODE.

Por: Michael Heredia Pérez
email: mherediap@unal.edu.co
Fecha: 2026-03-26

Universidad Nacional de Colombia
"""

import numpy as np
from scipy.linalg import solve
import matplotlib.tri as tri


# constantes que ayudaran en la lectura del codigo.
X, Y = 0, 1
NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8 = range(8)

# variables globales que se heredarán del programa principal.
xnod = None
LaG  = None


def compartir_variables(xnod_, LaG_):
    """Importa variables globales del programa principal a este módulo.

    Args:
        xnod_ (_type_): _description_
        LaG_ (_type_): _description_
    """
    
    '''
    Importa variables globales del programa principal a este módulo
    '''
    global xnod, LaG
    xnod, LaG = xnod_, LaG_
        
        
def t2ft_R4(xnod, lado, carga, espesor):
    """
    Función que convierte las fuerzas superficiales aplicadas a un elemento
    finito rectangular de 4 nodos a sus correspondientes cargas nodales 
    equivalentes ft
    
    Args:
        xnod:  coordenadas nodales del elemento finito de 4 nodos
          xnod = [ x1e y1e
                   x2e y2e
                   x3e y3e
                   x4e y4e ]

        lado:  arista en la que se aplica la carga, puede tomar los siguientes
               valores: 12, 23, 34, 41

        carga: fuerza distribuida en los nodos
        
               [ t1x t1y t2x t2y ]; % si carga se aplica sobre lado 12
               [ t2x t2y t3x t3y ]; % si carga se aplica sobre lado 23
               [ t3x t3y t4x t4y ]; % si carga se aplica sobre lado 34
               [ t4x t4y t1x t1y ]; % si carga se aplica sobre lado 41
    
        espesor: espesor del elemento
    """
    
    # se definen los indices de los lados
    if   lado == 12: idx = np.array([ 1, 2 ]) - 1
    elif lado == 23: idx = np.array([ 2, 3 ]) - 1
    elif lado == 34: idx = np.array([ 3, 4 ]) - 1
    elif lado == 41: idx = np.array([ 4, 1 ]) - 1
    else: 
        raise Exception('Únicamente se permiten los lados 12, 23, 34 o 41')

    nno = xnod.shape[0]
    if nno != 4:
        raise Exception('Solo para elementos rectangulares de 4 nodos')

    # parámetros para mejorar la lectura del código
    X, Y = 0, 1
   
    # se define el número de puntos de la cuadratura y se obtienen los puntos
    # de evaluación y los pesos de la cuadratura
    n_gl       = 5
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_gl)
    
    # se definen las funciones de forma unidimensionales y sus derivadas
    NN      = lambda xi: np.array([ (1-xi)/2, (1+xi)/2 ])
    dNN_dxi = lambda xi: np.array([     -1/2,      1/2 ])
       
    # se calcula el vector de fuerzas distribuidas en los nodos
    te = np.zeros(2*nno)
    te[np.c_[2*idx, 2*idx + 1].ravel()] = carga
    
    # cálculo de la integral:
    suma   = np.zeros((2*nno, 2*nno))
    N      = np.zeros(nno)
    dN_dxi = np.zeros(nno)

    for p in range(n_gl):

        # se evalúan las funciones de forma
        N[idx] = NN(x_gl[p])

        matN = np.hstack([ np.array([[N[i],    0],
                                     [    0, N[i]]]) for i in range(nno)])

        # se calcula el jacobiano
        dN_dxi[idx] = dNN_dxi(x_gl[p])
        dx_dxi      = np.dot(dN_dxi, xnod[:,X])
        dy_dxi      = np.dot(dN_dxi, xnod[:,Y])
        ds_dxi      = np.hypot(dx_dxi, dy_dxi)

        # y se calcula la sumatoria
        suma += matN.T @ matN * ds_dxi*w_gl[p]

    # Se retorna el vector de fuerzas nodales equivalentes
    return espesor * (suma @ te)


def generar_malla_Q4(L_total: float, h: float, nx: int, ny: int) -> dict:
    """
    Genera la malla estructural de elementos finitos cuadriláteros Q4
    para una viga rectangular 2D.

    Parámetros
    ----------
    L_total : float
        Longitud total de la viga [m].
    h : float
        Altura (canto) de la viga [m].
    nx : int
        Número de divisiones en la dirección x (horizontal).
    ny : int
        Número de divisiones en la dirección y (vertical).

    Retorna
    -------
    dict con las claves:
        'xnod' : np.ndarray, shape (nno, 2)
            Coordenadas (x, y) de cada nodo.
        'LaG'  : np.ndarray, shape (nef, 4)
            Tabla de conectividad. Cada fila contiene los índices
            [n1, n2, n3, n4] de los nodos del elemento, en sentido
            antihorario desde la esquina inferior-izquierda.
        'nno'  : int
            Número total de nodos.
        'nef'  : int
            Número total de elementos.
    """
    # --- Coordenadas nodales ---
    x_vec = np.linspace(0, L_total, nx + 1)
    y_vec = np.linspace(0, h,       ny + 1)
    X, Y  = np.meshgrid(x_vec, y_vec)
    xnod  = np.column_stack((X.flatten(), Y.flatten()))
    nno   = xnod.shape[0]

    # --- Tabla de conectividad (LaG) ---
    # Numeración local del elemento Q4:
    #
    #   n4 ------- n3
    #    |          |
    #    |          |
    #   n1 ------- n2
    #
    LaG = np.array([
        [i     +  j   *(nx+1),
         i + 1 +  j   *(nx+1),
         i + 1 + (j+1)*(nx+1),
         i     + (j+1)*(nx+1)]
        for j in range(ny)
        for i in range(nx)
    ])
    nef = LaG.shape[0]

    print(f"Malla generada: {nno} nodos | {nef} elementos Q4  "
          f"[nx={nx}, ny={ny}, L={L_total:.2f} m, h={h:.2f} m]")

    return {"xnod": xnod, "LaG": LaG, "nno": nno, "nef": nef}


def ensamblar_rigidez_Q4_incompatible(
    xnod: np.ndarray,
    LaG:  np.ndarray,
    nno:  int,
    nef:  int,
    E:    float,
    nu:   float,
    t:    float
) -> dict:
    """
    Ensambla la matriz de rigidez global K para elementos Q4 con modos
    incompatibles (Wilson, Taylor et al.), usando condensación estática
    para eliminar los grados de libertad internos.

    Parámetros
    ----------
    xnod : np.ndarray, shape (nno, 2)
        Coordenadas (x, y) de cada nodo.
    LaG  : np.ndarray, shape (nef, 4)
        Tabla de conectividad.
    nno  : int
        Número total de nodos.
    nef  : int
        Número total de elementos.
    E    : float
        Módulo de elasticidad [Pa].
    nu   : float
        Coeficiente de Poisson [-].
    t    : float
        Espesor del elemento (esfuerzo plano) [m].

    Retorna
    -------
    dict con las claves:
        'K'            : np.ndarray, shape (2*nno, 2*nno)
            Matriz de rigidez global ensamblada (con condensación estática).
        'f'            : np.ndarray, shape (2*nno,)
            Vector de fuerzas externas inicializado en cero.
        'B_storage'    : list[np.ndarray]
            Matrices B de cada elemento en los 4 puntos de Gauss,
            shape por elemento: (2, 2, 3, 12).
        'inv_Kee_list' : list[np.ndarray]
            Inversas de las submatrices Kee (4×4) de cada elemento.
        'Ker_list'     : list[np.ndarray]
            Submatrices Ker (4×8) de acoplamiento de cada elemento.
    """
    # --- Matriz constitutiva (esfuerzo plano) ---
    De = (E / (1 - nu**2)) * np.array([
        [1,      nu,          0],
        [nu,     1,           0],
        [0,      0,  (1-nu)/2  ]
    ])

    # --- Funciones de forma para modos incompatibles (6 nodos: 4 reales + 2 burbuja) ---
    #
    #   Nodos reales (i = 0..3): funciones bilineales estándar Q4
    #   Modos incompatibles (i = 4, 5): modos burbuja cuadráticos
    #
    #   dN/dxi  → derivadas respecto a xi
    #   dN/deta → derivadas respecto a eta
    dN_dxi  = lambda xi, eta: np.array([
         eta/4 - 1/4,   # N1,xi
         1/4 - eta/4,   # N2,xi
         eta/4 + 1/4,   # N3,xi
        -eta/4 - 1/4,   # N4,xi
        -2*xi,          # M1,xi  (modo incompatible)
         0              # M2,xi  (modo incompatible)
    ])
    dN_deta = lambda xi, eta: np.array([
         xi/4 - 1/4,    # N1,eta
        -xi/4 - 1/4,    # N2,eta
         xi/4 + 1/4,    # N3,eta
         1/4 - xi/4,    # N4,eta
         0,             # M1,eta (modo incompatible)
        -2*eta          # M2,eta (modo incompatible)
    ])

    # --- Puntos y pesos de cuadratura de Gauss 2×2 ---
    x_gl, w_gl = np.polynomial.legendre.leggauss(2)

    # --- Inicialización ---
    K          = np.zeros((2*nno, 2*nno))
    f          = np.zeros(2*nno)
    B_storage  = []
    inv_Kee_list = []
    Ker_list     = []

    # --- Slices para condensación estática ---
    # 8 GDL reales (nodos Q4) | 4 GDL incompatibles (modos burbuja)
    r_ = slice(0, 8)
    e_ = slice(8, 12)

    # --- Bucle de ensamblaje ---
    for e in range(nef):
        xe, ye = xnod[LaG[e], 0], xnod[LaG[e], 1]

        # Matriz de rigidez local expandida (8 reales + 4 incompatibles = 12 GDL)
        Ke12  = np.zeros((12, 12))
        # Almacén de matrices B en los 4 puntos de Gauss del elemento
        B_el  = np.empty((2, 2, 3, 12))

        # --- Cuadratura de Gauss 2×2 ---
        for p in range(2):
            for q in range(2):
                xi, eta  = x_gl[p], x_gl[q]
                dxi      = dN_dxi(xi, eta)
                deta     = dN_deta(xi, eta)

                # Jacobiano de la transformación isoparamétrica
                Je   = np.array([
                    [np.dot(dxi[:4],  xe), np.dot(dxi[:4],  ye)],
                    [np.dot(deta[:4], xe), np.dot(deta[:4], ye)]
                ])
                detJ = np.linalg.det(Je)

                # Construcción de la matriz B (deformación-desplazamiento, 3×12)
                Bpq = np.zeros((3, 12))
                for i in range(6):
                    dN_dx = ( Je[1,1]*dxi[i]  - Je[0,1]*deta[i]) / detJ
                    dN_dy = (-Je[1,0]*dxi[i]  + Je[0,0]*deta[i]) / detJ
                    Bpq[:, [2*i, 2*i+1]] = [
                        [dN_dx,     0  ],
                        [0,      dN_dy ],
                        [dN_dy,  dN_dx ]
                    ]

                B_el[p, q] = Bpq
                Ke12 += Bpq.T @ De @ Bpq * detJ * t * w_gl[p] * w_gl[q]

        B_storage.append(B_el)

        # --- Condensación estática de los modos incompatibles ---
        #
        #   K_condensada = Krr - Kre · Kee⁻¹ · Ker
        #
        invKee = np.linalg.inv(Ke12[e_, e_])
        Ker    = Ke12[e_, r_]
        inv_Kee_list.append(invKee)
        Ker_list.append(Ker)

        Ke_cond = Ke12[r_, r_] - Ke12[r_, e_] @ invKee @ Ker

        # --- Índices globales del elemento (GDL x e y de los 4 nodos) ---
        nodes     = LaG[e]
        idx_total = np.array([
            2*nodes[0],   2*nodes[0]+1,
            2*nodes[1],   2*nodes[1]+1,
            2*nodes[2],   2*nodes[2]+1,
            2*nodes[3],   2*nodes[3]+1
        ])

        K[np.ix_(idx_total, idx_total)] += Ke_cond

    print(f"Ensamblaje completado: K global {K.shape} | "
          f"{nef} elementos procesados.")

    return {
        "K":            K,
        "f":            f,
        "B_storage":    B_storage,
        "inv_Kee_list": inv_Kee_list,
        "Ker_list":     Ker_list,
        "De":           De           
    }
    

def aplicar_condiciones_frontera(
    xnod:      np.ndarray,
    LaG:       np.ndarray,
    nno:       int,
    nef:       int,
    f:         np.ndarray,
    w_carga:   float,
    t:         float,
    h:         float,
    L:         float,
    x_inicio_carga: float = 0.0
) -> dict:
    """
    Aplica cargas nodales equivalentes y define las condiciones de
    frontera (apoyos) del modelo de viga 2D con elementos Q4.

    Esquema estructural
    -------------------
                    w_carga [N/m]
            ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
    A ▽─────────────────────────────── B ▽
      |← x_inicio_carga →|←  zona cargada  →|
      x=0                x=L_total

    Apoyos
    ------
        A (x=0,  y=0): restringido en X e Y  →  GDL [2·idA,   2·idA+1]
        B (x=L,  y=0): restringido solo en Y →  GDL [2·idB+1]

    Parámetros
    ----------
    xnod           : np.ndarray, shape (nno, 2)
        Coordenadas nodales.
    LaG            : np.ndarray, shape (nef, 4)
        Tabla de conectividad.
    nno            : int
        Número total de nodos.
    nef            : int
        Número total de elementos.
    f              : np.ndarray, shape (2*nno,)
        Vector de fuerzas externas (se actualiza in-place y se retorna).
    w_carga        : float
        Intensidad de la carga distribuida [N/m], positiva hacia arriba.
    t              : float
        Espesor del elemento [m].
    h              : float
        Altura de la viga [m] — identifica la cara superior (y = h).
    L              : float
        Posición x del apoyo B [m].
    x_inicio_carga : float, opcional
        Coordenada x desde la cual comienza la carga distribuida [m].
        Por defecto 0.0 (carga en toda la longitud).

    Retorna
    -------
    dict con las claves:
        'f'            : np.ndarray, shape (2*nno,)
            Vector de fuerzas actualizado con las cargas nodales equivalentes.
        'restringidos' : np.ndarray
            Índices de GDL con desplazamiento impuesto (= 0).
        'libres'       : np.ndarray
            Índices de GDL sin restricción cinemática.
        'id_A'         : int
            Índice del nodo de apoyo A.
        'id_B'         : int
            Índice del nodo de apoyo B.

    Dependencias externas
    ---------------------
    t2ft_R4(xnod_e, cara, cargas, t) → np.ndarray, shape (8,)
        Convierte carga distribuida en la cara de un Q4 a fuerzas
        nodales equivalentes por integración numérica.
        (debe estar definida en el entorno antes de llamar esta función)
    """

    # ------------------------------------------------------------------ #
    # 1. CARGAS NODALES EQUIVALENTES (cara superior, zona de carga)       #
    # ------------------------------------------------------------------ #
    n_elementos_cargados = 0

    for e in range(nef):
        nodes = LaG[e]

        # Nodos de la cara superior del elemento: n3 y n4 (índices 2 y 3 en LaG)
        nodos_top    = nodes[[2, 3]]
        y_top        = xnod[nodos_top, 1]
        x_top        = xnod[nodos_top, 0]

        # Condición: cara superior (y ≈ h) y dentro de la zona cargada (x ≥ x_inicio_carga)
        en_cara_top  = np.all(np.isclose(y_top, h))
        en_zona_carga = np.all(x_top >= x_inicio_carga)

        if en_cara_top and en_zona_carga:
            idx_e = np.array([
                2*nodes[0],   2*nodes[0]+1,
                2*nodes[1],   2*nodes[1]+1,
                2*nodes[2],   2*nodes[2]+1,
                2*nodes[3],   2*nodes[3]+1
            ])
            # Cara 34: nodos n3-n4 (cara superior en numeración local Q4)
            # Carga en dirección -y (gravedad): [0, -w, 0, -w]
            f[idx_e] += t2ft_R4(
                xnod[nodes],
                lado    = 34,
                carga   = [0, -w_carga, 0, -w_carga],
                espesor = t
            )
            n_elementos_cargados += 1

    print(f"Carga aplicada en {n_elementos_cargados} elementos "
          f"(cara superior, x >= {x_inicio_carga:.2f} m).")

    # ------------------------------------------------------------------ #
    # 2. IDENTIFICACIÓN DE NODOS DE APOYO                                 #
    # ------------------------------------------------------------------ #
    # Apoyo A: esquina inferior-izquierda (x=0, y=0)
    id_A = np.argmin(np.hypot(xnod[:, 0] - 0, xnod[:, 1] - 0))

    # Apoyo B: nodo inferior en x=L (x=L, y=0)
    id_B = np.argmin(np.hypot(xnod[:, 0] - L, xnod[:, 1] - 0))

    print(f"Apoyo A (articulado) → nodo {id_A}  {xnod[id_A]}  "
          f"| GDL restringidos: {[2*id_A, 2*id_A+1]}")
    print(f"Apoyo B (rodillo)    → nodo {id_B}  {xnod[id_B]}  "
          f"| GDL restringido:  {[2*id_B+1]}")

    # ------------------------------------------------------------------ #
    # 3. PARTICIÓN GDL LIBRES / RESTRINGIDOS                              #
    # ------------------------------------------------------------------ #
    # A: articulado → restringe Ax (x) y Ay (y)
    # B: rodillo    → restringe solo By (y)
    restringidos = np.array([2*id_A, 2*id_A+1, 2*id_B+1])
    libres       = np.setdiff1d(np.arange(2*nno), restringidos)

    print(f"GDL totales: {2*nno}  |  libres: {len(libres)}  |  "
          f"restringidos: {len(restringidos)}")

    return {
        "f":            f,
        "restringidos": restringidos,
        "libres":       libres,
        "id_A":         id_A,
        "id_B":         id_B
    }
    
    
def resolver_sistema(
    K:            np.ndarray,
    f:            np.ndarray,
    nno:          int,
    nef:          int,
    libres:       np.ndarray,
    LaG:          np.ndarray,
    inv_Kee_list: list,
    Ker_list:     list,
    B_storage:    list,
    De:           np.ndarray,
    id_A:         int,
    id_B:         int
) -> dict:
    """
    Resuelve el sistema de ecuaciones FEM, extrapola esfuerzos a nodos
    y calcula las reacciones en los apoyos.

    Parámetros
    ----------
    K            : np.ndarray, shape (2*nno, 2*nno)
        Matriz de rigidez global ensamblada.
    f            : np.ndarray, shape (2*nno,)
        Vector de fuerzas externas.
    nno          : int
        Número total de nodos.
    nef          : int
        Número total de elementos.
    libres       : np.ndarray
        Índices de GDL libres.
    LaG          : np.ndarray, shape (nef, 4)
        Tabla de conectividad.
    inv_Kee_list : list[np.ndarray]
        Inversas de las submatrices Kee por elemento (condensación estática).
    Ker_list     : list[np.ndarray]
        Submatrices Ker de acoplamiento por elemento.
    B_storage    : list[np.ndarray]
        Matrices B en los 4 puntos de Gauss por elemento, shape (2,2,3,12).
    De           : np.ndarray, shape (3, 3)
        Matriz constitutiva (esfuerzo plano).
    id_A         : int
        Índice del nodo de apoyo A (articulado).
    id_B         : int
        Índice del nodo de apoyo B (rodillo).

    Retorna
    -------
    dict con las claves:
        'a'    : np.ndarray, shape (2*nno,)
            Vector de desplazamientos nodales globales [m].
        'sx'   : np.ndarray, shape (nno,)
            Esfuerzo normal en x extrapolado a nodos [kPa].
        'sy'   : np.ndarray, shape (nno,)
            Esfuerzo normal en y extrapolado a nodos [kPa].
        'txy'  : np.ndarray, shape (nno,)
            Esfuerzo cortante extrapolado a nodos [kPa].
        'R_Ax' : float
            Reacción horizontal en apoyo A [kN].
        'R_Ay' : float
            Reacción vertical en apoyo A [kN].
        'R_By' : float
            Reacción vertical en apoyo B [kN].
    """

    # ------------------------------------------------------------------ #
    # 1. SOLUCIÓN DEL SISTEMA LINEAL (solo GDL libres)                    #
    # ------------------------------------------------------------------ #
    a          = np.zeros(2*nno)
    a[libres]  = solve(K[np.ix_(libres, libres)], f[libres])

    # ------------------------------------------------------------------ #
    # 2. EXTRAPOLACIÓN DE ESFUERZOS DE GAUSS A NODOS                      #
    # ------------------------------------------------------------------ #
    # Matriz de extrapolación estándar para Q4 (puntos de Gauss → nodos)
    #
    #   Orden de puntos de Gauss:  (p=0,q=0), (p=0,q=1),
    #                              (p=1,q=0), (p=1,q=1)
    #   Orden de nodos locales:    n1, n2, n3, n4
    #
    extrap = np.array([
        [ 1.866, -0.500, -0.500,  0.134],
        [-0.500,  0.134,  1.866, -0.500],
        [ 0.134, -0.500, -0.500,  1.866],
        [-0.500,  1.866,  0.134, -0.500]
    ])

    sx      = np.zeros(nno)
    sy      = np.zeros(nno)
    txy     = np.zeros(nno)
    num_adj = np.zeros(nno)

    for e in range(nef):
        nodes = LaG[e]
        idx_r = np.array([
            2*nodes[0],   2*nodes[0]+1,
            2*nodes[1],   2*nodes[1]+1,
            2*nodes[2],   2*nodes[2]+1,
            2*nodes[3],   2*nodes[3]+1
        ])

        # Desplazamientos reales + recuperación de GDL incompatibles
        # ae = [a_reales (8) | a_incompatibles (4)]
        ae = np.r_[a[idx_r],
                   -inv_Kee_list[e] @ Ker_list[e] @ a[idx_r]]

        # Esfuerzos en los 4 puntos de Gauss: σ = De · B · ae
        eg = np.array([
            De @ B_storage[e][p, q] @ ae
            for p in range(2)
            for q in range(2)
        ])                                      # shape (4, 3)

        # Extrapolación a los 4 nodos del elemento y acumulación
        sx[nodes]      += extrap @ eg[:, 0]
        sy[nodes]      += extrap @ eg[:, 1]
        txy[nodes]     += extrap @ eg[:, 2]
        num_adj[nodes] += 1

    # Promediado nodal y conversión Pa → kPa
    sx  /= (num_adj * 1000)
    sy  /= (num_adj * 1000)
    txy /= (num_adj * 1000)

    # ------------------------------------------------------------------ #
    # 3. REACCIONES EN APOYOS                                             #
    # ------------------------------------------------------------------ #
    # R = K·a − f  (vector de reacciones globales)
    reacciones = K @ a - f

    R_Ax = reacciones[2*id_A]     / 1000   # Horizontal en A  [kN]
    R_Ay = reacciones[2*id_A + 1] / 1000   # Vertical en A    [kN]
    R_By = reacciones[2*id_B + 1] / 1000   # Vertical en B    [kN]

    print(f"Desplazamiento máximo : {np.max(np.abs(a)):.6e} m")
    print(f"Reacciones → R_Ax = {R_Ax:.3f} kN | "
          f"R_Ay = {R_Ay:.3f} kN | R_By = {R_By:.3f} kN")

    return {
        "a":    a,
        "sx":   sx,
        "sy":   sy,
        "txy":  txy,
        "R_Ax": R_Ax,
        "R_Ay": R_Ay,
        "R_By": R_By
    }
    

def calcular_esf_def_V_M(
    xnod:  np.ndarray,
    sx:    np.ndarray,
    sy:    np.ndarray,
    txy:   np.ndarray,
    nno:   int,
    nx:    int,
    ny:    int,          # ← nuevo parámetro
    h:     float,
    t:     float,
    x_vec: np.ndarray,
    E:     float,       # ← nuevo
    nu:    float,       # ← nuevo
) -> dict:
    """
    Calcula los esfuerzos principales (círculo de Mohr), la triangulación
    para graficación y los diagramas de fuerza cortante (V) y momento
    flector (M) a partir de los esfuerzos nodales extrapolados.

    Parámetros
    ----------
    xnod  : np.ndarray, shape (nno, 2)
        Coordenadas nodales [m].
    sx    : np.ndarray, shape (nno,)
        Esfuerzo normal en x por nodo [kPa].
    sy    : np.ndarray, shape (nno,)
        Esfuerzo normal en y por nodo [kPa].
    txy   : np.ndarray, shape (nno,)
        Esfuerzo cortante por nodo [kPa].
    nno   : int
        Número total de nodos.
    nx    : int
        Número de divisiones en x (secciones transversales = nx + 1).
    h     : float
        Altura de la viga [m].
    t     : float
        Espesor del elemento [m].
    x_vec : np.ndarray, shape (nx + 1,)
        Coordenadas x de cada sección transversal [m].

    Retorna
    -------
    dict con las claves:
        's1'     : np.ndarray, shape (nno,)
            Esfuerzo principal mayor en cada nodo [kPa].
        's2'     : np.ndarray, shape (nno,)
            Esfuerzo principal menor en cada nodo [kPa].
        'tmax'   : np.ndarray, shape (nno,)
            Cortante máximo (radio de Mohr) en cada nodo [kPa].
        'ang_p1' : np.ndarray, shape (nno,)
            Ángulo de la dirección principal mayor [rad].
        'ang_p2' : np.ndarray, shape (nno,)
            Ángulo de la dirección principal menor [rad].
        'triang' : matplotlib.tri.Triangulation
            Objeto de triangulación para graficación sobre la malla.
        'V_diag' : np.ndarray, shape (nx + 1,)
            Diagrama de fuerza cortante en cada sección [kN].
        'M_diag' : np.ndarray, shape (nx + 1,)
            Diagrama de momento flector en cada sección [kN·m].
        'x_vec'  : np.ndarray, shape (nx + 1,)
            Coordenadas x de cada sección transversal [m].
    """

    # ------------------------------------------------------------------ #
    # 1. CÍRCULO DE MOHR: ESFUERZOS PRINCIPALES Y ÁNGULOS                 #
    # ------------------------------------------------------------------ #
    R_mohr = np.sqrt(((sx - sy) / 2)**2 + txy**2)
    prom   = (sx + sy) / 2

    s1     = prom + R_mohr          # Esfuerzo principal mayor  [kPa]
    s2     = prom - R_mohr          # Esfuerzo principal menor  [kPa]
    tmax   = R_mohr                 # Cortante máximo           [kPa]

    # Ángulos de las direcciones principales respecto al eje x
    ang_p1 = 0.5 * np.arctan2(2*txy, sx - sy)
    ang_p2 = ang_p1 + np.pi / 2

    # ------------------------------------------------------------------ #
    # 2. TRIANGULACIÓN PARA GRAFICACIÓN                                   #
    # ------------------------------------------------------------------ #
    triang = tri.Triangulation(xnod[:, 0], xnod[:, 1])

    # ------------------------------------------------------------------ #
    # 3. DIAGRAMAS DE FUERZA CORTANTE Y MOMENTO FLECTOR                   #
    # ------------------------------------------------------------------ #
    # Integración por sección transversal usando la regla del trapecio:
    #
    #   V(x) = t · ∫ τ_xy(x, y) dy       [kN]
    #   M(x) = t · ∫ σ_x(x, y)·(y − ȳ) dy   [kN·m]
    #
    # donde ȳ = h/2 es la fibra neutra de la sección rectangular.

    y_centro = h / 2
    V_diag   = np.zeros(nx + 1)
    M_diag   = np.zeros(nx + 1)

    for i in range(nx + 1):
        # Nodos de la sección transversal en x = x_vec[i], ordenados por y
        idx_col   = np.where(np.abs(xnod[:, 0] - x_vec[i]) < 1e-6)[0]
        sort_idx  = np.argsort(xnod[idx_col, 1])
        nodos_col = idx_col[sort_idx]
        y_loc     = xnod[nodos_col, 1]

        # V = t · ∫ τ_xy dy  (ya en kPa → resultado en kN)
        V_diag[i] = np.trapz(txy[nodos_col] * t, y_loc)

        # M = t · ∫ σ_x · (y − ȳ) dy  (resultado en kN·m)
        brazo     = y_loc - y_centro
        M_diag[i] = np.trapz(sx[nodos_col] * brazo * t, y_loc)

    print(f"V_max = {np.max(np.abs(V_diag)):.3f} kN  |  "
          f"M_max = {np.max(np.abs(M_diag)):.3f} kN·m")

    # ------------------------------------------------------------------ #
    # 4. DEFORMACIONES (Ley de Hooke, tensión plana)                      #
    # ------------------------------------------------------------------ #
    G   = E / (2*(1 + nu))          # Módulo de cortante [Pa]

    # sx, sy, txy están en kPa → convertir a Pa para aplicar Hooke
    sx_Pa  = sx  * 1000
    sy_Pa  = sy  * 1000
    txy_Pa = txy * 1000

    ex  = (1/E) * (sx_Pa - nu*sy_Pa)   # Deformación normal en x  [-]
    ey  = (1/E) * (sy_Pa - nu*sx_Pa)   # Deformación normal en y  [-]
    ez  = -(nu/E) * (sx_Pa + sy_Pa)    # Deformación normal en z  [-]
    gxy = txy_Pa / G                   # Deformación angular xy   [-]
    
    # ------------------------------------------------------------------ #
    # 5. RESHAPING A 2D PARA GRAFICACIÓN CON plot_esf_def                 #
    # ------------------------------------------------------------------ #
    # plot_esf_def(x, y, variable, ...) espera arrays 2D tipo meshgrid,
    # shape (ny+1, nx+1), coherente con la malla estructurada Q4.
    shape2D = (ny + 1, nx + 1)

    X   = xnod[:, 0].reshape(shape2D)
    Y   = xnod[:, 1].reshape(shape2D)
    SX  = sx.reshape(shape2D)
    SY  = sy.reshape(shape2D)
    TXY = txy.reshape(shape2D)
    S1  = s1.reshape(shape2D)
    S2  = s2.reshape(shape2D)
    TMAX   = tmax.reshape(shape2D)
    ANG_P1 = ang_p1.reshape(shape2D)
    ANG_P2 = ang_p2.reshape(shape2D)
    EX     = ex.reshape(shape2D)
    EY     = ey.reshape(shape2D)
    EZ     = ez.reshape(shape2D)
    GXY    = gxy.reshape(shape2D)
    
    return {
        # 1D
        "s1": s1, "s2": s2, "tmax": tmax,
        "ang_p1": ang_p1, "ang_p2": ang_p2,
        "V_diag": V_diag, "M_diag": M_diag, "x_vec": x_vec,
        "ex": ex, "ey": ey, "ez": ez, "gxy": gxy,
        # 2D
        "X": X, "Y": Y,
        "SX": SX, "SY": SY, "TXY": TXY,
        "S1": S1, "S2": S2, "TMAX": TMAX,
        "ANG_P1": ANG_P1, "ANG_P2": ANG_P2,
        "EX": EX, "EY": EY, "EZ": EZ, "GXY": GXY
    }
    
# Fin :)