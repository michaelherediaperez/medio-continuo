"""
Interpretación de gráficos de colores.

por: Michael Heredia Pérez
fecha: octubre 21/2022
email: mherediap@unal.edu.co

Traducción del código de:
- Diego Andrés Álvarez
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable     # Sólo para el tamaño de la colorbar
import os

def plot_esf_def(variable, titulo, nombre, angulo = None):
    """Función para graficar los esfuerzos y las deformaciones. Usar:
    
    Variables:
        var    : es la variable que se quiere graficar.
        titulo : título del gráfico.
        angulo : (opcional) inlcinaciones para los esf. ppls s1, s2 y tmax.
    """
    global x, y

    # Creo el lienzo.
    fig, ax = plt.subplots()
    
    # Coloreo el gráfico según la variable.
    im = ax.contourf(x, y, variable, 50, cmap = 'bwr', alpha=0.8)
    # Grafico las curvas de nivel.
    ax.contour(x, y, variable, 20, colors="black", linewidths=0.5, 
               linestyles="solid")
    # Creo un axes adicional a la derecha en donde va la escala de colores.
    divider = make_axes_locatable(ax)
    # Creo la escala de colores.
    cax = divider.append_axes("right", size="5%", pad=0.05)
    # Configuro la escala para mostrar 5 intervales/4 etiquetas.
    fig.colorbar(im, cax=cax).ax.locator_params(nbins=5)
        
    # Para los esfuerzos principales se grafican las líneas que indiquen las
    # direcciones de los esfuerzos en cada nodo de la malla
    if angulo is not None:
       if type(angulo) is np.ndarray: 
           angulo = [ angulo ]
       for ang in angulo:
           ax.quiver(x, y, variable*np.cos(ang), variable*np.sin(ang), 
                headwidth=0, headlength=0, headaxislength=0, pivot='middle', 
                linewidths=2
                )
    
    # Se especifican los ejes y el título, y se colocan los ejes iguales.
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title(titulo, fontsize=15)
    ax.set_aspect('equal')
    ax.autoscale(tight=True)    
    plt.tight_layout()

    #plt.show()

    # Almaceno el gráfico
    plt.savefig(f"graficos_colores/{nombre}.png")

# Variables y constantes.
# ------------------------------------------------------------------------------

# Propiedades geométricas de la viga
c = 0.50    # m, altura = 2c
L = 3.00    # m, luz    = 2L

# Se calcula la inercia I = bh^3/12, con b=1 y h=2c
I = (2*c**3)/3

# Carga aplicada
q = -10.0 # kN/m

# Propiedades del material
E  = 21e6         # kPa = 21 GPa, módulo de Young.
nu = 0.23         # ad          , coeficiente de Poisson.
G = E/(2*(1+nu))  # kPa         , módulo de cortante.

# Cálculos.
# ------------------------------------------------------------------------------

# Se crea la grilla de puntos donde se harán los cálculos.
nnds_x = 50
nnds_y = 20
x, y = np.meshgrid( np.linspace(-L, L, nnds_x), np.linspace(-c, c, nnds_y) )

# Se definen los esfuerzos (en tensión plana sz = txz = tyz = 0), eq (4.46)
sx  = -(q/(2*I))*(x**2*y - 2*y**3/3 + 2*c**2*y/5 - L**2*y)
sy  = -(q/(2*I))*(y**3/3 - c**2*y - 2*c**3/3)
txy = -(q/(2*I))*(c**2 - y**2)*x

# Se calculan las deformaiones, eq (4.36)
ex  = (1/E)*(sx - nu*sy)
ey  = (1/E)*(sy - nu*sx)
ez  = -(nu/E)*(sx + sy)
gxy = txy/G

# Se calculan los esfuerzos principales, los esfuerzos cortantes máximos y sus 
# ángulos.
tmax = np.sqrt( ((sx - sy)/2)**2 + txy**2 )
s1 = (sx + sy)/2 + tmax
s2 = (sx + sy)/2 - tmax
t1 = np.arctan2(2*txy, sx-sy)/2     
t2 = t1 + np.pi/2                   # t1 + 90°

# -----------------------------------------------------------------------------
# Carpeta para resultados.

# Check whether the specified path exists or not and if not create it.
path = "/graficos_colores"

if not os.path.exists(path):
    os.mkdir(path) 
    
print(f"\nGráficos y excel con resultados se guardarán en: {path}\n")

# Graficación.
# ------------------------------------------------------------------------------

# esfuerzos
plot_esf_def(sx,  r"$\sigma_x$ [Pa]", "sigma_x")
plot_esf_def(sy,  r"$\sigma_y$ [Pa]", "sigma_y")
plot_esf_def(txy, r"$\tau_{xy}$ [Pa]", "tau_xy")

# Deformaciones
plot_esf_def(ex,  r"$\epsilon_x$",        "epsilon_x")
plot_esf_def(ey,  r"$\epsilon_y$",        "epsilon_y")
plot_esf_def(ez,  r"$\epsilon_z$",        "epsilon_z")
plot_esf_def(gxy, r"$\gamma_{xy}$ [rad]", "gamma_xy")

# esfuerzos principales con sus orientaciones
plot_esf_def(s1,   r"$\sigma_1$ [Pa]",   "s1",    t1                     )
plot_esf_def(s2,   r"$\sigma_2$ [Pa]",   "s2",    t2                     )
plot_esf_def(tmax, r"$\tau_{máx}$ [Pa]", "tmax", [t1-np.pi/4, t1+np.pi/4])

# Fin :)