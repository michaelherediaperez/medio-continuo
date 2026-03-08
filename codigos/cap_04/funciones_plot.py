"""
En este archivo se presentan varias funciones que serán usadas para la 
graficación de los gráficos de colores de esfuerzos y deformaciones resultado 
del análisis de vigas sometidas a diferetnes configuraciones de apoyo y carga. 

Por: Michael Heredia Pérez
email: mherediap@unal.edu.co
Fecha: 2026-03-06

Universidad Nacional de Colombia sede Manizales 
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Para modificar el tamaño de la colorbar.
from mpl_toolkits.axes_grid1 import make_axes_locatable   

# Para centrar el color neutro en cero. 
from matplotlib.colors import TwoSlopeNorm

# Para gráficos 3D interactivos con Plotly.
# pip install plotly
import plotly.graph_objects as go


def save_image(fig, file_name, directory="results_plots", formats=("pdf","eps","png"), dpi=300):
    """
    Save the given Matplotlib figure in multiple formats.

    Args:
        fig (matplotlib.figure.Figure): 
            The Matplotlib figure to save.
        file_name (str): 
            The base name of the file (without extension).
        directory (str): Default 300
            The directory where the files will be saved.
        formats (tuple): Default ("pdf", "eps", "png").  
            The formats to save the figure.
        dpi (int): Default 300.
            The resolution of the saved images in dots per inch (DPI).
            
    Prints:
        Saving confirmation. "Figure saved to: <output_file>"
    """
    # Ensure the directory exists
    os.makedirs(directory, exist_ok=True)

    # Save the figure in all specified formats
    for fmt in formats:
        full_path = os.path.join(directory, f"{file_name}.{fmt}")
        fig.savefig(full_path, bbox_inches="tight", dpi=dpi)
        print(f"Figure saved to: {full_path}")
        

def plot_esf_def(x, y, variable, titulo, nombre=None, angulo=None, cmap="RdBu_r"):
    """
    Grafica esfuerzos/deformaciones sobre una malla (x,y) y ajusta el layout
    para que no quede apretado, sin tocar rcParams globales.

    Args:
        x (np.ndarray): 
            coordenadas x de la malla.
        y (np.ndarray): 
            coordenadas y de la malla.
        variable (np.ndarray): 
            campo a graficar (misma forma que x,y).
        titulo (str): 
            título del gráfico.
        nombre (str, default None): 
            Nombre de la figura al ser guardada (sin extensión). Si es None, 
            no se guarda la figura.
        angulo (np.ndarray | list[np.ndarray] | None): 
            ángulos para quiver.
        cmap (str): 
            colormap. 
            Recomendados: coolwarm, RdBu_r, vik (requiere ScientificColourMaps).
    """

    # Creo el lienzo.
    fig, ax = plt.subplots(figsize=(10, 3), constrained_layout=True)
    
    # Calcular rango real de los datos
    vmin = float(np.nanmin(variable))
    vmax = float(np.nanmax(variable))
    
    # Garantizar que vcenter=0 siempre quede ENTRE vmin y vmax,
    # incluso si todos los valores son del mismo signo.
    if vmin >= 0:
        variable = np.clip(variable, 0, None)  # Elimina residuos negativos espurios
        vmin = -1e-1
    if vmax <= 0:
        variable = np.clip(variable, None, 0)  # Elimina residuos positivos espurios
        vmax = 1e-1
        
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    # Coloreo el gráfico según la variable.
    im = ax.contourf(x, y, variable, 60, cmap=cmap, norm=norm)
    
    # Grafico las curvas de nivel.
    ax.contour(x, y, variable, 20, colors="black", linewidths=0.5,
           linestyles="solid", norm=norm)
    
    # Creo un axes adicional a la derecha en donde va la escala de colores.
    divider = make_axes_locatable(ax)
    
    # Creo la escala de colores.
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    # Configuro la escala para mostrar 5 intervales/4 etiquetas.
    fig.colorbar(im, cax=cax).ax.locator_params(nbins=3)
        
    # Para los esfuerzos principales se grafican las líneas que indiquen las
    # direcciones de los esfuerzos en cada nodo de la malla
    
    # Si se da un array para los ángulos. 
    if angulo is not None:
        # Se verifica que sea de tipo array numpy y se castea a lista. 
        if isinstance(angulo, np.ndarray):
            angulo = [angulo]
        for ang in angulo:
            ax.quiver(
                # Las coordenadas de la malla.
                x, y, 
                # Las direcciones.
                variable*np.cos(ang), variable*np.sin(ang), 
                # Configura cómo se grafican.
                headwidth=0, headlength=0, headaxislength=0, pivot='middle', 
                linewidths=2
                )
    
    # Se especifican los ejes y el título, y se colocan los ejes iguales.
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel(f"{titulo}\n$y$ [m]", rotation = 90)
    #ax.set_title(titulo)
    ax.set_aspect("equal", adjustable="box")   

    # Muestro la gráfica
    plt.show()

    # Almaceno el gráfico en formatos dados.
    if nombre:
        save_image(fig, nombre)
        

def plot_3d_esf_def(x, y, variable, titulo, nombre=None, cmap="RdBu_r"):
    """
    Grafica un campo escalar de esfuerzos/deformaciones como superficie 3D
    interactiva sobre la malla (x,y), con colormap centrado en cero.
    
    Args:
        x (np.ndarray): 
            coordenadas x de la malla.
        y (np.ndarray): 
            coordenadas y de la malla.  
        variable (np.ndarray): 
            campo escalar a graficar (misma forma que x,y).
        titulo (str): 
            título del gráfico y etiqueta del eje z.
        nombre (str, default None): 
            Nombre de la figura al ser guardada (sin extensión). Si es None, 
            no se guarda la figura.
        cmap (str): 
            colormap.
            Recomendados: coolwarm, RdBu_r, vik (requiere ScientificColourMaps).
            
    Notas:
        Requiere %matplotlib widget (ipympl) para interactividad en Jupyter.
    """

    # Normalización centrada en cero (mismo criterio que plot_esf_def)
    vmin = float(np.nanmin(variable))
    vmax = float(np.nanmax(variable))
    if vmin >= 0:
        vmin = -1e-10
    if vmax <= 0:
        vmax =  1e-10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    # Mapa de colores normalizado
    mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    colors   = mappable.to_rgba(variable)

    # Figura interactiva con %matplotlib widget (requiere ipympl)
    fig = plt.figure(figsize=(10, 5))
    ax  = fig.add_subplot(111, projection='3d')

    # Superficie
    ax.plot_surface(x, y, variable, facecolors=colors,
                    rstride=1, cstride=1, antialiased=True)

    # Colorbar
    mappable.set_array(variable)
    fig.colorbar(mappable, ax=ax, shrink=0.5, pad=0.1).ax.locator_params(nbins=3)

    # Etiquetas
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_zlabel(f"{titulo}")
    ax.set_title(titulo)

    plt.show()

    # Almaceno el gráfico en formatos dados.
    if nombre:
        save_image(fig, nombre)
    

def plot_3d_esf_def_interactive(x, y,variable, titulo, cmap="RdBu_r"):
    """
    Grafica un campo escalar de esfuerzos/deformaciones como superficie 3D
    interactiva usando Plotly, con colormap centrado en cero.
    
    Args:
        x (np.ndarray): 
            coordenadas x de la malla.
        y (np.ndarray): 
            coordenadas y de la malla.
        variable (np.ndarray): 
            campo escalar a graficar (misma forma que x,y).
        titulo (str): 
            título del gráfico y etiqueta del eje z.
        cmap (str): 
            colormap.
            Recomendados: RdBu_r, coolwarm, RdYlBu.
    """

    # Creo la superficie 3D con el colormap centrado en cero
    fig = go.Figure(data=[go.Surface(
        x=x, y=y, z=variable,
        colorscale=cmap,
        cmid=0,                 # Ancla el color neutro (blanco) en z=0
        colorbar=dict(nticks=4) # Número de divisiones en la escala de colores
    )])

    # Configuro el layout: título, etiquetas de ejes y tamaño de la figura
    fig.update_layout(
        title=titulo,
        scene=dict(
            xaxis_title="x [m]",
            yaxis_title="y [m]",
            zaxis_title=titulo,
        ),
        width=900, height=500
    )

    # Muestro la figura interactiva
    fig.show()
    

def plot_v_m(distancia, fuerza, titulo, nombre=None):
    """
    Grafica el diagrama unidimensional de fuerza cortante V(x) o 
    momento flector M(x) a lo largo de la viga.

    Args:
        distancia (np.ndarray): 
            Vector 1D con las coordenadas longitudinales x [m].
        fuerza (np.ndarray): 
            Vector 1D con los valores de la función interna evaluada 
            en cada punto de 'distancia'. 
        titulo (str): 
            Título del gráfico (admite sintaxis LaTeX).
        nombre (str, default None): 
            Nombre de la figura al ser guardada (sin extensión). Si es None, 
            no se guarda la figura.
    """
    
    fig, ax = plt.subplots(figsize=(10, 3), constrained_layout=True)
    
    ax.plot(distancia, fuerza, "-b")
    
    # ← Relleno en azul claro, sin distinción de signo
    ax.fill_between(distancia, fuerza, 0, color="b", alpha=0.15)
    
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel(r"$x$ [m]")
    ax.set_title(titulo)
    ax.grid(True)
    
    if nombre=="fuerza_cortante":
        ax.set_ylabel(r"$V(x)$ [kN]")
    if nombre=="momento_flector":
        ax.set_ylabel(r"$M(x)$ [kN$\cdot$m]")

    plt.show()
    
    # Almaceno el gráfico en formatos dados.
    if nombre:
        save_image(fig, nombre)