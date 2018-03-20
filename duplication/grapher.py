#!/usr/bin/python3

import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

def leer_ds(archivo_ds):
    """Cuenta la frecuencia de DS en el genoma.

    Argumentos:
    archivo_ds -- El archivo donde se encuentra el DS

    Regresa:
    freq -- Una tabla hash o diccionario con las frecuencia
    de cada uno de los DS's
    """
    
    freq = {}
    with open(archivo_ds, 'r') as f:
        for line in f:
            _, _, ds = line.split()
            freq[ds] = freq.get(ds, 0) + 1
    return freq

def plot(freq, title):
    """Gráfica las frecuencias de los DS. """

    # Ordena las frecuencias en orden ascdente, regresa una lista de tuplas.

    lists = sorted(freq.items())

    # la lista de pares las desempaquetamos en dos tuplas.

    x, y = zip(*lists)
    
    #plt.title(title)
    plt.xlabel("dS")
    plt.ylabel("Frequency")
    
    plt.plot(x[:], y[:], color="#3F5D7D")
    plt.savefig(title + ".png", transparent=False)

def run(path):
    print('Qué pedo con esto')
    titulo, _ = path.split('.')
    freq = leer_ds(path)
    plot(freq, titulo)

if __name__ == '__main__':
    archivo = sys.argv[1]
    titulo, _ = archivo.split('.')
    freq = leer_ds(archivo)
    plot(freq, titulo)
