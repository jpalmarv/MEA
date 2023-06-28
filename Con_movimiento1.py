import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 
import imageio


df = pd.read_csv('Posiciones4.txt', sep=" ", header=None).to_numpy().transpose() # bordes 
dg = pd.read_csv('mover4.txt', sep=" ", header=None).to_numpy().transpose()      # las pepitas
 
npepas = 4

Pepas = [[] for i in range(npepas)]

for i in range(len(dg[0])):
    Pepas[i%npepas].append([dg[0][i], dg[1][i]])

Pepas = np.array(Pepas)

################################# dim

for i in range(10000):
    if i % 200 == 0:
        plt.plot(Pepas[0][i][0], Pepas[0][i][1], '.')
        plt.plot(Pepas[1][i][0], Pepas[1][i][1], '.')
        plt.plot(Pepas[2][i][0], Pepas[2][i][1], '.')
        plt.plot(Pepas[3][i][0], Pepas[3][i][1], '.')
        plt.plot(df[0][npepas:], df[1][npepas:], '.', markersize=2)
        plt.xlim(-1, 602)
        plt.ylim(-1, 205)
        print(i)
        plt.savefig('gif/' + str(i) + '.png')
        plt.close()        

path = 'gif/'

def obtener_numero(nombre):
    numero = ''
    for caracter in nombre:
        if caracter.isdigit():
            numero += caracter
        else:
            break
    return int(numero)

archivos_ordenados = sorted(os.listdir(path), key=lambda x: obtener_numero(x))

img_array = []

for x in range(1, len(archivos_ordenados)):
    nomArchivo = archivos_ordenados[x]
    dirArchivo = path + str(nomArchivo)
    
    leer_imagen = imageio.v2.imread(dirArchivo)
    img_array.append(leer_imagen)

imageio.v2.mimwrite('Elgif.gif', img_array, 'GIF', duration=0.7) # duración original 1 

print(archivos_ordenados)
