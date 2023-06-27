import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 
import imageio


df = pd.read_csv('Posiciones4.txt', sep=" ",header=None).to_numpy().transpose() # bordes 
dg = pd.read_csv('mover4.txt', sep=" ",header=None).to_numpy().transpose()      # las pepitas
 
npepas = 4

Pepas = [[] for i in range(npepas)]

for i in range(len(dg[0])):
    Pepas[i%npepas].append([dg[0][i],dg[1][i]])
# print(Pepas[0])
Pepas = np.array(Pepas)

################################# dim

# plt.plot(df[0][2:],df[1][2:],'.',markersize = 2)
# plt.xlim(-1,602)
# plt.ylim(-1,205)
# for i in range(int(3*(len(dg[0]))/128/npepas)):
# for i in range(int((len(dg[0]))/npepas)):
for i in range(10000):
    if i%200 == 0:
        plt.plot(Pepas[0][i][0],Pepas[0][i][1],'.')
        plt.plot(Pepas[1][i][0],Pepas[1][i][1],'.')
        plt.plot(Pepas[2][i][0],Pepas[2][i][1],'.')
        plt.plot(Pepas[3][i][0],Pepas[3][i][1],'.')
        plt.plot(df[0][npepas:],df[1][npepas:],'.',markersize = 2)
        plt.xlim(-1,602)
        plt.ylim(-1,205)
        plt.savefig('gif/gif'+str(i)+'.png')
        plt.close()        


path = 'gif/'
archivo = sorted(os.listdir(path))
img_array = []

for x in range(0,len(archivo)):
    nomArchivo = archivo[x]
    dirArchivo = path + str(nomArchivo)
    
    leer_imagen = imageio.imread(dirArchivo)
    img_array.append(leer_imagen)

imageio.mimwrite('Elgif.gif',img_array, 'GIF' , duration = 1)

################################# Adim

# # plt.plot(df[0][2:],df[1][2:],'.',markersize = 2)
# # plt.xlim(-0.05,10.02)
# # plt.ylim(-0.1,3.5)

# # for i in range(int((len(dg[0])-100000)/npepas)):
# for i in range(int(len(dg[0])/npepas-500/npepas)):
#     if i%10 == 0:
#         plt.plot(Pepas[0][i][0],Pepas[0][i][1],'.')
#         plt.plot(Pepas[1][i][0],Pepas[1][i][1],'.')
#         # plt.plot(Pepas[2][i][0],Pepas[2][i][1],'.')
#         # plt.plot(Pepas[3][i][0],Pepas[3][i][1],'.')
#         plt.plot(df[0][npepas:],df[1][npepas:],'.',markersize = 2)
#         plt.xlim(-0.05,10.02)
#         plt.ylim(-0.1,3.5)
#         plt.show()
#         plt.close()