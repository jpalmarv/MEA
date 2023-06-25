import numpy as np

Lx, Ly=600, 205
Nb,Ns=1,1
N=Nb+Ns

Nsierras = 10  

H = Ly*0.2
L = Lx/Nsierras
d = Ly - 3*H*0.5
R0,M0 = 2.4*1.25,1000

cAbajo = int(np.sqrt(L*L+H*H)*0.5*0.5/R0)
cArriba = int(np.sqrt(L*L+H*H)*0.5/R0)
lAbajo = int(H*0.5*0.5/R0)-1
lArriba = int(H*0.5/R0)-1

Npepas=Nsierras*(cArriba+2*cAbajo)+(2*Nsierras-1)*lAbajo+(Nsierras-1)*lArriba
Ntot = N + Npepas + 2

print("const int Nsierras{",Nsierras,"}; // # sierras arriba",
       "\nconst double H{",H,"}; // Altura, longitud pico de arriba",
       "\nconst double L{",L,"};",
       "\nconst double d{",d,"}; // Separación picos de arriba y abajo",
       "\nconst double R0{",R0,"},M0{",M0,"}; // Radio y masa esferas de picos",
       "\n",
       "\nconst int cAbajo =",cAbajo,"; // # pepas en una diagonal de abajo",
       "\nconst int cArriba =",cArriba,"; // # pepas en una diagonal de arriba",
       "\nconst int lAbajo = ",lAbajo,"; // # pepas en una linea de abajo",
       "\nconst int lArriba = ",lArriba,"; // # pepas en una linea de arriba",
       "\nconst int Npepas=",Npepas,"; // # total de pepas en sierras",
       "\nconst int Ntot =",Ntot,"; // N pepas en líquido, Npepas pepas en sierras, 2 paredes)")

# //---- declarar constantes ---
# const double K=1.0e4;               // Constantes de simulación
# const double Lx=600, Ly=205;        // Tamaño de espacio de simulación (~medidas del paper)
# const int Nb=1, Ns=0, N=Nb+Ns;      // Nb = bolas grandes, Ns = bolas pequeñas, N = Total de coloides
# const double g=9.8, Gamma=20, Kcundall=500, mu=0.4;  /// Gravedad y otras cosas (***************Gamma, Kcundall, mu = ?????)

# // Número de pepas en sierras
# const int Nsierras{10};                   // # sierras arriba
# const double H{Ly*0.2};    // Altura, longitud pico de arriba
# const double L{Lx/Nsierras};
# const double d{Ly - 3*H*0.5};             // Separación picos de arriba y abajo
# const double R0{1.25},M0{1000.0};         // Radio y masa esferas de picos

# const int cAbajo = (sqrt(L*L+H*H)*0.5*0.5/R0);    // # pepas en una diagonal de abajo
# const int cArriba = (sqrt(L*L+H*H)*0.5/R0);       // # pepas en una diagonal de arriba
# const int lAbajo = (H*0.5*0.5/R0);                // # pepas en una linea de abajo
# const int lArriba = (H*0.5/R0);                   // # pepas en una linea de arriba

# const int Npepas=Nsierras*(cArriba+2*cAbajo)+(2*Nsierras-1)*lAbajo+(Nsierras-1)*lArriba; // # total de pepas en sierras

# const int Ntot = N + Npepas + 2; // N pepas en líquido, Npepas pepas en sierras, 2 paredes
