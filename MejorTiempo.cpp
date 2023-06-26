#include <iostream>     // General salida y entrada de datos
#include <cmath>        // Mates
#include "vector.h"     // Vectores
#include "Random64.h"   // Números aleatorios
#include <fstream>      // Salida a txt
using namespace std;

//---- declarar constantes ---
const double K=1.0e4;               // Constantes de simulación (Cambiar con viscosidad)
const double Lx=600, Ly=205;        // Tamaño de espacio de simulación (~medidas del paper)
const int Nb=2, Ns=2, N=Nb+Ns;      // Nb = bolas grandes, Ns = bolas pequeñas, N = Total de coloides
const double g=9.8, Gamma=150, Kcundall=200, mu=1.4;  /// Gravedad y otras cosas (***************Gamma, Kcundall, mu = ?????)
const double Reff = 0.119;
const double Delta = 1;

// Número de pepas en sierras

const int Nsierras{ 10 }; // # sierras arriba 
const double H{ 41.0 }; // Altura, longitud pico de arriba 
const double L{ 60.0 }; 
const double d{ 143.5 }; // Separación picos de arriba y abajo 
const double R0{ 3.0 },M0{ 10000 }; // Radio y masa esferas de picos 
 
const int cAbajo = 6 ; // # pepas en una diagonal de abajo 
const int cArriba = 12 ; // # pepas en una diagonal de arriba 
const int lAbajo =  2 ; // # pepas en una linea de abajo 
const int lArriba =  5 ; // # pepas en una linea de arriba 
const int Npepas= 323 ; // # total de pepas en sierras 
const int Ntot = Npepas + 2 + N; // N pepas en líquido, Npepas pepas en sierras, 2 paredes)

//---- Constantes de PEFRL ----
const double epsilon=0.1786178958448091e00; 
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

 

//--- declarar clases -----
class Cuerpo;                       // Clase de cuerpo (pepitas)
class Colisionador;                 // Clase de colisión (interacción)



//---- implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,I;             // r,V,F (F incluye el torque)
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
            double theta0,double omega0,double m0,double R0);       // Inicializa los circulos con velocidades y pocisiones dadas
  void BorreFuerza(){F.load(0,0,0);tau=0;};                         // Borra fuerzas e inicializa tau (tau maybe torque)
  void AdicioneFuerza(vector3D F0,double tau0){F+=F0;tau+=tau0;};   // Suma fuerzas y torques en el elemento
  void Mueva_r(double dt, double Coeficiente);                      // Cositas de PERFL para mover partículas
  void Mueva_V(double dt, double Coeficiente);
  double Getx(void){return r.x();}; //inline                        // Retornar posiciones
  double Gety(void){return r.y();}; //inline
  friend class Colisionador;                                        // Npi qué hace friends
};
//--- funciones de Cuerpo ---
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,  
		    double theta0,double omega0,double m0,double R0){       
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;                        // Da las posiciones y calcula el momento de inercia (esfera 3D)
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt); theta+=omega*(Coeficiente*dt);             // Cositas de PERFL para mover partículas
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); omega+=tau*(Coeficiente*dt/I);
}


//--- clase Colisionador ----
class Colisionador{
private:                                                            
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];                     // Cundall para las esferas, sold?******************** 
public:                                                             
  void Inicie(void);                                                // Npi
  void CalculeFuerzas(Cuerpo * Grano,double dt, int & value,double mu_al, double sigma_al,Crandom & RAN64);                   // Calcula fuerzas del grano en el tiempo dt
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,          // Fuerzas entre granos (Npi)
			  double & x_Cundall,double & s_old,double dt);
};

//--- funciones de Colisionador ---
void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
}

void Colisionador::CalculeFuerzas(Cuerpo * Grano,double dt,int & value, double FAx, double FAy,Crandom & RAN64){ 
  int i,j; double F0x, F0y;
  //Se crea una variable para que luego en ciclo for de despues calcule la fuerza aleatoria
  /* se crea mu_al y sigma_al*/
  //--- Borrar todas las fuerzas ---
  for(i=0;i<Ntot;i++)
    Grano[i].BorreFuerza();                                                 // Borra fuerzas de todo lado
  //--- Añadir fuerza de gravedad ---
  /*vector3D F0;
  for(i=0;i<N;i++){                                                         // No pone fuerzas en paredes (?)
    F0x = -Grano[i].V.x();
    if(a > 0){value*=-1;}
    F0y = -Grano[i].V.y()-value*pow((Grano[i].R/Reff),3);
    F0.load(F0x,F0y,0);
    Grano[i].AdicioneFuerza(F0,0);                                          // Fuerza de gravedad (toca cambiar por fuerza periódica)
  }*/
  //--- Calcular Fuerzas entre pares de planetas ---
  for(i=0;i<N;i++)                                                          // No calcula fuerzas en paredes (?)
    for(j=i+1;j<Ntot;j++)      
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);   // Calcula todas parejas pero solo pone las que se mueven
      
   //--- Calcular Fuerza aleatoria  ----
   vector3D FA;
   for (i=0; i<N; i++){
		FA.load(FAx, FAy,0);
		Grano[i].AdicioneFuerza(FA,0);
   }
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
				      double & x_Cundall,double & s_old,double dt){         // Fuerzas entre dos granos (?)
  //Cantidades generales para saber si hay contacto
  vector3D r21=Grano2.r-Grano1.r; double R1=Grano1.R, R2=Grano2.R;
  double d=r21.norm(); double s=R1+R2-d;                                    // Claro, constantes de contacto (vector entre ambos centers)
  
  if(s>0){//si hay contacto
    //Variables a calcular
    vector3D F2,F1,tau2,tau1;                                               // fuerzas y torques de ambos granos
    //Vectores unitarios
    vector3D n=r21*(1.0/d),t,k; t.load(n.y(),-n.x(),0); k.load(0,0,1);      // unitarios, n = fuerza central, t = torque uni, k = z
    //Velocidades relativas
    vector3D V21=Grano2.V-Grano1.V;
    vector3D Rw; Rw.load(0,0,R1*Grano1.omega+R2*Grano2.omega);              // Npi, buscar
    vector3D Vc=V21-(Rw^n); double Vn=Vc*n, Vt=Vc*t;                        // Npi, cosas raras maybe de elementos finitos
    //Fn (Fuerza de Hertz-Kuwabara-Kono) 
    double m12=(Grano1.m*Grano2.m)/(Grano1.m+Grano2.m);                     // Fuerza de entrar en grano
    double Fn=K*pow(s,1.5)-Gamma*m12*sqrt(s)*Vn; if(Fn<0) Fn=0;
    
    //Ft (Fuerza de Cundall)
    x_Cundall+=Vt*dt; double Ft=-Kcundall*x_Cundall; double Ftmax=mu*fabs(Fn); // Fricción? Toca ponerla? Maybe?
    if(fabs(Ft)>Ftmax) Ft=Ft/fabs(Ft)*Ftmax;
    
    //Calcula y Cargue las fuerzas
    F2=n*Fn+t*Ft; tau2=((n*(-R2))^F2); F1=F2*(-1); tau1=((n*R1)^F1);        // Carga las fuerzas
    Grano2.AdicioneFuerza(F2,tau2*k);   Grano1.AdicioneFuerza(F1,tau1*k);
  }

  if(s_old>=0 && s<0) x_Cundall=0;                                          // Npi
    
  s_old=s;                                                                  // S antiguo, makes sense, qué tan dentro están
}

void Tiempo(double & t, int & val,double Tmas){
    if (t<Tmas){val = 1;}
    else if (t<2*Tmas+Delta){val = -1;}
    else {t -= 2*Tmas+Delta, val = 1;}
}



// Main 
int main(void){
    Cuerpo Molecula[Ntot];
    Colisionador Hertz; Hertz.Inicie();
    Crandom ran64(1);
    Crandom RAN64(1);
    double m0{M0/100},r0{4.2};
    double Omega0,OmegaMax=8.0;
    double cuadros=10,t,dt=1e-3,tmax=cuadros*sqrt(Ly/g),tcuadro=tmax/(4*cuadros),tdibujo,tPaso;
    
    Molecula[0].Inicie(570.0,100.0,0,0,0,0,m0,r0);
    Molecula[1].Inicie(350.0,100.0,0,0,0,0,m0,3.0);
    Molecula[2].Inicie(250.0,100.0,0,0,0,0,m0,r0);
    Molecula[3].Inicie(50.0,100.0,0,0,0,0,m0,3);

    int i{0},j{0};
    double x0 = 2*R0*L/sqrt(L*L+H*H);
    double pendiente = H/L;
    for(i = 0; i < 2*Nsierras; i++){ // Se multiplica Nsierra por 2 porque estamos abajo
        
        for(j = 0; j < cAbajo;j++){
            Molecula[N+i*cAbajo+j].Inicie(x0*j+i*L*0.5,pendiente*x0*j,0,0,0,0,M0,R0);
        }

        for(j = 0; j < lAbajo;j++){
            if(i!=2*Nsierras-1)
            Molecula[N+cAbajo*2*Nsierras+cArriba*Nsierras+j+lAbajo*i].Inicie(i*L*0.5+L*0.5,2*R0*j+2*R0,0,0,0,0,M0,R0);
        }
        if(i<Nsierras){
            for(j = 0; j < cArriba;j++){
                Molecula[N+cAbajo*2*Nsierras+i*cArriba+j].Inicie(x0*j+i*L,pendiente*x0*j+d+H*0.5,0,0,0,0,M0,R0);
            }
            for(j = 0; j < lArriba;j++){
                if(i!=Nsierras-1)
                Molecula[N+cAbajo*2*Nsierras+cArriba*Nsierras+lAbajo*(2*Nsierras-1)+j+lArriba*i].Inicie(i*L+L,2*R0+2*R0*j+d+H*0.5,0,0,0,0,M0,R0);
            }
        }
        
    }
    
    //Inicializar las paredes
    double Rpared=10000*Lx;
    //------------------(  x0,       y0,Vx0,Vy0,theta0,omega0,m0,R0)
    Molecula[Ntot-2].Inicie(Lx+Rpared,Ly/2,  0,  0,     0,     0,3*M0,Rpared); //Pared derecha
    Molecula[Ntot-1].Inicie(  -Rpared,Ly/2,  0,  0,     0,     0,3*M0,Rpared); //Pared izquierda

    ofstream myfile;
    myfile.open ("Posiciones4.txt");
    for (i = 0; i< Ntot; i++){
        myfile << Molecula[i].Getx() << " " << Molecula[i].Gety() << '\n';
    }
    myfile.close();
   
    myfile.open ("mover4.txt");

    
    int value{-1};
    double Tmas{0.5};
    
    for(t=0,tdibujo=tcuadro +1,tPaso = Tmas ; t<2*tmax; t+=dt,tdibujo+=dt,tPaso+=dt){
        if(tdibujo>tcuadro){
            for (i = 0; i< N; i++){
            myfile << Molecula[i].Getx() << " " << Molecula[i].Gety() << '\n';
            }
        }

        Tiempo(tPaso, value,Tmas);
        
        //--- Muevase por PEFRL ---
        for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);
        Hertz.CalculeFuerzas(Molecula,dt,value,RAN64.gauss(0,2000),RAN64.gauss(0,2000),RAN64);
        for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
        for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
        Hertz.CalculeFuerzas(Molecula,dt,value,RAN64.gauss(0,2000),RAN64.gauss(0,2000),RAN64);
        for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
        for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chiepsilon);
        Hertz.CalculeFuerzas(Molecula,dt,value,RAN64.gauss(0,2000),RAN64.gauss(0,2000),RAN64);
        for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
        for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
        Hertz.CalculeFuerzas(Molecula,dt,value,RAN64.gauss(0,2000),RAN64.gauss(0,2000),RAN64);
        for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
        for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon); 
   
    }
    myfile.close();


    return 0;
}





