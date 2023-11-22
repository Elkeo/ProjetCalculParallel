#include <stdlib.h>
#include <stdio.h>
#include <math.h>




int Nx=25,Ny=25;
double Lx=1,Ly=1,D=1;




double *produit(double alpha,double beta,double gamma,double*x){

int n=Nx*Ny;
double *prod;
prod=(double*)malloc((n)*sizeof(double));



for (int i = 0; i < Nx; i++) {
    if (i == 0) {
        prod[i] = alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    } else if (i == Nx - 1) {
        prod[i] = beta * x[i - 1] + alpha * x[i] + gamma * x[i + Nx];
    } else {
        prod[i] = beta * x[i - 1] + alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    }
}

for (int i = Nx; i < n - Nx; i++) {
    if (i % Nx == 0) {
        prod[i] = gamma * x[i - Nx] + alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    } else if (i % Nx == Nx-1) {
        prod[i] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i] + gamma * x[i + Nx];
    } else {
        prod[i] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    }
}

for (int i = n - Nx; i < n; i++) {
    if (i == n - Nx) {
        prod[i] = gamma * x[i - Nx] + alpha * x[i] + beta * x[i + 1];
    } else if (i == n - 1) {
        prod[i] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i];
    } else {
        prod[i] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i] + beta * x[i + 1];
    }
}

return prod;
}


double * Sousvect(double *x,double *y,int n)
{
    double *s;
     s=(double*)malloc(n*sizeof(double));
     for(int i=0;i<n;i++){
         s[i]=x[i]-y[i];
     }
     return s;
}


double * Sommevect(double *x,double *y,int n)
{
    double *s;
     s=(double*)malloc(n*sizeof(double));
     for(int i=0;i<n;i++){
         s[i]=x[i]+y[i];
     }
     return s;
}


double produit_scalaire(double *X, double *Y,int N)
{
    double ps;
    ps=0;
    for(int i=0;i<N;i++)
    {
        ps=ps+X[i]*Y[i];


    }
    return ps;
}


double * multi_scalaire(double a,double *x,int n)
{
    double *s;
    s=(double*)malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
        s[i]=a*x[i];
    }
    return s;
}


double * gradc(double alpha, double beta, double gama, double *source, double eps, int kmax, int Nx, int Ny)
{
    double *x, *r, *z, *p, *r0;
    double a, re, b, y;
    int k, n = Nx * Ny;
    r = (double*)malloc((n) * sizeof(double));
    z = (double*)malloc((n) * sizeof(double));
    p = (double*)malloc((n) * sizeof(double));
    r0 = (double*)malloc((n) * sizeof(double));
    x = (double*)malloc((n) * sizeof(double));

    for(int i = 0; i < n; i++)
    {
        x[i] = 0;
    }

    r = Sousvect(source, produit(alpha, beta, gama, x), n);
    p = r;
    re = sqrt(produit_scalaire(r, r, n));
    k = 0;

    while(k < kmax && re > eps)
    {
        re = sqrt(produit_scalaire(r, r, n));
        z = produit(alpha, beta, gama, p);
        a = produit_scalaire(r, r, n) / produit_scalaire(z, p, n);
        
        x = Sommevect(x, multi_scalaire(a, p, n), n);
        r0 = r;
        r = Sousvect(r0, multi_scalaire(a, z, n), n);
        b = produit_scalaire(r, r, n) / produit_scalaire(r0, r0, n);
        p = Sommevect(r, multi_scalaire(b, p, n), n);
        k = k + 1;
    }
   printf("La pente est %f",a);
    return x;
}

double f(double x,double y,double t, int cas)
{
  if (cas == 1) {
    return 2*(y-pow(y,2)+x-pow(x,2));
  }

  else if (cas == 2) {
     return sin(x)+cos(y);
  }

  else if (cas == 3) {
    double pi=4*atan(1.);
       return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos((pi*t)/2);
  }
}

double g(double x,double y,double t, int cas)
{
  if (cas == 1) {
     return 0;
  }

  else if (cas == 2) {
    return sin(x)+cos(y);
  }

  else if (cas == 3) {
     return 0;
  }
}


double h(double x,double y,double t, int cas)
{
  if (cas == 1) {
      return 0;
  }

  else if (cas == 2){

     return sin(x)+cos(y);
  }

  else if (cas == 3) {
     return 1;
  }
}
double norme_infinie(double *w, int n) {
    double norme = 0.0;
    for (int i = 0; i < n; i++) {
        double valeur_absolue = fabs(w[i]);
        if (valeur_absolue > norme) {
            norme = valeur_absolue;
        }
    }
    //printf("%lf\n",norme);
    return norme;
}






double erreur(double *w,double *v,int n){

    double *z;
    z = (double*)malloc(n * sizeof(double));

    z=Sousvect(w,v,n);
    

return norme_infinie(z,n)/norme_infinie(v,n);


}
double* uex(int cas,double dx,double dy){
    int n = Nx * Ny;
    double *uex;
    uex = (double*)malloc(n * sizeof(double));

   if(cas==2){
   for(int i=1;i<=Nx;i++) {
        for(int j=1;j<=Ny;j++) {
            int I = (j-1)*Nx + i-1;
            uex[I]=sin(i*dx)+cos(j*dy);
        }        
    }
   }
   else if(cas==1){

       for(int i=1;i<=Nx;i++) {
        for(int j=1;j<=Ny;j++) {
            int I = (j-1)*Nx + i-1;
            uex[I]=(i*dx)*(1-i*dx)*j*dy*(1-j*dy);
        }    
     

    } 



   }
   else{
           for(int i=1;i<=Nx;i++) {
        for(int j=1;j<=Ny;j++) {
            int I = (j-1)*Nx + i-1;
            uex[I]=0;
        }   
   }
   }
    return uex;
}

double* Source(double dt, double t, double *u, int cas)
{
    double *b;

    double dx, dy;
    dx = Lx / (Nx + 1);
    dy = Ly / (Ny + 1);
    int n = Nx * Ny;
    b = (double*)malloc(n * sizeof(double));
    int i, j, I;

    for (j = 1; j <= Ny; ++j) {
        if (j == 1) {
            for (i = 1; i <= Nx; ++i) {
                I = (j - 1) * Nx + i - 1;
                b[I] = u[I] + dt * (f(i * dx, j * dy, t, cas) + D / (dy * dy) * g(i * dx, 0, t, cas));

                if (i == 1)
                    b[I] += dt * (D / (dx * dx) * h(0, j * dy, t, cas));
                if (i == Nx)
                    b[I] += dt * (D / (dx * dx) * h(Lx, j * dy, t, cas)); 
            }
        }
        if (j >= 2 && j < Ny) {
            for (i = 1; i <= Nx; ++i) {
                I = (j - 1) * Nx + i - 1;
                b[I] = u[I] + dt * f(i * dx, j * dy, t, cas);

                if (i == 1)
                    b[I] += dt * (D / (dx * dx) * h(0, j * dy, t, cas));
                if (i == Nx)
                    b[I] += dt * (D / (dx * dx) * h(Lx, j * dy, t, cas));
            }
        }
        if (j == Ny) {
            for (i = 1; i <= Nx; ++i) {
                I = (j - 1) * Nx + i - 1;
                b[I] = u[I] + dt * (f(i * dx, j * dy, t, cas) + D / (dy * dy) * g(i * dx, Ly, t, cas));

                if (i == 1)
                    b[I] += dt * (D / (dx * dx) * h(0, j * dy, t, cas));
                if (i == Nx)
                    b[I] += dt * (D / (dx * dx) * h(Lx, j * dy, t, cas));
            }
        }
    }

    return b;
}



int main(int argc, char** argv) {
    int cas=2,nb_iterations=10,kmax=100000;

    double *w,*u,*solexacte;
    double *z;
    int n=Nx*Ny;
    double eps=0.000001,dt,tn,tf=2;
    z=(double*)malloc((n)*sizeof(double));
    w=(double*)malloc((n)*sizeof(double));
    u=(double*)malloc((n)*sizeof(double));



    double dx,dy;
    dx=Lx/(Nx+1);
    dy=Ly/(Ny+1);
    solexacte=uex(cas,dx,dy);
    dt=tf/nb_iterations;
    printf("dt=%f\n",dt);

    double alpha,beta,gama;
    alpha=1+dt*D*((2./(pow(dx,2))+2./pow(dy,2)));
    beta=(-1./pow(dx,2))*D*dt;
    gama=(-1./pow(dy,2))*D*dt;

    printf("alpha=%f\n",alpha);
    for(int i=1;i<=Nx;i++) {
        for(int j=1;j<=Ny;j++) {
            int I = (j-1)*Nx + i-1;
            u[I]=sin(i*dx)+cos(j*dy);
        }        
    }

    
 
    for(int ni=0; ni<nb_iterations; ni++) {
        tn=(ni+1)*dt;
        z=Source(dt,tn,u,cas);
        w=gradc(alpha,beta,gama,z,eps,kmax,Nx,Ny);
        u=w;
        printf("L'erreur à la %d ieme itération est %e\n",ni,erreur(w,solexacte,n));
      


        FILE* fichier;
        char file_name[256];
        sprintf(file_name,"sol.%d.dat",ni);
        fichier=fopen(file_name, "w");
        for (int I=0; I<n; I++) {
            int i=I%Nx+1;
            int j=I/Nx+1;
            fprintf(fichier, "%f %f %f %f \n",i*dx, j*dy, w[I],solexacte[I]);
        }
        fclose(fichier);
    }

    return 0;
}
   