#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include <string.h>





void charge_c(int me,int n,int np,int *ibeg,int *iend)
{
  int r=n%np;
  if(me<r){
    *ibeg=(me)*(n/np+1);
    *iend=(me+1)*(n/np+1)-1;
  }
  else{
    *ibeg=r+me*(n/np);
    *iend=*ibeg+n/np-1;
  }
}

double *fusion12(double *u,double *v,int n,int m,int ibeg,int iend)
{
    double *z;
    z=(double*)calloc((n*m),sizeof(double));
    for(int i=ibeg;i<=n+iend;i++)
    {
        if(i<=iend)
        {
            z[i]=u[i-ibeg];
        }
        else if(i>=iend+1 && i<=n+iend)
        {
            z[i]=v[i-iend-1];
        }
    }
 
return z;

}




double *fusion22(double *u,double *v,int n,int m,int ibeg,int iend)
{
    double *z;
    z=(double*)calloc((n*m),sizeof(double));
    for(int i=ibeg-n;i<=iend;i++)
    {
        
         if(i>=ibeg-n && i<ibeg)
        {
            z[i]=v[i-ibeg+n];
        }
        else if(i>=ibeg)
        {
            z[i]=u[i-ibeg];
        }
    }

    return z;

}






double* fusion32(double* u, double* v, double* w, int n, int m, int ibeg, int iend) {
    double* z = calloc(n*m,sizeof(double));

    for (int i = ibeg - n; i <= iend + n; i++) {
        if (i < ibeg)
            z[i] = v[i - ibeg + n];
        else if (i <= iend)
            z[i] = u[i - ibeg];
        else if (i <= n + iend)
            z[i] = w[i - iend - 1];
    }
    
    return z;
}




double *prodloc(double alpha,double beta,double gamma,double*x,int Nx,int Ny,int ibeg,int iend){

int n=Nx*Ny;
double *prod;
prod=(double*)malloc((iend-ibeg+1)*sizeof(double));

 for(int i=ibeg;i<=iend;i++)

 {

if( i >= 0 && i < Nx) {
    if (i == 0) {
        prod[i-ibeg] = alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    } else if (i == Nx - 1) {
        prod[i-ibeg] = beta * x[i - 1] + alpha * x[i] + gamma * x[i + Nx];
    } else {
        prod[i-ibeg] = beta * x[i - 1] + alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    }
}

if ( i >= Nx && i < n - Nx) {
    if (i % Nx == 0) {
        prod[i-ibeg] = gamma * x[i - Nx] + alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    } else if (i % Nx == Nx-1) {
        prod[i-ibeg] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i] + gamma * x[i + Nx];
    } else {
        prod[i-ibeg] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i] + beta * x[i + 1] + gamma * x[i + Nx];
    }
}

if(i >= n - Nx && i < n) {
    if (i == n - Nx) {
        prod[i-ibeg] = gamma * x[i - Nx] + alpha * x[i] + beta * x[i + 1];
    } else if (i == n - 1) {
        prod[i-ibeg] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i];
    } else {
        prod[i-ibeg] = gamma * x[i - Nx] + beta * x[i - 1] + alpha * x[i] + beta * x[i + 1];
    }
}


}

return prod;


}











double *produitmat(double a,double b,double c,double *y,int nproc,int rank,int ibeg,int iend,int Nx,int Ny){



double *v,*w,*x,*z,*x1,*x2,*z1,*z2;





if(rank==0){
//ne pas utiliser N dernier N premier pointer sur le premier element direct

z=(double*)malloc((Nx)*sizeof(double));

w=(double*)malloc((Nx+iend-ibeg+1)*sizeof(double));
//w=(double*)malloc((Nx*Ny)*sizeof(double));
v=(double*)malloc((iend-ibeg+1)*sizeof(double));

if(nproc>1){


MPI_Send(&y[iend-ibeg+1-Nx],Nx, MPI_DOUBLE,rank+1,0, MPI_COMM_WORLD);
MPI_Recv(z,Nx, MPI_DOUBLE,rank+1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);



w=fusion12(y,z,Nx,Ny,ibeg,iend);
v=prodloc(a,b,c,w,Nx,Ny,ibeg,iend);



}
else {
v=prodloc(a,b,c,y,Nx,Ny,ibeg,iend);
}





free(w);
free(z);


return v;
}
else if (rank ==nproc-1){

z=(double*)malloc((Nx)*sizeof(double));


z1=(double*)malloc((Nx)*sizeof(double));
z2=(double*)malloc((Nx)*sizeof(double));
if(rank%2==0){

MPI_Send(&y[0],Nx, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD);
MPI_Recv(z,Nx, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
else {

MPI_Recv(z,Nx, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Send(&y[0],Nx, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD);
}

w=(double*)malloc((Nx+iend-ibeg+1)*sizeof(double));
//w=(double*)malloc((Nx*Ny)*sizeof(double));
v=(double*)malloc((iend-ibeg+1)*sizeof(double));

w=fusion22(y,z,Nx,Ny,ibeg,iend);
v=prodloc(a,b,c,w,Nx,Ny,ibeg,iend);




free(w);
free(z);

free(z1);
free(z2);


return v;
}
else{

 z=(double*)malloc((Nx)*sizeof(double));


z1=(double*)malloc((Nx)*sizeof(double));
z2=(double*)malloc((Nx)*sizeof(double));
    
//pour les maillages fin commbloquantes soit ca soit les free.
    if(rank%2==0){

   MPI_Send(&y[iend-ibeg+1-Nx],Nx, MPI_DOUBLE,rank+1, 0, MPI_COMM_WORLD);
   MPI_Send(&y[0],Nx, MPI_DOUBLE,rank-1, 0, MPI_COMM_WORLD);
   MPI_Recv(z1,Nx, MPI_DOUBLE,rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(z2,Nx, MPI_DOUBLE,rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
}
else{
 
   MPI_Recv(z1,Nx, MPI_DOUBLE,rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(z2,Nx, MPI_DOUBLE,rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Send(&y[iend-ibeg+1-Nx],Nx, MPI_DOUBLE,rank+1, 0, MPI_COMM_WORLD);
   MPI_Send(&y[0],Nx, MPI_DOUBLE,rank-1, 0, MPI_COMM_WORLD);
}

w=(double*)malloc((2*Nx+iend-ibeg+1)*sizeof(double));
//w=(double*)malloc((Nx*Ny)*sizeof(double));
v=(double*)malloc((iend-ibeg+1)*sizeof(double));
         
       w=fusion32(y,z1,z2,Nx,Ny,ibeg,iend);//fusion3 marche
       v=prodloc(a,b,c,w,Nx,Ny,ibeg,iend);


free(w);
free(z);
free(z1);
free(z2);
return v;




}

}




double * Sousvect(double *x,double *y,int n)
{
    double *h;
     h=(double*)malloc(n*sizeof(double));
     for(int i=0;i<n;i++){
         h[i]=x[i]-y[i];
     }
       
     return h;
   
    
}



double * Sommevect(double *x,double *y,int n)
{
    double *h;
     h=(double*)malloc(n*sizeof(double));
     for(int i=0;i<n;i++){
         h[i]=x[i]+y[i];
     }

     return h;
     

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
  

return 100*norme_infinie(z,n)/norme_infinie(v,n);




}


double* uex2(int cas,double dx,double dy,int Nx,int Ny,int ibeg,int iend){
    int n = Nx * Ny;
    double *uex;
    uex = (double*)malloc((iend-ibeg+1) * sizeof(double));

for(int k=0;k<iend-ibeg+1;++k){
        int i=(k+ibeg)%Nx+1;
        int j=(k+ibeg)/Nx+1;
    
        if (cas==2){
            
               uex[k]= sin(i*dx)+cos(j*dy);
          
        }
        else if(cas==1){
               uex[k]=(i*dx)*(1-i*dx)*j*dy*(1-j*dy);

        }
}
  return uex;
  
   
}


double produit_scalaire(double *X, double *Y,int N)
{
    double s;
    s=0;
    for(int i=0;i<N;i++)
    {
        s=s+X[i]*Y[i];


    }
    return s;
}



double * multi_scal(double a,double *x,int n)
{
    double *h;
    h=(double*)malloc(n*sizeof(double));
    for(int i=0;i<n;i++)
    {
        h[i]=a*x[i];
    }
       
    return h;
 
   
}



//test
double *gradC_par(double alpha,double beta,double gamma,double *source,double eps,int kmax,int Nx,int Ny,int ibeg,int iend)
{
   
    double *x, *r, *z, *p, *r0;
    double a, re, b, y,sum,s,s1,sum1,s2,sum2;
    int k, n = Nx * Ny,err,rank,nproc;
    r = (double*)malloc((iend-ibeg+1) * sizeof(double));
    z = (double*)malloc((iend-ibeg+1) * sizeof(double));
    p = (double*)malloc((iend-ibeg+1) * sizeof(double));
    r0 = (double*)malloc((iend-ibeg+1) * sizeof(double));
    x = (double*)calloc((iend-ibeg+1),sizeof(double));
   err=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   err=MPI_Comm_size(MPI_COMM_WORLD,&nproc);


    memcpy(r, source, (iend - ibeg + 1) * sizeof(double)); 

    memcpy(p, r, (iend - ibeg + 1) * sizeof(double)); 
    sum=0;
    sum1=0;
    sum2=0;
    s=produit_scalaire(r,r,iend-ibeg+1);
    err=MPI_Allreduce(&s,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    re=sqrt(sum);
    
    
    k=0;
    while(k<kmax && re>eps)
    {
        s=produit_scalaire(r,r,iend-ibeg+1);
        err=MPI_Allreduce(&s,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        re=sqrt(sum);
         z=produitmat(alpha,beta,gamma,p,nproc,rank,ibeg,iend,Nx,Ny);
        s1=produit_scalaire(z, p,iend-ibeg+1);
        err=MPI_Allreduce(&s1,&sum1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        a=sum/sum1;   
         x = Sommevect(x, multi_scal(a, p, iend-ibeg+1), iend-ibeg+1);
         memcpy(r0, r, (iend - ibeg + 1) * sizeof(double)); 
       
        r = Sousvect(r0, multi_scal(a, z, iend-ibeg+1), iend-ibeg+1);
        free(r0);
        s2=produit_scalaire(r,r,iend-ibeg+1);
        err=MPI_Allreduce(&s2,&sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        b=sum2/sum;
        p = Sommevect(r, multi_scal(b, p, iend-ibeg+1), iend-ibeg+1);
        k = k + 1;       
        
    }



    free(r);
    free(z);
    free(p);
   
   

    return x;

   
}


double f(double x,double y,double t,double Lx,double Ly, int cas)
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


double g(double x,double y,double t,double Lx,double Ly, int cas)
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



double h(double x,double y,double t,double Lx,double Ly, int cas)
{
  if (cas == 1) {
    return 0;
  }

  else if (cas == 2) {
    return sin(x)+cos(y);
  }

  else if (cas == 3) {
    return 1;
  }
}



double* Source2(double D, double dt, double t, int Nx, int Ny, double* x, double* y, double* u, int cas, int ibeg, int iend) {
    double* b;
    b = (double*)malloc((iend - ibeg + 1) * sizeof(double));
    int i, j, I, k;
    double Lx = 1;
    double Ly = 1;
    double dx, dy;
    dx = Lx / (Nx + 1);
    dy = Ly / (Ny + 1);

    for (k = 0; k < iend - ibeg + 1; ++k) {
        i = (k + ibeg) % Nx + 1;
        j = (k + ibeg) / Nx + 1;

        if (j == 1) {
            b[k] = u[k] + dt * (f(x[i], y[j], t, Lx, Ly, cas) + D / (dy * dy) * g(x[i], 0, t, Lx, Ly, cas));

            if (i == 1)
                b[k] += dt * (D / (dx * dx) * h(0, dy, t, Lx, Ly, cas));
            if (i == Nx)
                b[k] += dt * (D / (dx * dx) * h(dx * (Nx + 1), dy, Lx, Ly, t, cas));
        }

        if (j >= 2 && j < Ny) {
            b[k] = u[k] + dt * f(i * dx, j * dy, t, Lx, Ly, cas);

            if (i == 1)
                b[k] += dt * (D / (dx * dx) * h(0, j * dy, t, Lx, Ly, cas));

            if (i == Nx)
                b[k] += dt * (D / (dx * dx) * h(dx * (Nx + 1), j * dy, t, Lx, Ly, cas));
        }

        if (j == Ny) {
            b[k] = u[k] + dt * (f(i * dx, j * dy, t, Lx, Ly, cas) + D / (dy * dy) * g(i * dx, (Ny + 1) * dy, t, Lx, Ly, cas));

            if (i == 1)
                b[k] += dt * (D / (dx * dx) * h(0, j * dy, t, Lx, Ly, cas));
            if (i == Nx)
                b[k] += dt * (D / (dx * dx) * h(dx * (Nx + 1), j * dy, t, Lx, Ly, cas));
        }
    }

    return b;
}


int main(int argc, char* argv[]) {

    int Nx, Ny, cas, nb_iterations, kmax, rank, ibeg, iend, err, nproc;
    double *v, *b, *w, *z, *x, *y, *h, *u, *solexacte, *vect;
    double epsilon;
    double alpha, beta, gama, tf, tn, dt;
    
    cas = 2;   // Paramètres à modifier
    Nx = 150;
    Ny = 150;

    double start_time = MPI_Wtime();
    epsilon = 0.00000001;
    kmax = 100000000;
    dt = 0.1;
    double Lx = 1, Ly = 1, D = 1;
    tf = 2;
    double dx, dy;
    
    dx = Lx / (Nx + 1);
    dy = Ly / (Ny + 1);
    int n = Nx * Ny, l;
    alpha = 1 + dt * D * ((2 / pow(dx, 2)) + 2 / pow(dy, 2));
    beta = (-1 / pow(dx, 2)) * D * dt;
    gama = (-1 / pow(dy, 2)) * D * dt;

    x = (double*)malloc((Nx + 2) * sizeof(double));
    y = (double*)malloc((Ny + 2) * sizeof(double));
    for (int i = 0; i < Nx + 2; i++) {
        x[i] = i * dx;
    }
    for (int i = 0; i < Ny + 2; i++) {
        y[i] = i * dy;
    }
    
    nb_iterations = ceil(tf / dt);
    dt = tf / nb_iterations;
    
    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    err = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    charge_c(rank, n, nproc, &ibeg, &iend);
    
    h = (double*)malloc((iend - ibeg + 1) * sizeof(double));
    u = (double*)malloc((iend - ibeg + 1) * sizeof(double));
    w = (double*)malloc((iend - ibeg + 1) * sizeof(double));
    
    for (int i = 0; i < iend - ibeg + 1; i++) {
        u[i] = 0; 
    }
    
    solexacte = uex2(cas, dx, dy, Nx, Ny, ibeg, iend);
    
    for (int ni = 0; ni < nb_iterations; ni++) {
        tn = (n + 1) * dt;
        h = Source2(D, dt, tn, Nx, Ny, x, y, u, cas, ibeg, iend);
        w = gradC_par(alpha, beta, gama, h, epsilon, kmax, Nx, Ny, ibeg, iend);
        memcpy(u, w, (iend - ibeg + 1) * sizeof(double)); 
        free(h);
        free(w);
        
        FILE* fichier; 
        char file_name[256]; 
        sprintf(file_name, "./sol.%d.%d.dat", ni, rank);
        fichier = fopen(file_name, "w");
        for (int I = 0; I < iend - ibeg + 1; I++) {
            int i = (I + ibeg) % Nx + 1;
            int j = (I + ibeg) / Nx + 1;
            
            fprintf(fichier, "%f %f %e \n", x[i], y[j], u[I], solexacte[I]);
        }
        
        fclose(fichier);
    }
    
    //printf("Sur le proc %d   L'erreur relative est de %e \n", rank, erreur(u, solexacte, iend - ibeg + 1));

    double end_time = MPI_Wtime();
    printf("Temps d'exécution pour le processus %d : %f secondes\n", rank, end_time - start_time);
    
    free(u);
    free(x);
    free(y);
    
    MPI_Finalize();
    
    return 0;
}
