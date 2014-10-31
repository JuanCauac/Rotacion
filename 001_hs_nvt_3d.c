#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <time.h>
//### Variables
        int n;
        int imc, nmc, nsub, nsubmc, nmcold, nmctot;
        int ngr, ngrold, ngrtot;
        int acept;
        int ndivr;
        int ncx, ncy, ncz; 

        char newprog, orden, orientacion_orden;

        float x[1000], y[1000], z[1000];
        float xnew, ynew, znew;
        float dr, dx, dy, dz;
        float ax, ay, az;

        float deltar;
        float gmax;
        unsigned long g[5000];

        float lx, ly, lz;
        float lx2, ly2, lz2;
        float volumen;

        float rho;

        float sigma,sigma2,sigma_m;
        float dpmax;
        
        float A1[1000], A2[1000], A3[1000];
        
// ######### FUNCIONES
        void volumenes(void);
        void sigmas(void);
        void entrada(void);
        void verificacion(void);
        void initconf(void);
        void leeconf(void);
        void leenmc(void);
        void leerdf(void);
        void informacion(void);
        void mc(void);
        void condicionesperiodicas(double,double,double);
        void density(void);
        void rdf(void);
        void mathematica(void);
	// nuevas
	void rotacion(double,double,double,int);
	void orientacion_init(void);
	void orientacion_read(void);
	void orientacion_write(void);
	
        int baby(double,double,double,int);
        int traslape(double,double,double,int);

        double mindist(double,double,double, double,double,double);

//##################################

main()
{
        printf("     Empezamos\n");
        printf("\n");

        entrada();
        sigmas();
        volumenes();
        /*getchar();*/

/*      generacion semilla (seed)                   */
        struct timeval tv;
        gettimeofday(&tv, NULL);
        srand(tv.tv_sec * tv.tv_usec);

        switch(newprog){
        case 'n':
          printf("Empieza configuracion nueva\n");
          orientacion_init();
          initconf();
        
          nmcold=0;
          ngrold=0;
          printf("\n");
          /*getchar();*/
          break;
        case 'o':
          printf("Continua con la configuracion anterior\n");
          leenmc();
          leeconf();
          leerdf();
          orientacion_read();
          printf("\n");	
          /*getchar();*/
          break;
        default:
          break;
        }
        nmctot=nmc+nmcold;
        verificacion();
        /*getchar();*/

        informacion();
        mc();
        verificacion();
        rdf();
        mathematica();

        printf("terminamos\n");
} //Termina main

//##################################

void entrada()
{
        FILE *entra;
        entra = fopen("input","r");
        printf("Va a leer el archivo de entrada\n");
           fscanf(entra,"%d %d",&nmc,&nsub);
           fscanf(entra,"%d %f",&n,&rho);
           fscanf(entra,"%f",&sigma);
           fscanf(entra,"%f ",&dpmax);
           fscanf(entra,"%s %s ",&newprog,&orden);
           close(entra);
        printf("Termino de leer el archivo de entrada\n");
        printf("\n");
} //Termina entrada

//##################################

void sigmas()
{
        sigma2=sigma*sigma;
        sigma_m=sigma/2.0;
}// Termina sigmas

//##################################

void volumenes()
{
        volumen=n/rho;
        lx=pow(volumen,1.0/3.0);
        ly=lx;
        lz=lx;
        lx2=lx/2.0;
        ly2=ly/2.0;
        lz2=ly/2.0;

        deltar=0.05; 
        ndivr=floor((lx2)/deltar);

        FILE *caja;
        caja = fopen("caja","w");

        fprintf(caja,"%f %f %f \n",-lx2, -ly2, -lz2);
        fprintf(caja,"%f %f %f \n",-lx2,  ly2, -lz2);
        fprintf(caja,"%f %f %f \n", lx2,  ly2, -lz2);
        fprintf(caja,"%f %f %f \n", lx2, -ly2, -lz2);
        fprintf(caja,"%f %f %f \n",-lx2, -ly2, -lz2); //repetido
        fprintf(caja,"%f %f %f \n",-lx2, -ly2,  lz2);
        fprintf(caja,"%f %f %f \n",-lx2,  ly2,  lz2);
        fprintf(caja,"%f %f %f \n", lx2,  ly2,  lz2);
        fprintf(caja,"%f %f %f \n", lx2, -ly2,  lz2);
        fprintf(caja,"%f %f %f \n",-lx2, -ly2,  lz2); //
        close(caja);
}// Termina volumenes

//##################################

void initconf()
{
        printf("Inicia la configuracion de las particulas\n");
        int i,j,k,l,born;
        double xt,yt,zt;
        double xo1,yo1,zo1,xo2,yo2,zo2;
        double xo3,yo3,zo3,xo4,yo4,zo4;

        FILE *fp;
        fp = fopen("initconf", "w");

        FILE *fmi;
        fmi = fopen("initconf.nb", "w");
        fprintf(fmi,"Graphics3D[{");

        i=1;
        born=1;

        switch(orden){
        case 'd':
          printf("Se utiliza una configuracion inicial desordenada\n");
          while (i<=n) {

            xt = (2.0*(double)rand()/(double)RAND_MAX-1.0)*lx2;
            yt = (2.0*(double)rand()/(double)RAND_MAX-1.0)*ly2;
            zt = (2.0*(double)rand()/(double)RAND_MAX-1.0)*lz2;

/*          printf("coordenadas tentativas:");
            printf("x( %d )= %f ",i,xt);
            printf("y( %d )= %f ",i,yt);
            printf("z( %d )= %f \n",i,zt);*/

            born=baby(xt,yt,zt,i);

            if (born==1) {
               x[i]=xt;
               y[i]=yt;
               z[i]=zt;
               fprintf(fp,"%f \t %f \t %f \t %d \n",x[i],y[i],z[i],i);
               if (i == n){
                  fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
                  fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]}",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);}

               else if (i != n){
                  fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
                  fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]},",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);}
               printf("nace particula %d\n",i);
               i++;
            }  
          }
          /*getchar();*/
          break;
        case 'o':
          printf("Se utiliza una configuracion inicial ordenada\n");
          printf("el numero de particulas puede ser 4, 32, 108, 256, 500, 804 o 1372 \n");
          printf("si el numero de particulas no es alguno de estos puede haber errores en el programa \n");
          switch(n){
          case 1:
            ncx=1;
            break;
          case 32:
            ncx=2;
            break;
          case 108:
            ncx=3;
            break;
          case 256:
            ncx=4;
            break;
          case 500:
            ncx=5;
            break;
          case 804:
            ncx=6;
            break;
          case 1372:
            ncx=7;
            break;
          default:
            ncx=(int)pow((double)n/4.0,1.0/3.0);
            break;
          }
          ncy=ncx;
          ncz=ncx;

          ax=lx/ncx;
          ay=ly/ncy;
          az=lz/ncz;

          xo1=-lx2+ax/4.0;
          yo1=-ly2+ay/4.0;
          zo1=-lz2+az/4.0;

          xo2=-lx2+3.0*ax/4.0;
          yo2=-ly2+3.0*ax/4.0;
          zo2=-lz2+az/4.0;

          xo3=-lx2+3.0*ax/4.0;
          yo3=-ly2+ay/4.0;
          zo3=-lz2+3.0*az/4.0;

          xo4=-lx2+ax/4.0;
          yo4=-ly2+3.0*ax/4.0;
          zo4=-lz2+3.0*az/4.0;

          printf("origen 1  %f %f %f \n",xo1,yo1,zo1);
          printf("origen 2  %f %f %f \n",xo2,yo2,zo2);
          printf("origen 3  %f %f %f \n",xo3,yo3,zo3);
          printf("origen 4  %f %f %f \n",xo4,yo4,zo4);
          /*getchar();*/

          while (i<=n) {
            for (k=1;k<=ncx;k++){
              for (j=1;j<=ncy;j++){
                for (l=1;l<=ncz;l++){
                  x[i]=xo1 + (j-1)*ax;
                  y[i]=yo1 + (k-1)*ay;
                  z[i]=zo1 + (l-1)*az;
                  printf("nace particula %d\n",i);
                  fprintf(fp,"%f\t%f\t%f\t%d\n",x[i],y[i],z[i],i);
                  fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
                  fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]},",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);
                  i++;
                  x[i]=xo2 + (j-1)*ax;
                  y[i]=yo2 + (k-1)*ay;
                  z[i]=zo2 + (l-1)*az;
                  printf("nace particula %d\n",i);
                  fprintf(fp,"%f\t%f\t%f\t%d\n",x[i],y[i],z[i],i);
                  fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
	          fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]},",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);
                  i++;
                  x[i]=xo3 + (j-1)*ax;
                  y[i]=yo3 + (k-1)*ay;
                  z[i]=zo3 + (l-1)*az;
                  printf("nace particula %d\n",i);
                  fprintf(fp,"%f\t%f\t%f\t%d\n",x[i],y[i],z[i],i);
                  fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
                  fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]},",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);
                  i++;
                  x[i]=xo4 + (j-1)*ax;
                  y[i]=yo4 + (k-1)*ay;
                  z[i]=zo4 + (l-1)*az;
                  printf("nace particula %d\n",i);
                  fprintf(fp,"%f\t%f\t%f\t%d\n",x[i],y[i],z[i],i);
                  if (i != n){
                     fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
                     fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]},",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);}
                  else if (i == n){
                     fprintf(fmi,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
                     fprintf(fmi,"{White, Sphere[{%f, %f, %f}, %f]}",x[i]+A1[i]*sigma_m,y[i]+A2[i]*sigma_m,z[i]+A3[i]*sigma_m,sigma_m/3);}
                  i++;
                }
              }
            }
          }
          break;
        default:
          break;
        }

        /*fprintf(fmi,"Black,Line[{{-%f,%f},{%f,%f},{%f,-%f},{-%f,-%f},{-%f,%f}}]}]",lx2,ly2,lx2,ly2,lx2,ly2,lx2,ly2,lx2,ly2);*/
        fprintf(fmi,"}]");

        fclose(fp);
        fclose(fmi);
        printf("Termina la configuracion de las particulas\n");
} // Termina initconf

//##################################

int baby(double xt, double yt, double zt, int i)
{
        int j,u;

/*        printf("subrutina baby \n");*/
        u=1;

        if(i>1){
           for(j=1;j<=i-1;j++) {
              dr=mindist(xt,yt,zt,x[j],y[j],z[j]);

              if (dr<=sigma2) {
                 u=0;
/*                 printf("traslape de particula %d ",i);
                 printf(" con particula %d\n",j);
                 printf("x( %d",i);
                 printf(" )= %f",xt);
                 printf(" y( %d ",i);
                 printf(" )= %f \n",yt);
                 printf(" z( %d ",i);
                 printf(" )= %f \n",zt);
                 printf("x( %d ",j);
                 printf(" )= %f ",x[j]);
                 printf(" y( %d ",j);
                 printf(" )= %f \n",y[j]);
                 printf(" z( %d ",j);
                 printf(" )= %f \n",z[j]);
                 printf(" dr= %f \n",dr);
                 printf(" sigma2= %f \n",sigma2);
                 getchar();*/
                 break;
                 return u;                                 
              } 
              else {
                 u=1;
/*                 printf(" no hay traslape de particula %d ",i);
                 printf(" con particula %d \n",j);
                 printf(" dr= %f \n",dr);
                 printf(" sigma2= %f \n",sigma2);
                 getchar();*/
              }
           }
        }
        return u;
} //Termina baby

//##################################

void informacion()
{
        printf("Informacion general de la simulacion");
        printf("\n");
        printf("Numero de particulas %d \n",n);
        printf("Un sitio de asosiación en (1,0,0) \n"); // <- modificado
        printf("densidad %f \n",rho);
        printf("Longitudes de la caja de simulacion\n");
        printf("lx = %f\n",lx);
        printf("ly = %f\n",ly);
        printf("lz = %f\n",lz);
        printf("Volumen de la caja de simulacion\n");
        printf("volumen=  %f\n",volumen);
        if (orden=='o'){
           printf("Se inicio con una configuracion ordenada\n");
           printf("Numero de celdas en la caja de simulacion\n");
           printf("numero de celdas en x %d\n",ncx);
           printf("numero de celdas en y %d\n",ncy);
           printf("numero de celdas en z %d\n",ncz);
           printf("parametro de red en x %f\n",ax);
           printf("parametro de red en y %f\n",ay);
           printf("parametro de red en z %f\n",az);
           printf("Hay cuatro particulas por celda\n");
        }
        printf("Numero de divisiones para calcular la funcion de distribucion radial\n");
        printf("ndivr = %d\n",ndivr);
        printf("\n");
        printf("El numero de configuraciones Monte Carlo anteriores es de %d\n",nmcold);
        printf("El numero de configuraciones de este programa son %d\n",nmc);
} //Termina informacion

//##################################

void leenmc()
{
        FILE *fnmc;
        fnmc = fopen("nmc", "r");

        fscanf(fnmc,"%d %d",&nmcold,&ngrold);
     
        fclose(fnmc);
} // Termina leenmc

//##################################

void leerdf()
{
        int j;
        float trash;

        FILE *frdf;
        frdf = fopen("rdf", "r");

        for (j=0; j<ndivr; j++)
            fscanf(frdf,"%f %lu %f",&trash,&g[j],&trash);
     
        fclose(frdf);
} // Termina leerdf

//##################################

void leeconf()
{
        int j;

        FILE *f;
        f = fopen("finalconf", "r");

        for (j=1;j<=n;j++) {
        
        fscanf(f,"%f %f %f",&x[j],&y[j],&z[j]);
        /*printf("x y z %f %f %f\n",x[j],y[j],z[j]);
        getchar();*/
        }
        close(f);
} // Termina leeconf

//##################################

void verificacion()
{
        int k,j;

        for(k=1;k<n-1;k++) {
          for(j=k+1;j<=n;j++){
            dr=mindist(x[k],y[k],z[k],x[j],y[j],z[j]);

            if(dr<=sigma2) {
              printf("Error en verificacion:\n");
              printf("traslape de particula %d ",k);
              printf("con particula %d \n",j);
              printf("x( %d",k);
              printf(" )= %f",x[k]);
              printf(" y( %d ",k);
              printf(" )= %f ,",y[k]);
              printf(" z( %d ",k);
              printf(" )= %f \n",z[k]);
              printf("x( %d ",j);
              printf(" )= %f ",x[j]);
              printf(" y( %d ",j);
              printf(" )= %f ,",y[j]);
              printf(" z( %d ",j);
              printf(" )= %f \n",z[j]);
              printf(" dr= %f \n",dr);
              printf(" sigma2= %f \n",sigma2);
            }
          }
        }
}// Termina verificacion

//##################################

void mc()
{
        printf("\n Inicia Monte Carlo\n");
        printf("\n");

        int i,j,ir,in,mov;

        double xold,yold,zold;
        double dr;

        FILE *ff;
        ff = fopen("finalconf", "w");

        FILE *fnmc;
        fnmc = fopen("nmc", "w");

        acept=0;
        nsubmc=0;
        ngr=0;

        for (imc=1;imc<=nmc;imc++){
           for (in=1;in<=n;in++){
              i = ((double)rand()/(double)RAND_MAX)*n + 1;
/*             printf("Intento de mover particula %d\n",i);
              printf("\n");*/
              mov=1;
/*              printf("mov %d \n",mov);*/
              xold=x[i];
              yold=y[i];
              zold=z[i];
/*             printf("coordenadas viejas %f %f\n",xold,yold,zold);
              printf("\n");*/
              xnew = xold + (2.0*(double)rand()/(double)RAND_MAX-1.0)*dpmax;
              ynew = yold + (2.0*(double)rand()/(double)RAND_MAX-1.0)*dpmax;
              znew = zold + (2.0*(double)rand()/(double)RAND_MAX-1.0)*dpmax;
/*              printf("coordenadas nuevas %f %f\n",xnew,ynew,znew);
              printf("\n");*/
              condicionesperiodicas(xnew,ynew,znew);
/*              printf("condiciones periodicas %f %f \n",xnew,ynew,znew);
              printf("\n");*/
              mov=traslape(xnew,ynew,znew,i);
              if (mov==1){
                 x[i]=xnew;
                 y[i]=ynew;
                 z[i]=znew;
		rotacion(A1[i],  A2[i], A3[i], i);
                 acept++;
                 if (abs(x[i]) > lx2 || abs(y[i])>ly2 || abs(z[i])>lz2){
                    printf("hubo aceptacion %d \n",mov);
                    printf("x = %f\n",x[i]); 
                    printf("y = %f\n",y[i]); 
                    printf("z = %f\n",z[i]); 
                    printf("acept = %d\n",acept);
                    printf("coordenadas viejas\n");
                    printf("xold = %f\n",xold);                   
                    printf("yold = %f\n",yold);                   
                    printf("zold = %f\n",zold);                   
                    printf("numero mc = %d\n",imc);
                    getchar();
                 }
              }
              for (j=1;j<=n;j++) {
                 if (j!=i){
                    dr=mindist(x[j],y[j],z[j],x[i],y[i],z[i]);
                    dr=sqrt(dr);
                    if (dr<=lx2){
                       ir=floor(dr/deltar);
                       g[ir]=g[ir]+1;
                    }
                 }
              }
              ngr++;
           }
           if ((imc % nsub) == 0){
              nsubmc++;
              printf("nmc %d, nsubrutina %d, numero de rdf %d \n",imc,nsubmc,ngr);
              printf("aceptacion %f\n",(double)acept/(imc*n));
           }
        }
        for (j=1;j<=n;j++) {
            fprintf(ff,"%f\t%f\t%f\n",x[j],y[j],z[j]);
        }
	orientacion_write();
	
        ngrtot=ngr+ngrold;
        fprintf(fnmc,"%d %d\n",nmctot,ngrtot);

        close(ff);
        close(fnmc);
        printf("Termina subrutina Monte Carlo \n");
}// termina mc

//##################################

void condicionesperiodicas(double xx, double yy, double zz) // mejorar
{
        int signox,signoy,signoz;

        if (xx > lx2) signox = -1;
        if (xx < lx2) signox =  1;
        if (yy > ly2) signoy = -1;
        if (yy < ly2) signoy =  1;
        if (zz > lz2) signoz = -1;
        if (zz < lz2) signoz =  1;
        if((xx > -lx2) && (xx < lx2)) signox = 0;
        if((yy > -ly2) && (yy < ly2)) signoy = 0;
        if((zz > -lz2) && (zz < lz2)) signoz = 0;

        xnew = xx + signox*lx;
        ynew = yy + signoy*ly;
        znew = zz + signoz*lz;
} // termina condiciones periodicas

//##################################

double mindist(double x1, double y1, double z1, double x2, double y2, double z2)
{
        dx = fabs(x1 - x2);
        if (dx > lx2){
           dx = dx - lx;
        }

        dy = fabs(y1 - y2);
        if (dy > ly2){
             dy = dy - ly;
        }

        dz = fabs(z1 - z2);
        if (dz > lz2){
             dz = dz - lz;
        }

/*         printf("dx= &f dy= %f dz= %f dr= %f \n",dx,dy,dx*dx+dy*dy+dz*dz);
        getchar();*/
        return dx*dx + dy*dy + dz*dz;

} // Termina mindist

//##################################

int traslape(double xt, double yt, double zt, int i)
{
        int j,u;

/*      printf("subrutina traslape\n");                  */
        u=1;

        for (j = 1; j <= n;j++){
           if (j != i){
              dr = mindist(xt, yt, zt, x[j], y[j], z[j]);
              if (dr<=sigma2) {
                 u=0;
/*                 printf("traslape de %d ",i);
                 printf(" con %d\n",j);                   */
                 break;
              }
              else {
              u=1;
/*              printf(" no hay traslape de %d ",i); 
              printf(" con %d \n",j);                     */
              }
           }
        }
        return u;
} // Termina traslape

void rdf(){

        int j;
        float r,volcapa,deltar2pi,gg;

        FILE *fgs;
        fgs = fopen("rdf", "w");

        FILE *fg;
        fg = fopen("rdf.nb", "w");
        fprintf(fg,"ListLinePlot[{");

        deltar2pi=4.0*deltar*deltar*deltar*M_PI*rho*ngrtot/3.0;
        gmax=0.0;

        for (j=0;j<=ndivr-1;j++){
            r=((float)j+0.5)*deltar;
            volcapa=((j+1)*(j+1)*(j+1)-(j)*(j)*(j))*deltar2pi;
            gg=(double)g[j]/volcapa;
            if (j!=ndivr-1)
               fprintf(fg,"{%f , %f },",r,gg);
            else 
               fprintf(fg,"{%f , %f }",r,gg);
            fprintf(fgs,"%f %ld %f\n ",r,g[j],gg);
            if (gg > gmax)
               gmax = gg;
        }
        fprintf(fg,"},PlotRange->{{0, %f },{0, %f }}]",lx2,gmax);
        fclose(fg);
} // Termina rdf

//##################################

void mathematica()
{
        int i;

        FILE *fm;
        fm = fopen("finalconf.nb", "w");
        fprintf(fm,"Graphics3D[{");

	double 	s=sigma_m;

        for (i=1; i<n; i++){ 
            fprintf(fm,"{Red, Sphere[{%f, %f, %f}, %f]},",x[i],y[i],z[i],sigma_m);
            fprintf(fm,"{White, Sphere[{%f, %f, %f}, %f]},",x[i]+A1[i]*s,y[i]+A2[i]*s,z[i]+A3[i]*s,s/3);
        }
            fprintf(fm,"{Red, Sphere[{%f, %f, %f}, %f]},",x[n],y[n],z[n],sigma_m);
            fprintf(fm,"{White, Sphere[{%f, %f, %f}, %f]}",x[n]+A1[n]*s,y[n]+A2[n]*s,z[n]+A3[n]*s,s/3);

        fprintf(fm,"}]");
        fclose(fm);
}// Termina mathematica


void orientacion_init(){
// dota de orientacion inicial a las particulas
	int i;
	if (orientacion_orden=='o'){
		for(i=1; i<=n; i++){
			A1[i]=1;
			A2[i]=0;
			A3[i]=0;
		}//for
	}//if
	else{
		for(i=1; i<=n; i++){
			double ty= (double)rand()/RAND_MAX*M_PI;
			double tz= (double)rand()/RAND_MAX*M_PI;
			A1[i]=cos(ty)*cos(tz);
			A2[i]=-cos(ty)*sin(tz);
			A3[i]=sin(ty);
		}//for
	}//else



}//termina orientacion


void rotacion(double a1, double a2, double a3, int i){
	double tx=(double)rand()/RAND_MAX*M_PI/2;///180; // a lo mucho 1 grado
	double ty=(double)rand()/RAND_MAX*M_PI/2;///180; // se puede optimizar  RAND_MAX*M_PI/180 -> 0.0174532925
	double tz=(double)rand()/RAND_MAX*M_PI/2;///180;
	
	double Cx=cos(tx);
	double Sx=sin(tx);
	
	double Cy=cos(ty);
	double Sy=sin(ty);
	
	double Cz=cos(tz);
	double Sz=sin(tz);
	
	double a1New=a1*(Cy*Cz)+a2*(Sx*Sy*Cz+Cx*Sz)+a3*(-Cz*Sy*Cz+Sx*Sz);
	double a2New=a1*(-Cy*Sz)+a2*(Cx*Cz-Sx*Sy*Sz)+a3*(Sx*Cz+Cx*Sy*Sz);
	double a3New=a1*(Sy)+a2*(-Cy*Sx)+a3*(Cx*Cy);
	
	double norm=sqrt(a1New*a1New+a2New*a2New+a3New*a3New);
	
	A1[i]=a1New/norm;
	A2[i]=a2New/norm;
	A3[i]=a3New/norm;

}//termina rotacion

void orientacion_write(){
	int j;
 	FILE *f;
        f = fopen("orientacion","w");
         for (j=1;j<=n;j++) {
            fprintf(f,"%f %f %f \n",A1[j], A2[j], A3[j]);
	close(f);
        }
        
}//termina orientacion_write

void orientacion_read(){
 	int j;
        FILE *f;
        f = fopen("orientacion", "r");
        for (j=1;j<=n;j++) {
        fscanf(f,"%f %f %f",&A1[j],&A2[j],&A3[j]);
        }
        close(f);
}//termina orientacion_read
