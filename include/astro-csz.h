/********************************************
 *    Chen songzhan @ 2008                  *
 *    for ARGO-YBJ and LHAASO data analysis *
 *                                          *
 *                                          *
 *******************************************/
#ifndef _ASTROCSZ_
#define _ASTROCSZ_
//#include <iostream>
//#include <fstream.h>
//#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TNtupleD.h"
#include "TF1.h"
#include "TRandom.h"
#include "slalib.h"
#include "slamac.h"
#include "Sun.h"
#include "Moon.h"
#include "TFeldmanCousins.h"
#define hNDEG 10
#define hNRA 1400
#define hNDEC 1000
#define hPI 3.14159265358979312
#define YBJLA 0.512387992517690072 //29.357669444444447 degree daocheng
#define YBJLO 1.74775162355199343 //100.138791666666677 degree
#define  deg_rad 0.017453293      /*   pai/180  */
#define  rad_deg 57.29577951      /*   180/pai  */
#define NP 10


double  read_time(char *file){
   FILE *readfile;
   double time;
   int i;
   double dm;
   char str[10];
   readfile=fopen(file,"r");
   fscanf(readfile,"%s",str);
   fscanf(readfile,"\n");
   for(i=0;i<NP;i++){
       fscanf(readfile,"%lf\t",&dm);
       }
   fscanf(readfile,"\n%lf",&time);
   return time;
   fclose(readfile);

}



static float nTotalEff[NP],TimeOff,TimeOn;
int readeff(char *file,short int **nBin)
{
   int i,j,k,n,m,l;
   double total[NP];
   FILE *efffile;
   double dm;
      efffile=fopen(file,"r");
      if(!efffile){
              printf("There is no effective files!\n");	      
              exit(3);
      }
       printf("nTotalEff\n");
      for(i=0;i<NP;i++){
              fscanf(efffile,"%lf\t",&dm);nTotalEff[i]=dm; 
              printf("%.1lf\t",dm);
              total[i]=0;
      }
      printf("\n");
      fscanf(efffile,"\n%*d\t%*d\t%*d\n");
      fscanf(efffile,"%f\t%*f\t%*f\t%*f\n",&TimeOff);
      printf("%f\n",TimeOff);  
      {
            for(i=0;i<hNRA;i++){
                    for(j=0;j<hNDEC;j++){
                          for(k=0;k<NP;k++){
                               fscanf(efffile,"%d\t",&m);
                               total[k]+=m;
                               nBin[i*hNDEC+j][k]=m; 
                           }
			   fscanf(efffile,"\n");
                    }
            }
      }
      fclose(efffile);
      printf("nTotal\n");
      for(i=0;i<NP;i++){
              printf("%lf\t",total[i]);
              nTotalEff[i]=total[i];
       }
       printf("\n");
       return 0;
}
void init0(short int **nBin)
{
      int i,j,k;
      for(i=0;i<hNRA;i++){
          for(j=0;j<hNDEC;j++){
              for(k=0;k<NP;k++)  nBin[i*hNDEC+j][k]=0;
          }
      }
}
double Li_Ma(double aa,double ddon, double ddoff)
{
     double off,on,ss,all,aa1;
         all = ddon + ddoff;
         aa1 = aa + 1.;
	 if(ddon>0.5&&aa<0.2&&aa>0){
           if(fabs(ddon-ddoff*aa)>0.4){
            on  = ddon*log(aa1/aa*ddon/all);
            off = ddoff*log(aa1*ddoff/all);
            ss  = sqrt(2.)*sqrt(on+off);
	    if(ddon<aa*ddoff)ss=-ss;
           }
           else ss=(ddon-ddoff*aa)/sqrt(ddoff*aa+aa*ddon);
	 }
	 else ss=0;
         return ss;
}
//from J2000 to Jnow degree
double precession2New(double mjd,double *hra, double *ddec){
   double T=(mjd-51544.0)/365.24/100.;
   double M=1.2812323*T+0.0003879*T*T+0.0000101*T*T*T;
   double N=0.5567530*T-0.0001185*T*T-0.0000116*T*T*T;
   double ra,dec;
   ra = *hra+(M+N*sin(*hra*hPI/180.)*tan(*ddec*hPI/180.));
   dec=  *ddec+N*cos(*hra*hPI/180.);
   if(ra<0) ra+=360;
   if(ra>=360)ra-=360;
   *hra=ra;
   *ddec=dec;
   return 1;
}
//from Jnow to J2000 degree
double precession2Old(double mjd,double *hra, double *ddec){
   double T=(mjd-51544.0)/365.24/100.;
   double M=1.2812323*T+0.0003879*T*T+0.0000101*T*T*T;
   double N=0.5567530*T-0.0001185*T*T-0.0000116*T*T*T;
   double ra,dec;
   ra = *hra-(M+N*sin(*hra*hPI/180.)*tan(*ddec*hPI/180.));
   dec=  *ddec-N*cos(*hra*hPI/180.);
   if(ra<0) ra+=360;
   if(ra>=360)ra-=360;
   *hra=ra;
   *ddec=dec;
   return 1;
}
int src_position(const char *name, double mjd,double *ra1,double *dec1)
{
   double ra,dec;	
   if(strcmp(name,"MOON")==0)//(name=="MOON")
           moon_orbit(mjd, &ra, &dec);
   if(strcmp(name,"CRAB")==0){
	   ra=83.63;
           dec=22.02;
           //from J2000 to Jnow
           precession2New(mjd,&ra,&dec);
   }

   if(strcmp(name,"SUN")==0)
           sun_orbit(mjd, &ra, &dec);
   *ra1  = ra;
   *dec1 = dec;
   return 0;
}
//mjd to sidereal time in arc
double slaLst(double mjd, double LONGITUDE, double *lst )
{
  double t;
  t = 0.671262 + 1.002737909*(mjd - 40000)+LONGITUDE/D2PI;
  t = t - (int)t;
  *lst = D2PI * t;
  return *lst;
}


//Horizon to hour angle:  Az,El to HA,Dec
int h2eh_argo(double mjd, double zen, double az, double *hra, double *ddec)
{
    double z,a,tra,tdec;
    //a=360-(az+18.04+90.);  //ARGO-YBJ
   //a=360-(az+180.+18.04);   //KM2A@YBJ
    a=360-(az+180.-0.61);  //KM2A@Daoc
    if(a>=360) a = a-360;
    if(a<0)    a = a+360;
    z=90-zen;
    a=a/DR2D;
    z=z/DR2D;  
    slaDh2e(a, z, YBJLA,&tra, &tdec);
    *hra=-tra*DR2D;
    *ddec=tdec*DR2D;
    return 0;
}
//Horizon to equatorial coordinates:  Az,El to HA,Dec
int h2e_argo(double mjd,double zen, double az, double *hra, double *ddec)
{
    double z,a,tra,tdec,tside;
     h2eh_argo(mjd,zen,az,&tra,&tdec); 
     slaLst(mjd,YBJLO, &tside );
     tra=tside*DR2D+tra;
     if(tra<0)tra+=360.;
     if(tra>=360)tra-=360;  
    *hra=tra;
    *ddec=tdec;
    return 0;
}

//Horizon to equatorial coordinates in solar time:  Az,El to HA,Dec
int h2es_argo(double mjd,double zen, double az, double *hra, double *ddec)
{
    double z,a,tra,tdec,tsolar;
     h2eh_argo(mjd,zen,az,&tra,&tdec);
     tsolar=(mjd-(int)mjd)*360;
     tra=tsolar+tra;
     if(tra<0)tra+=360.;
     if(tra>=360)tra-=360;
    *hra=tra;
    *ddec=tdec;
    return 0;
}
// equatorial to Horizon  coordinates:  Az,El to HA,Dec
int e2h_argo(double mjd,double hra, double ddec,double *zen, double *az)
{
    double z,a,tra,tdec,tside;
     slaLst(mjd,YBJLO,&tside);
     tra=tside*DR2D - hra;
  //   if(tra<-180)tra+=360.;
  //   if(tra>=180)tra-=360;
      tra=tra/DR2D; tdec=ddec/DR2D;
     slaDe2h(tra,tdec,YBJLA,&a, &z);
     a=a*DR2D;
      a = 360. - a;
     // a = a - 108.04;
    // a = a - (180.+18.04);
      a = a - (180.-0.61);//DaoC
      if(a<0) a = a+360.;
      if(a>=360) a = a-360;
     *az=a;
     *zen=90-z*DR2D;
    return 0;
}

int e2g(double tra, double tdec, double *tlong, double *tlat)
{
  double ra,dec,l,b; 
     ra=tra*DR2D;
     dec=tdec*DR2D;
     slaEqgal(ra,dec,&l,&b);
     *tlong=l*DR2D;
     *tlat=b*DR2D;
     return 1;
}
int g2e(double tlong, double tlat, double *tra, double *tdec)
{
  double ra,dec,l,b;
     l=tlong*DR2D;
     b=tlat*DR2D;
     slaGaleq(l, b,&ra,&dec);
     *tra=ra*DR2D;
     *tdec=dec*DR2D;
     return 1;
}
//for panel draw
void slaAitoff(double ra,double dec,double *rra ,double *ddec)
{
 ra=ra*DPI/180.;
 dec=dec*DPI/180.;
 double x, y;
 double alpha2 = (ra-DPI)/2.;
 double delta = dec;
 double r2 = sqrt(2.);
 double f = 2*r2/DPI;
 double cdec = cos(delta);
 double denom = sqrt(1+cdec*cos(alpha2));
 x = cdec*sin(alpha2)*2*r2/denom;
 y = sin(delta)*r2/denom;
 x /= f;
 y /= f;
 x += DPI;
 *rra = x*180./DPI;
 *ddec = y*180./DPI;
}

double GetSpace_h(double decc1, double raa1, double decc2, double raa2)
{
    double dr,dl0,dm0,dz0,dl1,dm1,dz1,ra,dec; 
  
    ra=raa1*hPI/180.;
    dec=(90.-decc1)*hPI/180.;              
       dl0= sin(dec)*cos(ra);
       dm0= sin(dec)*sin(ra);
       dz0= cos(dec);
     ra=raa2*hPI/180.;
     dec=(90.-decc2)*hPI/180.;
       dl1 = sin(dec)*cos(ra);
       dm1 = sin(dec)*sin(ra);
       dz1 = cos(dec);
     dr = acos(dl0*dl1+dm0*dm1+dz0*dz1)*180./hPI;
    return dr;
}

double GetSpace(double decc1, double raa1, double decc2, double raa2)
{
    double dr,dl0,dm0,dz0,dl1,dm1,dz1,ra,dec;

    ra=raa1*hPI/180.;
    dec=(decc1)*hPI/180.;
       dl0= sin(dec)*cos(ra);
       dm0= sin(dec)*sin(ra);
       dz0= cos(dec);
     ra=raa2*hPI/180.;
     dec=(decc2)*hPI/180.;
       dl1 = sin(dec)*cos(ra);
       dm1 = sin(dec)*sin(ra);
       dz1 = cos(dec);
     dr = acos(dl0*dl1+dm0*dm1+dz0*dz1)*180./hPI;
    return dr;
}

double position( double RAS, double DEC ,double ras, double dec ){

  double sd, sD, sf, cd, cD, cf;
  double stheta, ctheta;

  sd = sin( dec * deg_rad );
  sD = sin( DEC * deg_rad );
  sf = sin( ( ras - RAS ) * deg_rad );
  cd = cos( dec * deg_rad );
  cD = cos( DEC * deg_rad );
  cf = cos( ( ras - RAS ) * deg_rad );

  stheta = cd * sf;
  ctheta = cD * sd - sD * cd * cf;
  return ( atan2( stheta, ctheta ) * rad_deg );

}

double FeldmanCousinsUpper(double Non, double Nb,double upper=0.9)//Nb<50 Non<50 libPhysics
{
   TFeldmanCousins f;
   f.SetCL(upper);
   return f.CalculateUpperLimit(Non, Nb);
}
double FeldmanCousinsLower(double Non, double Nb,double upper=0.68)//Nb<50 Non<50 libPhysics
{
   TFeldmanCousins f;
   f.SetCL(upper);
   return f.CalculateLowerLimit(Non, Nb);
}

#endif
