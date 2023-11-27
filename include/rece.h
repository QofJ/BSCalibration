#define PI 3.14159265358979312
double getNKGdensity(double age,double size, double r){
    double cs,ro,rm=130.; //the molliere radius is fixed to be 130m
    cs = TMath::Gamma(4.5-age)/(TMath::Gamma(age-2.5+2.)*TMath::Gamma(4.5+2.5-2.0-2.*age));
    cs = cs/(2*PI*rm*rm);
    ro=cs*size*pow(r/rm,age-2.5)*pow(1+r/rm,age-4.5);
    return ro;
}
double recer50new3(double age,double size,double theta){
   double p2[12]={0.02503, 0.02361, 0.02489, 0.02357, 0.02149, 0.02303, 0.01852, 0.01409, 0.01289, 0.02556, 0.03406, 0.03635,};
   double p1[12]={0.96475, 0.94221, 0.91961, 0.90142, 0.88660, 0.86839, 0.86159, 0.85272, 0.84218, 0.82113, 0.79659, 0.77589};
   double p0[12]={1.74558, 1.75789, 1.77675, 1.80182, 1.83308, 1.86808, 1.91001, 1.95468, 1.99553, 2.04191, 2.09133, 2.14474};
   double r50,re;
   int k;
   if(age<0.6||age>2.4)return -1;
    if(theta<0||theta*57.3>51)return -1;
    k=int((1./cos(theta)-1)/0.05);
    if(k>11)k=11;
    if(k<0)k=0;
    r50=log10(getNKGdensity(age,size,50.));
    return pow(r50,2)*p2[k] + r50*p1[k] + p0[k];
}
int CutKM2AErecR50(double erec,double ratio){
 int i,n,m;
   double ll;
   double cut[12]={-5.11, -5.24, -5.95, -6.08, -2.34, -2.35, -2.36, -2.36, -2.36, -2.36,-2.36,-2.36};
   if(erec<1)return -1;
     n=int((erec-1)/0.2);
     if(n>11)n=11;
     if(ratio<cut[n])return n;
     else return -1;
}
double GetSpaceLC(double zen1, double phi1, double zen2, double phi2)
{
    if(zen1<0||zen2<0)return -1;

    double dr,dl0,dm0,dz0,dl1,dm1,dz1,cs;

       dl0= sin(zen1)*cos(phi1);
       dm0= sin(zen1)*sin(phi1);
       dz0= cos(zen1);

       dl1 = sin(zen2)*cos(phi2);
       dm1 = sin(zen2)*sin(phi2);
       dz1 = cos(zen2);
       cs=dl0*dl1+dm0*dm1+dz0*dz1;
       if(cs>1.0)cs=1.0;
     dr = acos(cs)*57.295780;
    return dr;
}
