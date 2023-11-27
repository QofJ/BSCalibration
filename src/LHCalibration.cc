// -----------------------------------------------------------------------
//                                LHAASO
//                       --- LHCalibration ---
//                        class implementation file
// -----------------------------------------------------------------------
//  Description:
//
//    Managers of the LHAASO detector calibration and status
//
// -----------------------------------------------------------------------
// History:
//
//    chen songzhan      -- 2021.4 created
//           
//========================================================================
#include "LHCalibration.h"
#include <iostream>
#include <fstream>
#include <algorithm>   
#include <iterator>
#include "math.h"
#include <string.h>
#include "EOSopen.h"
using namespace std;

LHCalibration * LHCalibration::m_myself = 0;
G4KM2A_Geometry * LHCalibration::geom=0;
LHCalibration::LHCalibration()
{
    geom = G4KM2A_Geometry::GetInstance(7);
    int i;
    fNED=0;
    fNMD=0;
    fNWCDA=0;
    fNWCDA2=0;
    for(i=0;i<NNED;i++){
         fEDdt.push_back(0);
         fEDda.push_back(1);
         fEDratio.push_back(56);
         fEDstatus.push_back(0);
         fEDcount.push_back(0);
    }
    for(i=0;i<NNMD;i++){
         fMDdt.push_back(120);
         fMDda.push_back(1);
         fMDratio.push_back(103);
         fMDstatus.push_back(0);
         fMDcount.push_back(0);
    }
    for(i=0;i<NNWCDA;i++){
         fWCDAdt.push_back(0);
         fWCDAda.push_back(1);
         fWCDAratio.push_back(50);
         fWCDAstatus.push_back(0);
    }
    for(i=0;i<NNWCDA;i++){
         fWCDA2dt.push_back(0);
         fWCDA2da.push_back(1);
         fWCDA2ratio.push_back(50);
         fWCDA2status.push_back(0);
    }
}
LHCalibration::~LHCalibration()
{

}
int LHCalibration::SetTDCcalibrationKM2A(const char filename[100] )
{
    int i,j,n,m,id;
    float dtt;
    char name[20];
    FILE *calfile;
    for(i=0;i<NNED;i++){  fEDdt[i]=-1000;   }
    for(i=0;i<NNMD;i++){  fMDdt[i]=-1000;   }
    n=0;m=0;
    calfile=fopen(filename,"r");
    if(!calfile){printf("%s not exist or can't read! procedure end up!\n",filename); return 1;}
    KM2ATDCfilename=filename;
    printf("using TDC calibration file %s\n",filename);
    while(!feof(calfile)){
         fscanf(calfile, "%s %d %f\n", name,&id, &dtt);
         if(strcmp(name,"ED")==0){ fEDdt[id]=dtt; n++;}
         else if (strcmp(name,"MD")==0){ fMDdt[id]=dtt; m++; }
    }
    fNED=n;
    fNMD=m;
    printf("TDC %d ED, %d MD\n",n,m);
    fclose(calfile);
 //   SetFinalStatusKM2A();
    return 0; 
}
int LHCalibration::SetADCcalibrationKM2A(const char filename[100] ) 
{
    int i,j,n,m,id;
    float daa,drr;
    char name[20];
    FILE *calfile;
    for(i=0;i<NNED;i++){  fEDda[i]=0;  fEDratio[i]=0;  }
    for(i=0;i<NNMD;i++){  fMDda[i]=0;  fMDratio[i]=0;  }
    n=0;m=0;
    calfile=fopen(filename,"r");
    if(!calfile){printf("%s not exist or can't read! procedure end up!\n",filename); return 1;}
    KM2AADCfilename=filename;
    printf("using ADC calibration file %s\n",filename);
    while(!feof(calfile)){
        fscanf(calfile, "%s %d %f %*f %f %*f\n",name, &id, &daa,&drr);
        if(strcmp(name,"ED")==0){
            if(daa>10&&daa<50){
                fEDda[id]=28.24/daa;
                if(drr>10&&drr<150)fEDratio[id]=drr;
                else fEDratio[id]=0;
                n++;
            }
            else{
                  fEDda[id]=0;
                  fEDratio[id]=0;
            }
        }
        else if(strcmp(name,"MD")==0){
            if(daa>15&&daa<150){
                fMDda[id]=77.66/daa;
                if(drr>10&&drr<150)fMDratio[id]=drr;
                else fMDratio[id]=0;
                  m++;
            }   
            else{
                fMDda[id]=0;
                fMDratio[id]=0;
            }   
              
        }
     }
     fclose(calfile); 
     fNED=n;
     fNMD=m; 
     printf("ADC calibration set: %d ED, %d MD\n",n,m);
 //    SetFinalStatusKM2A();
     return 0;
}
int LHCalibration::SetMDcalibration(const char filename[100] )
{
    int i,j,n,m,id;
    float ta,tt,tr;
    FILE *calfile;
    float MDdt2[NNMD];
    if(fNMD<10){
       for(i=0;i<NNMD;i++){  
          fMDdt[i]=-1000;
          fMDda[i]=0; 
          fMDratio[i]=0;
       }
    }
    for(i=0;i<NNMD;i++) MDdt2[i]=-1000;

    n=0;
    calfile=fopen(filename,"r");
    if(!calfile){printf("%s not exist or can't read! procedure end up!\n",filename); return 1;}
    MDfilename=filename;
    printf("using MD calibration file %s\n",filename);
    fscanf(calfile,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n");
    while(!feof(calfile)){
         fscanf(calfile,"%*s %d %f %*lf %*lf %*lf %f %*lf %f %*lf %*lf\n",&id,&ta,&tt,&tr);
         //fMDda[id]=ta;
         MDdt2[id]=tt; 
         //fMDratio[id]=tr;
         n++;
    }
    fNMD=n;
    printf("Set %d MD TDC\n",n);
    fclose(calfile);
    //reset MD dt
    n=0;
    double all=0.;
    for(i=0;i<NNMD;i++){
          if(MDdt2[i]>-900&&fMDdt[i]>-900){
               n++;
               all += MDdt2[i]-fMDdt[i];
          }      
    }
    if(n>100) all = all/n;
    else all=0.;
    printf("Reset MD TDC %.2lf ns\n",all); 
    for(i=0;i<NNMD;i++){
         if(MDdt2[i]>-900){
               fMDdt[i]=MDdt2[i]-all;
          }      
    }
 //   SetFinalStatusKM2A();
    return 0;
}
int LHCalibration::SetStatusKM2A(const char filename[100],int cflag=0){
    KM2AStatusfilename = filename;
    TFile *infile=EOSopen(filename,0);
    if(!infile){
        printf("open file error %s !!\n",filename);
        infile->Close();
        return 1;
    }
    TNtupleD *nt = (TNtupleD *)infile->Get("nt1");
    if(!nt){
        printf("get nt1 error %s !!\n",filename);
        infile->Close();
        return 1;
    }
    int i,n=0,m=0;
    for(i=0;i<NNED;i++){  fEDstatus[i]=-1; }
    for(i=0;i<NNMD;i++){  fMDstatus[i]=-1; }
    //mode:id:Frate:Foccu:Foccu2:rate:occu:dtpeak:dtmed:qamax:n:mjd:nevent
    double fmode,fid,fFrate,fFoccu,fFoccu2,dtmed,nevent;
    int mode,id,Frate,Foccu,Foccu2;
    nt->SetBranchAddress("mode",&fmode);
    nt->SetBranchAddress("id",&fid);
    nt->SetBranchAddress("Frate",&fFrate);
    nt->SetBranchAddress("Foccu",&fFoccu);
    nt->SetBranchAddress("Foccu2",&fFoccu2);
    nt->SetBranchAddress("dtmed",&dtmed);
    nt->SetBranchAddress("nevent",&nevent);
    int nTotalentry = nt->GetEntriesFast();
    for(i=0;i<nTotalentry;i++){
        nt->GetEntry(i);
        mode  =int(fmode);
        id    =int(fid);
        Frate =int(fFrate);
        Foccu =int(fFoccu);
        Foccu2=int(fFoccu2);
        if(mode==0){
             if(Frate==0&&Foccu==0){ 
                fEDstatus[id]=0;
//                if(fabs(dtmed)>5)printf("ED dt: %d %lf\n",id,dtmed);
                n++;  
             }
        }
        else if(mode==1){
             if(Frate==0&&Foccu==0){
                fMDstatus[id]=0; 
                if(cflag==1&&fMDdt[id]>-900){
                   fMDdt[id]=fMDdt[id]+(dtmed-1);  
                }
//                if(fabs(dtmed)>10)printf("MD dt: %d %lf\n",id,dtmed);
                m++;
             }
        }        
    }
    infile->Close();
    fNED=n;
    fNMD=m;
    printf("KM2A Status set: %d good ED, %d good MD\n",n,m);
    //SetFinalStatusKM2A();
    return 0;
}
int LHCalibration::SetFinalStatusKM2A(){
   int i,j,n=0,m=0;
   double x,y,z;
   for(i=0;i<NNED;i++){ 
      if(fEDstatus[i]==0&&(fEDdt[i]<-900||fEDda[i]==0))fEDstatus[i] = 1;
      if(geom->Getxyz(i,x,y,z,1,"ED")<0)fEDstatus[i] = 2;
      if(fEDstatus[i]==0)n++;
   }
   for(i=0;i<NNMD;i++){  
      if(fMDstatus[i]==0&&(fMDdt[i]<-900||fMDda[i]==0))fMDstatus[i] = 1;
      if(geom->Getxyz(i,x,y,z,1,"MD")<0)fMDstatus[i] = 2;  
      if(fMDstatus[i]==0)m++;
   }
   fNED = n;
   fNMD = m;
   printf("KM2A Final Status set: %d good ED, %d good MD\n",n,m);
   return 0;
}
int LHCalibration::ApplyKM2AEvent(KM2AEvent *event,LHEvent *lhevent){
    int i,j,k,id,idt,tag,mode,nflag;
    double x,y,z,t,dt,mjd,pe,pd,peda,pedd;
    TClonesArray *Hits;
    KM2AHit *Hit;
    lhevent->Initcsz();
    Hits = event->Hitlist();
    k = event->Nmd()+event->NHit();
    nflag=0;//to check half event
    for(j=0; j<k; j++){
         Hit = (KM2AHit *)((*Hits)[j]);
         id=Hit->Id();
         t=Hit->Time();
         if(t<4000)nflag++;
         mode=Hit->Mode();

         pe=Hit->Qa();
         pd=Hit->Qd();
         peda=Hit->Peda();
         pedd=Hit->Pedd();
         tag=Hit->Tag();
         if(mode==0){
             //if(geom->Getxyz(id,x,y,z,1,"ED")<0)continue;
             //if(fEDda[id]==0)continue;
             fEDcount[id]++;
             if(fEDstatus[id]!=0)continue;
             if(pe>0){
                 pe=pe-peda/8.;
                 pd=pd-pedd/8.;
                 if(pe>3700&&fEDratio[id]>10&&fEDratio[id]<150) pe=pd*fEDratio[id];
                 pe = pe*fEDda[id];
             }
             t=t-fEDdt[id];

             if(fEDdt[id]>-900&&fEDda[id]!=0&&pe>=0){
                 lhevent->AddHitE(id,t,pe,0);
             }
          }
          else if(mode==1){
              //if(geom->Getxyz(id,x,y,z,1,"MD")<0)continue;
              //if(fMDda[id]==0)continue;
              fMDcount[id]++;
              if(fMDstatus[id]!=0)continue;
              if(pe>0){
                   pe=(pe-peda/8.*304)*0.0319;
                   pd=(pd-pedd/8.*304)*0.0319;
                   if(tag>0&&fMDratio[id]>10&&fMDratio[id]<150) pe=pd*fMDratio[id];
                   pe = pe*fMDda[id];
              }
              t=t-fMDdt[id];
              if(fMDdt[id]>-900&&fMDda[id]!=0&&pe>=0)  {
                  idt=id;
                  if(event->Mjd()<58978.4327986){ //temporay for MD DAQ id error
                    if(id==450) idt=1236; 
                    if(id==1028)idt=1058;
                    if(id==1058)idt=1028;
                  }
                  if(event->Mjd()<59418.5202){//error corrected on 2021-7-23 run117916,
                     if(id==78)idt=103;
                     if(id==103)idt=78;
                  } 
                  lhevent->AddHitM(idt,t,pe,0);
            }
       }
    }
    return nflag;
}
int LHCalibration::GetEDcount(){
   int i,n=0;
   for(i=0;i<NNED;i++){
      if(fEDcount[i]>0)n++;
   }
   return n;
}
int LHCalibration::GetEDcountGood(){
   int i,n=0;
   for(i=0;i<NNED;i++){
      if(fEDcount[i]>0&&fEDstatus[i]==0)n++;
   }
   return n;
}
int LHCalibration::GetMDcount(){
   int i,n=0;
   for(i=0;i<NNMD;i++){
      if(fMDcount[i]>0)n++;
   }
   return n;
}
int LHCalibration::GetMDcountGood(){
   int i,n=0;
   for(i=0;i<NNMD;i++){
      if(fMDcount[i]>0&&fMDstatus[i]==0)n++;
   }
   return n;
}
