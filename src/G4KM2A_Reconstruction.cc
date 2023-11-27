/*  main reconstruction code for LHAASOFM
 *  If you find any bug please send email to
 *  chensz@ihep.ac.cn
 *    2016-04-18
 *    update 2020-02-04 
 */
#include "G4KM2A_Reconstruction.h"
#include  <map>
#define DEBUG 0
#define PETH 0.2
using namespace std;

TMinuit *G4KM2A_Reconstruction::minuit = new TMinuit(4);
G4KM2A_Reconstruction * G4KM2A_Reconstruction::m_myself = 0;
G4KM2A_Geometry * G4KM2A_Reconstruction::geom=0;
float G4KM2A_Reconstruction::nED[]={0};
float G4KM2A_Reconstruction::nMD[]={0};
float G4KM2A_Reconstruction::nWCDA[]={0};
double G4KM2A_Reconstruction::_np=0;
double G4KM2A_Reconstruction::_dirl=0;
double G4KM2A_Reconstruction::_dirm=0;
double G4KM2A_Reconstruction::_dirn=0;
char G4KM2A_Reconstruction::_style[20]="ED";

G4KM2A_Reconstruction::G4KM2A_Reconstruction(int tFlag)
{
    int ned,nmd,nwc,i;
    geom = G4KM2A_Geometry::GetInstance(tFlag);
    ned=geom->GetNED();
    nmd=geom->GetNMD();
    nwc=geom->GetNWCDA();
 
    for(i=0;i<ned;i++)  nED[i]=0;
    for(i=0;i<nmd;i++)  nMD[i]=0;
    for(i=0;i<nwc;i++)  nWCDA[i]=0;  

}

int G4KM2A_Reconstruction::SetDetetorStatus(LHCalibration *cal)
{
    int i,mode,id,stat,ned,nmd;
    double x,y,z;

    ned=geom->GetNED();
    nmd=geom->GetNMD();
    for(i=0;i<ned;i++){
        id=geom->GetEDId(i);
        stat=cal->GetEDstatus(id);
        if(stat==0)nED[i]=0;
        else nED[i]=-2;
     }
    for(i=0;i<nmd;i++){
        id=geom->GetMDId(i);
        stat=cal->GetMDstatus(id);
        if(stat==0)nMD[i]=0;
        else nMD[i]=-2;
    }
    return 0; 
}
int G4KM2A_Reconstruction::SetDetetorStatus(const char file[500],int out[4])
{
   int i,mode,id,stat,ned,nmd;
   double x,y,z;
    FILE *statfile;
    statfile=fopen(file,"r");
    if(statfile)printf("using detector status file %s\n",file);
    else printf("Error: Cannot open file %s\n",file);
    
    ned=geom->GetNED();
    nmd=geom->GetNMD();
    for(i=0;i<ned;i++)  nED[i]=-2;
    for(i=0;i<nmd;i++)  nMD[i]=-2;
    for(i=0;i<4;i++) out[i]=0;
    fscanf(statfile, "%*s %*s %*s %*s %*s %*s\n");
    while(!feof(statfile)){
       i=fscanf(statfile, "%d %d %d %*lf\n", &mode,&id, &stat);
       //printf("%d %d %d\n",mode,id,stat);
       if(i>2){
        if(mode==0){
             if(geom->Getxyz(id,x,y,z,1,"ED")<0)continue;
             if(stat==0){ nED[geom->GetEDId2(id)]=0; out[0]++;}
             else { nED[geom->GetEDId2(id)]=-1; out[1]++; }
        }
        else if(mode==1){
             if(geom->Getxyz(id,x,y,z,1,"MD")<0)continue;
             if(stat==0){ nMD[geom->GetMDId2(id)]=0; out[2]++;}
             else { nMD[geom->GetMDId2(id)]=-1; out[3]++;}
        }
       }
   }
   fclose(statfile);
   return 0;
} 

G4KM2A_Reconstruction::~G4KM2A_Reconstruction()
{

}
int G4KM2A_Reconstruction::Init(LHRecEvent *trec){
    trec->rec_theta=-10;trec->rec_phi=-10;
    trec->rec_Etheta_p=-10;  trec->rec_Ephi_p=-10;
    trec->rec_Wtheta_p=-10;  trec->rec_Wphi_p=-10;
    trec->rec_Etheta_c=-10;  trec->rec_Ephi_c=-10;
    trec->rec_Wtheta_c=-10;  trec->rec_Wphi_c=-10;

    trec->rec_x=-10000  ; trec->rec_y=-10000;
    trec->rec_Ex=-10000 ; trec->rec_Ey=-10000;  trec->rec_Ez=-10000;
    trec->rec_Wx=-10000 ; trec->rec_Wy=-10000;  trec->rec_Wz=-10000;

    trec->rec_Esize=-10; trec->rec_Wsize=-10;
    trec->rec_Eage=-10;  trec->rec_Wage=-10;
    trec->rec_sigma=-1;
    trec->rec_a=0; 

    trec->NtrigE=-1 ;  trec->NtrigW=-1;
    trec->NfiltE=-1 ;  trec->NfiltM=-1;  trec->NfiltW=-1;
    trec->NpW=-1;
    trec->NpE1=-1; trec->NpE2=-1;trec->NpE3=-1;
    trec->NuM1=-1; trec->NuM2=-1;trec->NuM3=-1; trec->NuM4=-1;
    trec->NuW1=-1; trec->NuW3=-1;trec->NuW3=-1;
    trec->NfiltE =-1; trec->NfiltM =-1; trec->NfiltW =-1;
   return 0;
}
int G4KM2A_Reconstruction::eventrecline(LHEvent *tevent,LHRecEvent *trec,int nflag=10){
     int debug=DEBUG;
    TClonesArray *HitsE, *HitsM, *HitsW;
    //LHHit *tHit;
    //set parameters for this class
    //event=tevent;
    //rec=trec;
    HitsE = tevent->GetHitsE();
    HitsM = tevent->GetHitsM();
    HitsW = tevent->GetHitsW();
    //init the trec values
    Init(trec);
    //to exclude the half showers
    if(nflag<2)return 1;

    int ned=geom->GetNED();
    //mask the bad detector hit
    SetBadHit(*HitsE,tevent->GetNhitE(),"ED");
    SetBadHit(*HitsM,tevent->GetNhitM(),"MD");
   
    if(debug==1)     printf("timefilter\n");
    trec->NtrigE =  trigger(*HitsE,tevent->GetNhitE(),400,"ED",0);
    trec->NhitE =tevent->GetNhitE();
    trec->NhitM =tevent->GetNhitM();
    //trec->NhitW =tevent->GetNhitW();
    //trec->NtrigW=tevent->GetNtrigW();
    trec->NtrigW= trigger(*HitsE,tevent->GetNhitE(),400,"ED",1);
    /* for mc
    trec->dt    =tevent->GetDt();
    trec->id    =tevent->GetId();
    trec->E     =tevent->GetE();
    trec->phi   =tevent->GetPhi();
    trec->theta =tevent->GetTheta();
    trec->corex =tevent->GetCorex();
    trec->corey =tevent->GetCorey();
   //2020-02-14 wangly simulation ooo> World x y z zmax= 674.594 667.843 18.8665 27.189
   //chensz special for wangly simulation
    double ta=tan(trec->theta)*(27.189-18.8665);
    trec->corex =tevent->GetCorex()+ta*cos(trec->phi);
    trec->corey =tevent->GetCorey()+ta*sin(trec->phi);
    trec->pNpE  =tevent->GetNpE();
    trec->pNpW  =tevent->GetNpW();
    trec->pNuM  =tevent->GetNuM();
    trec->pNuW  =tevent->GetNuW();
    trec->pNuW2 =tevent->GetNuW2();
   */
    int i,j,id,tWind,rWind;
    double aa,dr;
    
     tWind=400;rWind=100;
    if(debug==1)printf("setnparticle\n");
    setnparticle(*HitsE,trec->NhitE,40.); //40.  set 0.2 ne to 0
    setnparticle(*HitsM,trec->NhitM,75.); //75.  set 0.2 nu to 0

    if(trec->NhitE>4*ned){
        printf("-----> %d  resetED   <-------\n",trec->NhitE);
        resetED(*HitsE,trec->NhitE,1);
    }
    //filter the hits using a time windows of tWind ns(using a step of 50ns) and radius windows of rWind m
    if(debug==1)printf("spacetimefilter\n");
    //printf("\t nhitE=%d \n",trec->NhitE);
    trec->NfiltW = spacetimefilter(*HitsE,trec->NhitE,tWind,rWind,"ED",0);  //return maxt
    if(debug==1)printf("%d %d\n",trec->NtrigE,trec->NfiltW);
    //reconstruct the direction using planar fit method
    if(debug==1)printf("planarfit\n");
    if(trec->NfiltW>=3&&planarfit(*HitsE,trec->NhitE,trec,"ED")==0){
        trec->rec_theta=trec->rec_Etheta_p;
        trec->rec_phi=trec->rec_Ephi_p;
        trec->rec_c0=trec->rec_Ec0_p;
        trec->rec_sigma=trec->rec_Esigma_p;
        trec->rec_a =0;

        if(debug==1)printf("core_centre2\n");
        //reconstruct core using centrial method 
        core_centre2(*HitsE,trec->NhitE,rWind,trec,"ED");
            dr=trec->rec_Ez*tan(trec->rec_theta);
            trec->rec_x=trec->rec_Ex+dr*cos(trec->rec_phi);
            trec->rec_y=trec->rec_Ey+dr*sin(trec->rec_phi);
            trec->rec_Wx=trec->rec_Ex+dr*cos(trec->rec_phi);
            trec->rec_Wy=trec->rec_Ey+dr*sin(trec->rec_phi);
            if(DEBUG)printf("core: %lf %lf\n",trec->rec_x,trec->rec_y);
        //reconstruct the direction using coincial fit method, fix a=0.03
        if(debug==1)printf("conicalfit\n");
        
        if(conicalfit(*HitsE,trec->NhitE,trec,0.035,"ED")==0){
            trec->rec_theta=trec->rec_Etheta_c;
            trec->rec_phi=trec->rec_Ephi_c;
            trec->rec_a =trec->rec_Ea;
            trec->rec_c0 =trec->rec_Ec0_c;
            trec->rec_sigma =trec->rec_Esigma_c;

            //filter out the noise using core and direction information
            trec->NfiltE = noisefilter(*HitsE,trec->NhitE,-50,100,rWind+100,trec,"ED"); 
            //reconstruct the core and size, age using the NKG and likelihood
            if(trec->NfiltE>=10){
               if(trec->NfiltE>200)noisefilter(*HitsE,trec->NhitE,-50,100,rWind+300,trec,"ED");
               core_likelihood(*HitsE,trec->NhitE,trec,"ED",0);
               trec->rec_x=trec->rec_Ex;
               trec->rec_y=trec->rec_Ey;
             }
            if(debug==1)printf("conicalfit 2\n");
            if(trec->NfiltE>200)aa=-0.035; //free a
            else    aa=0.035;
            if(conicalfit(*HitsE,trec->NhitE,trec,aa,"ED")==0){
                trec->rec_theta=trec->rec_Etheta_c;
                trec->rec_phi=trec->rec_Ephi_c;
                trec->rec_a =trec->rec_Ea;
                trec->rec_c0 =trec->rec_Ec0_c;
                trec->rec_sigma =trec->rec_Esigma_c;
            }
            //get sum particle number from radius r= 0- 100m, dt:-50,100
            noisefilter(*HitsE,trec->NhitE,-30,50,100,trec,"ED");
            trec->NpE1 = getnp(*HitsE,trec->NhitE,0,trec,"ED");
            //get sum particle number from radius 40m-100m, dt:-50,100
            trec->NpE2 = getnp(*HitsE,trec->NhitE,40,trec,"ED");
            //get sum particle number from radius 20m-100m, dt:-50,100
            noisefilter(*HitsE,trec->NhitE,-30,50,200,trec,"ED");
            trec->NpE3 = getnp(*HitsE,trec->NhitE,0,trec,"ED");

            //Float_t NuM1;  //number of muons of MD with r=15-200, dt:-50,100 
            trec->NfiltM = noisefilter(*HitsM,trec->NhitM,-30,50,200,trec,"MD");
            trec->NuM1 = getmuM(*HitsM,trec->NhitM,15,trec); //now the value 15m is set according to MC 
 
            noisefilter(*HitsM,trec->NhitM,-30,50,100,trec,"MD");
            trec->NuM2 = getmuM(*HitsM,trec->NhitM,40,trec);
  
            trec->NfiltM = noisefilter(*HitsM,trec->NhitM,-30,50,200,trec,"MD");
            trec->NuM3 = getmuM(*HitsM,trec->NhitM,40,trec);
                                    
            noisefilter(*HitsM,trec->NhitM,-30,50,400,trec,"MD");
            trec->NuM4 = getmuM(*HitsM,trec->NhitM,15,trec);

            noisefilter(*HitsM,trec->NhitM,-30,50,800,trec,"MD");
            trec->NuW1 = getmuM(*HitsM,trec->NhitM,15,trec);
      
         //   noisefilter(*HitsM,trec->NhitM,-30,50,400,trec,"MD");
         //   trec->NuW2 = getmuM(*HitsM,trec->NhitM,5,trec);

            //get the live MD number within 15-50m from core
            trec->NhitW=getliveMD(trec,15,50);

            trec->NpW=geom->GetCoreEdge(trec->rec_x,trec->rec_y);
            //chensz for test 2020-04-03
          /*
            noisefilter(*HitsE,trec->NhitE,-30,50,400,trec,"ED");
            setresidual(*HitsE,trec->NhitE,trec,"ED"); 
            noisefilter(*HitsM,trec->NhitM,-30,50,400,trec,"MD");
            setresidual(*HitsM,trec->NhitM,trec,"MD");
         */
        }
    }
    return 0;
}
//get the live MD within distance range
int G4KM2A_Reconstruction::getliveMD(LHRecEvent *trec,double rmin,double rmax){
    int i,j,k,nmd,ngood=0;
    nmd=geom->GetNMD();
    double dirl,dirm,dirn,dx,dy,dz,dr,x,y,z;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);

    for(i=0;i<nmd;i++){
        if(nMD[i]<0)continue;   
        if(geom->Getxyz(i,x,y,z,0,"MD")<0)continue;     
        dx = x-trec->rec_x;
        dy = y-trec->rec_y;
        dz = z;
        dr  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
        if(dr>rmin&&dr<rmax) ngood++;
    }
    return ngood;
}
//set the hit status to -1 for bad detector
int G4KM2A_Reconstruction::SetBadHit(TClonesArray &tHits,int np,const char *style){
    int i,id;
    LHHit *tHit;
    if (strcmp(style, "ED") == 0) {
        for(i=0; i<np; i++){
            tHit = (LHHit *)((tHits)[i]);
            id = tHit->GetId();
            if(nED[geom->GetEDId2(id)]<0)tHit->SetStatus(-2);
        }
    }
    else if(strcmp(style, "MD") == 0){
        for(i=0; i<np; i++){
            tHit = (LHHit *)((tHits)[i]);
            id = tHit->GetId();
            if(nMD[geom->GetMDId2(id)]<0)tHit->SetStatus(-2);
        }
    }
    return 0;
}
//get the sum number of particles for ED array
float G4KM2A_Reconstruction::getnp(TClonesArray &tHits,int np, int rcut, LHRecEvent *trec, const char *style){
    if(np<4)return 0.;
    LHHit *tHit;
    int i;
    float dn=0.;
    double dirl,dirm,dirn,dx,dy,dz,dr,x,y,z;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>2.5){
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)>0){
              dx = x-trec->rec_x;
              dy = y-trec->rec_y;
              dz = z;
              dr  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
              if(dr>rcut) dn += tHit->GetPe();
            }
        }
   }
   return dn;
}

//reset event
int G4KM2A_Reconstruction::resetED(TClonesArray &tHits,int np, int ncut){
    LHHit *tHit,*tHit2;
    int i,j,id;
    double t,pe;
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        id = tHit->GetId();
        t = tHit->GetTime();
        pe = tHit->GetPe();
        if(t<4500||t>9000||pe<ncut)continue; 
        for(j=i; j<np; j++){
            tHit2 = (LHHit *)((tHits)[j]);
            if(id== tHit2->GetId()&&tHit2->GetPe()>0.2){
               if(fabs(tHit2->GetTime()-t-430)<50){
                  pe += tHit2->GetPe();
                  tHit->SetPe(pe);
                  tHit2->SetPe(0);
               } 
            }
        }
   }
   return 0;
} 
//get the muon number of MD detector
float G4KM2A_Reconstruction::getmuM(TClonesArray &tHits,int np,int rcut,LHRecEvent *trec){
    LHHit *tHit;
    int i;
    float dn=0;
    double dirl,dirm,dirn,dx,dy,dz,dr,x,y,z;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta); 
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
        if(geom->GetMDxyz(tHit->GetId(),x,y,z,1)<0)continue;
        dx = x-trec->rec_x;
        dy = y-trec->rec_y;
        dz = z;
        dr  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
        if(dr>rcut) dn += tHit->GetPe();
   }
   return dn;
}
//get the muon number of WCDA detector
float G4KM2A_Reconstruction::getmuW(TClonesArray &tHits,int np,int rcut,LHRecEvent *trec){
    LHHit *tHit;
    int i;
    double dirl,dirm,dirn,dx,dy,dr,x,y,z,dn=0;

    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);

    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<2.5)continue;
        if(tHit->GetPe()>16){
            if(geom->GetWCDAxyz(tHit->GetId(),x,y,z,1)<0)continue;
            dx = x-trec->rec_x;
            dy = y-trec->rec_y;
            
            dr  = sqrt(dx*dx+dy*dy+z*z-pow(dx*dirl+dy*dirm-z*dirn,2.));
            if(dr>rcut) dn += tHit->GetPe()/45.;//from MC fit
        }  
   }
   return dn;
}
int G4KM2A_Reconstruction::setnparticle(TClonesArray &tHits,int np, double pe){
    LHHit *tHit;
    int i;
    double ne;
    for(i=0;i<np;i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>-1){
            ne= tHit->GetPe()/pe;
            if(ne<PETH)ne=0;
            tHit->SetPe(ne);
         }
    }
    return 1;
}

int G4KM2A_Reconstruction::trigger(TClonesArray &tHits,int np,int twind,const char *style,int OFlag=0){
    int i,j,t,flag,flag2,tflag,maxid,maxt,rmax;
    double a,b,x,y,z,r;
    LHHit *tHit;
    int id,id2;
    int Nbin,Nwind,Nstep,Nstar,tdc[10000],status[6000];
    std::map<int, int> TriggerED;
    Nwind=twind;
    if (strcmp(style, "WCDA") == 0) { Nbin = 2200 - 400; Nstep = 1; Nstar = 600; }
    else {Nbin=10000-1000; Nstep=1;Nstar=4000;}
    for(i=0; i<Nbin; i++)tdc[i]=0;
    for(i=0;i<np;i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>-1){
            t = int(tHit->GetTime());
            if(t>=0&&t<Nbin)tdc[t]++;
        }
        //else printf("EDid=%d\n",tHit->GetId());
    }
    rmax=0;
    for(tflag=Nstar;tflag<Nbin-Nwind;tflag+=Nstep){
        while(tdc[tflag]==0&&tflag<Nbin-Nwind)tflag++;
        flag=0;
        for(j=tflag;j<tflag+Nwind;j++)flag +=tdc[j];
        if(flag>rmax&&flag>=5) {
            TriggerED.clear(); 
            for(i=0;i<np;i++){
               tHit = (LHHit *)((tHits)[i]);
               t = int(tHit->GetTime());
               if(tHit->GetStatus()>-1&&t>=tflag&&t<tflag+Nwind){
                   id= (tHit->GetId());
                   if(OFlag==1){
                      if(geom->Getxyz(id,x,y,z,1,style)<0)continue;
                      if(sqrt(x*x+y*y)>575)continue;
                   }
                   if(TriggerED.count(id)) TriggerED[id] ++;
                   else   TriggerED[id] = 1;
               }
            }
          }
          if(TriggerED.size()>rmax)
          rmax=TriggerED.size();
    }
    return rmax;
}

int G4KM2A_Reconstruction::timefilter(TClonesArray &tHits,int np,int twind,const char *style, int Flag){
    int i,j,t,id,flag,tflag,maxid,maxt,rmax;
    double a,b,x,y,z,r;
    LHHit *tHit;
    int Nbin,Nwind,Nstep,Nstar,tdc[1000];
    Nwind=twind/10;
    if (strcmp(style, "WCDA") == 0) { Nbin = 220 - 40; Nstep = 1; Nstar = 60; }
    else {Nbin=900; Nstep=5;Nstar=400;}

    for(i=0; i<Nbin; i++)tdc[i]=0;
    for(i=0;i<np;i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>-1){
            if(tHit->GetPe()<1.e-3)tHit->SetStatus(0);
            else{
               tHit->SetStatus(5);
               t = int(tHit->GetTime()/10);
               if(t>=0&&t<Nbin)tdc[t]++;
            }
         }
    }
    //search the max time window
    rmax=0;
    for(tflag=Nstar;tflag<Nbin-Nwind;tflag+=Nstep){
        while(tdc[tflag]==0&&tflag<Nbin-Nwind)tflag++;
        flag=0;
        for(j=tflag;j<tflag+Nwind;j++)if(j<Nbin-Nwind)flag +=tdc[j];
        if(flag>rmax) {rmax=flag;maxt=tflag;}
    }
    //filter the noise
    for(i=0;i<np;i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>0){
            if(tHit->GetTime()<maxt*10||tHit->GetTime()>=(maxt+Nwind)*10)tHit->SetStatus(0);
        }
    }
    return rmax;
}

//filter out the noise hits using detector position and TDC time
int G4KM2A_Reconstruction::spacetimefilter(TClonesArray &tHits,int np,int twind,int rwind,const char *style, int Flag){
    int i,j,t,id,tflag,nflag,maxid,maxt;
    double a,b,x,y,z,r;
    LHHit *tHit;
    int Nbin,Nstep,Nstar,Nall=0;
    double flag,rmax,tdc[10000];
    if (strcmp(style, "WCDA") == 0) {
         Nbin=2200-400;Nstep=10;Nstar=600;
    }
    else {
         Nbin=10000-1000; Nstep=50;Nstar=4000;
         //if(np> 5*(geom->GetNED())){
            //Nbin=10000-4600;  
         //   rwind=1000;
         //}
    }
    for(i=0; i<Nbin; i++)tdc[i]=0;
    for(i=0;i<np;i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>-1){
            if(tHit->GetPe()<1.e-3)tHit->SetStatus(0);
            else {
               tHit->SetStatus(5);
               t = int(tHit->GetTime());
               if(t>=0&&t<Nbin){ tdc[t]++; Nall++;}
               //if(t>=0&&t<Nbin){ tdc[t]+=tHit->GetPe(); Nall++;}
            }
         }
    }
    if(Nall<4)return 0;
    //to search for the maximum circle with radius rwind
    rmax=0;
    for(tflag=Nstar;tflag<Nbin-twind;tflag+=Nstep){
        while(tdc[tflag]==0&&tflag<Nbin-twind-1)tflag++;
       // if(rmax>3){
       //    flag=0;
       //    for(j=tflag;j<tflag+twind;j++)if(j<Nbin-twind)flag +=tdc[j];
       //    if(flag<rmax)continue;
       // }
        flag=0;
        nflag=0;
        for(i=0;i<np;i++){
            tHit = (LHHit *)((tHits)[i]);
            if(tHit->GetStatus()>0&&tHit->GetTime()>=tflag&&tHit->GetTime()<tflag+twind){
                if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
                a=x;
                b=y;
                id=tHit->GetId();
                flag=0; 
                nflag=0;
                for(j=0;j<np;j++){
                    tHit = (LHHit *)((tHits)[j]);
                    if(tHit->GetStatus()>0&&tHit->GetTime()>=tflag&&tHit->GetTime()<tflag+twind){
                        if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
                        r=sqrt(pow(x-a,2)+pow(y-b,2));
                        if(r<rwind){
                           flag += tHit->GetPe();
                           nflag++;
                         }
                    }
                }
                if(flag>rmax&&nflag>10) {rmax=flag;maxid=id;maxt=tflag;}
                if(nflag>rmax) {rmax=nflag;maxid=id;maxt=tflag;}   
            }
        }//end of np
    }//end of tflag
    //filter out noise
    /*
    if(rmax<4){
            printf("ramx %.2lf %d %d rwind=%d twin=%d \n",rmax,maxt,maxid,rwind,twind);
            for(j=0;j<np;j++){
                  tHit = (LHHit *)((tHits)[j]);
                  printf("   %d %.2lf %.4lf %d\n",tHit->GetId(),tHit->GetTime(),tHit->GetPe(),tHit->GetStatus()); 
                }
     }
    */
    if(maxid>0&&geom->Getxyz(maxid,x,y,z,1,style)>0){
      a=x;
      b=y;
     for(i=0;i<np;i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>0){
            if(tHit->GetTime()<maxt||tHit->GetTime()>=maxt+twind)tHit->SetStatus(0);
            else{
                if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)>0){
                    r=sqrt(pow(x-a,2)+pow(y-b,2));
                    if(r>rwind)tHit->SetStatus(0);
                }
            }
         }
     }
    }
    //return maxt;
    return (int)rmax;
}
//reconstruct the direction of ED or WCDA array using planar fit
int G4KM2A_Reconstruction::planarfit(TClonesArray &tHits,int np,LHRecEvent *trec,const char *style){
    if(np<4)return -1;
    int i,j,k,nus,Flag=0;
    double dt,w,c[3],rr[3]={0,0,0};
    double x,y,z,the,phi,sigma,sigmasum,costheta=0.94,widt=10.; 
    TMatrixD a(3,3); 
    TVectorD b(3),r(3);  
    LHHit *tHit;
    for(int iter=0;iter<15;iter++){
        for(i=0;i<3;i++){
            for(j=0;j<3;j++)a[i][j]=0;
            b[i]=0;
        }
        nus=0; sigmasum=0;
        for(i=0; i<np; i++){
            tHit = (LHHit *)((tHits)[i]);
            if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue; //not use the noise hits
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue; 
            if(iter==0)w=1.;
            else{
                dt = tHit->GetTime()-(rr[0]*x+rr[1]*y+rr[2]-costheta*z)/0.2998; //ns
                w  = exp(-0.5*pow(dt/widt,2)); //sig=10. is better than 5,15,20
                sigmasum +=dt*dt; 
                if(iter>8){
                    if(fabs(dt)>3*sigma&&fabs(dt)>30)tHit->SetStatus(1);
                }
            }
            nus++;
            //w*=sqrt(tHit->GetNp()); not good
            //w=sqrt(w); there is a very little change for this weight
            c[0] = x*w;
            c[1] = y*w;
            c[2] = 1*w;
            for(k=0;k<3;k++){
                for(j=0;j<3;j++) a[j][k] +=c[j]*c[k];
                b[k] += ((tHit->GetTime())*0.2998+costheta*z)*w*c[k];
            }
        }
        //solve the equation
        //if(nus<3)printf("%d\n",nus);
        if(nus<3)break; //the minimu required hit number is 4
        //printf("%d\n",nus);
        TDecompSVD svd(a);
        Bool_t ok;
        r = svd.Solve(b,ok);
        if(DEBUG)printf(" PlanarFit : %d %d %lf %lf %lf\n",iter,nus,r[0],r[1],r[2]);
        if(!ok)break;
        if((r[0]*r[0] + r[1]*r[1])>1){
             //printf(" PlanarFit >1: %d %d %lf %lf %lf\n",iter,nus,r[0],r[1],r[2]); 
             break;
        }
         if(ok){ 
           sigma =sqrt(sigmasum/(nus-2.));
           Flag=1;
           costheta=sqrt(1-(r[0]*r[0] + r[1]*r[1])); 
           rr[0]=r[0];
           rr[1]=r[1];
           rr[2]=r[2]; 
         }
   }
   //get the direction
   phi = -10.;  the = -10.;
   dt = rr[0]*rr[0] + rr[1]*rr[1];
   if(Flag>0.5 && dt<=1. && dt>=0.){
       the = asin(sqrt(dt));
       phi = atan2(rr[1],rr[0]);
       if(phi<0.) phi += PI*2.;
       if(DEBUG)printf("%lf %lf\n",the*57.3,phi*57.3);
   }
   
   //set the result
   if (strcmp(style, "ED") == 0) {
       trec->rec_Etheta_p = the;
       trec->rec_Ephi_p = phi;
       trec->rec_Ec0_p = rr[2];
       trec->rec_Esigma_p = sigma;
   }
   else if (strcmp(style, "WCDA") == 0) {
       trec->rec_Wtheta_p = the;
       trec->rec_Wphi_p   = phi;
       trec->rec_Wc0_p    = rr[2];
       trec->rec_Wsigma_p = sigma;
   }
   if(the>-1)return 0;
   else return -1;
}
//reconstruct the direction of ED or WCDA array using conical fit;
//if alpha>0 fix alpha, else fit the alpha
int G4KM2A_Reconstruction::conicalfit(TClonesArray &tHits,int np,LHRecEvent *trec,float alpha,const char *style){
    if(np<4)return -1;
    int i,j,k,nus,AA,Flag=0;
    double dt,w,the,phi,dirl,dirm,dirn,dr,c[4];
    double x,y,z,dx,dy,dz,sigma;
    if(alpha>0)AA=1;
    else AA=0; //to fit the alpha
    TMatrixD a(4,4);
    TVectorD b(4),r(4);
    LHHit *tHit;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(int iter=0;iter<30;iter++){
        for(i=0;i<4;i++){
            for(j=0;j<4;j++)a[i][j]=0;
            b[i]=0;
        }
        nus=0; sigma=0;
        for(i=0; i<np; i++){
            tHit = (LHHit *)((tHits)[i]);
            if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue; //not use the noise hit filter by planar fit
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
            dx = x-trec->rec_x; 
            dy = y-trec->rec_y;
            dz = z;
            dr  = (dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
            if(dr<0)dr=0;
            dr =sqrt(dr);
            w=1.;
            if(iter>0){
               // dt = tHit->GetTime()-(r[0]*x+r[1]*y-dirn*z+(r[2]+alpha*AA)*dr+r[3])/0.2998; //ns
                dt = tHit->GetTime()-(r[0]*x+r[1]*y-dirn*z+r[3])/0.2998; //ns please use this dt 
                if(dt>0||dt<-20) w  = w* exp(-0.5*pow(dt/10.,2)); //sig=10ns is better than 5,15,20 
                sigma +=dt*dt;
            } 
            nus++;
            //w=sqrt(w); has a little worse for this 
            //w*=sqrt(tHit->GetNp()); not good
            c[0] = x*w;
            c[1] = y*w;
            c[2] = dr*(1-AA)*w;
            c[3] = 1*w;
            for(k=0;k<4;k++){
                for(j=0;j<4;j++) a[j][k] +=c[j]*c[k];
                b[k] += ((tHit->GetTime())*0.2998-AA*alpha*dr+dirn*z)*w*c[k];
            }
            //if(DEBUG)printf("%d %lf %lf %lf %lf,%lf\n",i,c[0],c[1],tHit->GetTime(),dr,dirn*z);
        }
        //solve the equation
        if(nus<4)break;
        TDecompSVD svd(a);
        Bool_t ok;
        r = svd.Solve(b,ok);
        //if(DEBUG)printf(" ConicalFit : %d %d %lf %lf %lf\n",iter,nus,r[0],r[1],r[3]);
        if(!ok) break;
        if((r[0]*r[0]+r[1]*r[1])>1){
             //printf(" ConicalFit >1: iter=%d Nhit=%d %lf %lf %lf\n",iter,nus,r[0],r[1],r[3]);
             break;
        }
          if(ok){
             Flag=1;
             sigma =sqrt(sigma/(nus-3.));
             dirl=r[0];
             dirm=r[1];
             dirn=sqrt(1-(r[0]*r[0]+r[1]*r[1])); 
          }
    }
    //get the direction
    phi = -10.;  the = -10.;
    dt = r[0]*r[0] + r[1]*r[1];
    if(Flag>0.5 && dt<=1. && dt>=0.){
        the = asin(sqrt(dt));
        phi = atan2(r[1],r[0]);
        if(phi<0.) phi += PI*2.;
    }
    //set the result
    if (strcmp(style, "ED") == 0) {
        trec->rec_Etheta_c = the;
        trec->rec_Ephi_c = phi;
        trec->rec_Ea = r[2] * (1 - AA) + AA * alpha;
        trec->rec_Ec0_c = r[3];
        trec->rec_Esigma_c = sigma;
    }
    else if (strcmp(style, "WCDA") == 0) {
        trec->rec_Wtheta_c = the;
        trec->rec_Wphi_c   = phi;
        trec->rec_Wa       = r[2]*(1-AA)+AA*alpha;
        trec->rec_Wc0_c    = r[3]; 
        trec->rec_Wsigma_c = sigma;
    }
    if(the>=0) return 0;
    else return -1;
}
//reconstruct the core of ED or WCDA array using centroid method 
void G4KM2A_Reconstruction::core_centre(TClonesArray &tHits,int np,LHRecEvent *trec,const char *style){        
    int i,j,nus=0;
    double x,y,z,w,mx=0,my=0,mz=0,total=0;
    LHHit *tHit;
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
        nus++;
        if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
            // w = pow(tHit->GetNp(),2); //resolution is better than tHit->GetNp()
            w = pow(tHit->GetPe(),2); 
            if(sqrt(pow(x,2)+pow(y,2))>CORE_R) w *=4.; 
            mx    += x*w;
            my    += y*w;
            mz    += z*w;
            total += w;    
    }
    //set result 
    if (strcmp(style, "ED") == 0) {
        trec->rec_Ex = mx / total;
        trec->rec_Ey = my / total;
        trec->rec_Ez = mz / total;
    }
    else if (strcmp(style, "WCDA") == 0) {
        trec->rec_Wx = mx/total;
        trec->rec_Wy = my/total;
        trec->rec_Wz = mz/total;
    }
}
//reconstruct the core of ED or WCDA array using optimized centroid method 
int G4KM2A_Reconstruction::core_centre2(TClonesArray &tHits,int np,int nr,LHRecEvent *trec,const char *style){
    if(np<4)return -1;
    int i,j;
    double x,y,z,w,hx,hy,hz,mx,my,mz,total,dr,hx2=0.,hy2=0.;
    LHHit *tHit;
   // G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    for(int iter=0;iter<51;iter++){
        mx=0;my=0;total=0;mz=0; 
        for(i=0; i<np; i++){
            tHit = (LHHit *)((tHits)[i]);
            if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
            if(iter==0){ w=1;}
            else{
                dr = sqrt(pow(x-hx,2)+pow(y-hy,2));
                w  = exp(-0.5*pow(dr/15.,2)); //why use 15m?? can we optimized it???
            }
            //w *=(tHit->GetNp()); //Np is better than Np**2
            w *=(tHit->GetPe()); 
            //enhance the guard ring weight according the ED density 
            if(sqrt(pow(x,2)+pow(y,2))>CORE_R) w *=4.;
            mx    += x*w;
            my    += y*w;
            mz    += z*w;
            total += w;
        }
        hx = mx/total;
        hy = my/total;
        hz = mz/total;
        dr = sqrt(pow(hx-hx2,2)+pow(hy-hy2,2));
        if(dr<0.01&&iter>10)break;
        hx2 = hx;  hy2 = hy;
    }
    //set the result 
    if (strcmp(style, "ED") == 0) {
        trec->rec_Ex = hx;
        trec->rec_Ey = hy;
        trec->rec_Ez = hz;
    }
    else if (strcmp(style, "WCDA") == 0) {
        trec->rec_Wx = hx;
        trec->rec_Wy = hy;
        trec->rec_Wz = hz;
    }
    return 0;
}
//filter out noisze using direction and core result 
int G4KM2A_Reconstruction::noisefilter(TClonesArray &tHits,int np,int twind,int twinu,int rwin,LHRecEvent *trec,const char *style){
    int i,nus=0;
    double dt,dirl,dirm,dirn,dr,tt;
    double x,y,z,dx,dy,dz;
    LHHit *tHit;
    if (strcmp(style, "WCDA") == 0)tt = trec->dt;
    else tt=0;
    dirl=sin(trec->rec_Etheta_c)*cos(trec->rec_Ephi_c);
    dirm=sin(trec->rec_Etheta_c)*sin(trec->rec_Ephi_c);
    dirn=cos(trec->rec_Etheta_c);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<0||tHit->GetPe()<1.e-3)continue;
        if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
        dx = x-trec->rec_x;
        dy = y-trec->rec_y;
        dz = z;
        dr  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
        dt = tHit->GetTime()-tt-(dirl*x+dirm*y-dirn*z+trec->rec_Ea*dr+trec->rec_Ec0_c)/0.2998; 
        tHit->SetNp(int(dt*100.));
        if(dt>twind&&dt<twinu){
            if(dr<rwin){
                tHit->SetStatus(5);
                nus ++;
            }
            else tHit->SetStatus(2); 
        } 
        else{
            if(dr<rwin)tHit->SetStatus(1);
            else tHit->SetStatus(0);
        }
   }
   return nus;
}
int G4KM2A_Reconstruction::setresidual(TClonesArray &tHits,int np,LHRecEvent *trec,const char *style){
    int i,nus=0;
    double dt,dirl,dirm,dirn,dr,tt;
    double x,y,z,dx,dy,dz;
    LHHit *tHit;
    if (strcmp(style, "WCDA") == 0)tt = trec->dt;
    else tt=0;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        //if(tHit->GetStatus()<0||tHit->GetPe()<1.e-3)continue;
       // printf("\n%d %d ",i,tHit->GetId());
        if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
        dx = x-trec->rec_x;
        dy = y-trec->rec_y;
        dz = z;
        dr = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
        dt = tHit->GetTime()-tt-(dirl*x+dirm*y-dirn*z+trec->rec_a*dr+trec->rec_c0)/0.2998;
        //dt = tHit->GetTime()-tt-(dirl*x+dirm*y-dirn*z+trec->rec_c0)/0.2998; 
        tHit->SetNp(int(dr*100.));
        tHit->SetTime(dt);
       // printf("%lf ",dt);
   }
   //printf("\n");
   return nus;
}
//init the position of each detector
/*void G4KM2A_Reconstruction::init(const char *style){
    int i,j,k,ned,nmd,nwc;
    if(strcmp(style,"ED")==0){
       ned=geom->GetNED();
       for(i=0;i<ned;i++) nED[i]=0;
    }
    else if(strcmp(style,"WCDA")==0){
       nwc=geom->GetNWCDA();
       for(i=0;i<nwc;i++) nWCDA[i]=0;
    } 
}
*/
//NKG function1 for likelihood1 taken into account the detector with no hit 
void G4KM2A_Reconstruction::NKGfunction1(int &npar,double *gin,double &f,double *par,int iflag){
    int i,ned,nwc;
    //parameters use FengYouliang 2020-01-20 
   double x,y,z,r,u,dx,dy,dz,cs,A,rm=130.,sum=0; //the molliere radius is fixed to be 130m
    cs = TMath::Gamma(4.5-par[3])/(TMath::Gamma(par[3]-2.5+2.)*TMath::Gamma(4.5+2.5-2.0-2.*par[3]));
    cs = cs/(2*PI*rm*rm);
    ned=geom->GetNED();
    nwc=geom->GetNWCDA();
    if(strcmp(_style,"ED")==0){
        A=1.0*_dirn;
        for(i=0;i<ned;i++){
            if(nED[i]>-0.5){
                geom->Getxyz(i,x,y,z,0,_style);
                dx = x-par[0];
                dy = y-par[1];
                dz = z;
                r  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*_dirl+dy*_dirm-dz*_dirn,2.));
                if(r<0.2)r=0.2;
                u  = A*cs*par[2]*pow(r/rm,par[3]-2.5)*pow(1+r/rm,par[3]-4.5);  
                if(u>0)sum += (nED[i])*log(u)-u; 
            }
        }
    }
    else if(strcmp(_style,"WCDA")==0){
        A=25.0*_dirn;
        for(i=0;i<nwc;i++){
            if(nWCDA[i]>-0.5){
                geom->Getxyz(i,x,y,z,0,_style);
                dx = x-par[0];
                dy = y-par[1];
                dz = z;
                r  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*_dirl+dy*_dirm-dz*_dirn,2.));
                if(r<0.1)r=0.1;
                u  = A*cs*par[2]*pow(r/rm,par[3]-2.5)*pow(1+r/rm,par[3]-4.5);
                if(u>0)sum += nWCDA[i]*log(u)-u;
            }
        }
    }
    f = -sum;
}
//NKG function2 for likelihood2 only taken into account the detector with  hit 
void G4KM2A_Reconstruction::NKGfunction2(int &npar,double *gin,double &f,double *par,int iflag){
    int i; 
    double x,y,z,r,u,dx,dy,dz,cs,A,rm=130.,sum=0;
    if(strcmp(_style,"ED")==0)A=1.0*_dirn;
    else if(strcmp(_style,"WCDA")==0)A=25.0*_dirn;
    cs = TMath::Gamma(4.5-par[3])/(TMath::Gamma(par[3]-2.5+2)*TMath::Gamma(4.5+2.5-2.-2.*par[3]));
    cs = cs/(2*PI*rm*rm);
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
        
    TClonesArray* tHits=(TClonesArray*)minuit->GetObjectFit();
    LHHit *tHit;
    for(i=0; i<_np; i++){
        tHit = (LHHit *)((*tHits)[i]);
        if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
        geom->Getxyz(tHit->GetId(),x,y,z,1,_style);
        dx = x-par[0];
        dy = y-par[1];
        dz = z;
        r  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*_dirl+dy*_dirm-dz*_dirn,2.));
        //r  = sqrt(dx*dx+dy*dy-pow(dx*_dirl+dy*_dirm,2.));
        if(r<0.2)r=0.2;
        u  = A*cs*par[2]*pow(r/rm,par[3]-2.5)*pow(1+r/rm,par[3]-4.5);
        if(u>0)sum += (tHit->GetPe())*log(u)-u;
    }   
    f = -sum; 
}

// likelihood using NKG function1 taken into account the detector with and without hit
// using NKG function2 only taken into account the detector with  hit
int G4KM2A_Reconstruction::core_likelihood(TClonesArray &tHits,int np,LHRecEvent *trec,const char *style,int method){
    if(np<5)return -1;
    int  i,flag,ned,nwc;
    double size,x,y,sigma,esize,ex,ey,esigma,arglist[10];
    x=-1.e4;
    y=-1.e4;
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    ned=geom->GetNED();
    nwc=geom->GetNWCDA();
    if(np>4){
        minuit->mncler();
        //init the value will be used for NKGfunction
        _np=np;
        size=0;
        strcpy(_style,style);
        _dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
        _dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
        _dirn=cos(trec->rec_theta);
        //init all the detector hit value 
        if(strcmp(_style,"ED")==0){
            for(i=0;i<ned;i++){ if(nED[i]>-0.5)nED[i]=0;}
            LHHit *tHit;
            for(i=0; i<_np; i++){
                tHit = (LHHit *)((tHits)[i]);
                if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
                if(nED[geom->GetEDId2(tHit->GetId())]>-0.5)nED[geom->GetEDId2(tHit->GetId())] +=tHit->GetPe();
                size += tHit->GetPe();
            }
        }
        else if(strcmp(_style,"WCDA")==0){
            for(i=0;i<nwc;i++)    { if(nWCDA[i]>-0.5)nWCDA[i]=0;}
            LHHit *tHit;
            for(i=0; i<_np; i++){
                tHit = (LHHit *)((tHits)[i]);
                if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
                if(nWCDA[geom->GetWCDAId2(tHit->GetId())]>-0.5)nWCDA[geom->GetWCDAId2(tHit->GetId())] +=tHit->GetPe();
                size += tHit->GetPe();
            }
        }
        if(size<4)return -1;
        arglist[0]=-1;
        minuit->mnexcm("SET PRINT",arglist,1,flag);//1:standard; 2: try to improve minuim (slower)
        minuit->mnparm(0,"corex",(double)trec->rec_x,20.,-1000.,1000.,flag);
        minuit->mnparm(1,"corey",(double)trec->rec_y,20.,-1000.,1000.,flag);
        minuit->mnparm(2,"size", size*160, 100,size,size*1.e7,flag);
        minuit->mnparm(3,"age",1.2,0.2,0.5,2.5,flag); //for NKG
        //flag=0 if no problems,  >0 if MNPARM unable to implement definition
        minuit->FixParameter(3);
        if(method==0)minuit->SetFCN(NKGfunction1);
        else minuit->SetFCN(NKGfunction2);
        minuit->SetObjectFit(&tHits); 
        minuit->SetPrintLevel(-1); // -1  quiet (also suppresse all warnings),0  normal,1  verbose

        arglist[0]=500; //loop
        arglist[1]=0.1; //tolerance
        flag=0;
        minuit->mnexcm("SIMPLEX",arglist,0,flag);
        minuit->GetParameter(0,x,ex);
        minuit->GetParameter(1,y,ey);

        arglist[0]=1000; //loop
        arglist[1]=0.01; //tolerance
        minuit->FixParameter(0);
        minuit->FixParameter(1);
        minuit->Release(2);
        minuit->Release(3);
        minuit->mnexcm("MIGRAD",arglist,0,flag);
        minuit->GetParameter(2,size,esize);
        minuit->GetParameter(3,sigma,esigma);
    }
    if (strcmp(style, "ED") == 0) {
        trec->rec_Esize = float(size);
        trec->rec_Ex = float(x);
        trec->rec_Ey = float(y);
        trec->rec_Eage = float(sigma);
    }
    else if (strcmp(style, "WCDA") == 0) {
        trec->rec_Wsize =float(size);
        trec->rec_Wx = float(x);
        trec->rec_Wy = float(y);
        trec->rec_Wage = float(sigma);
    }
    return 0;
}
