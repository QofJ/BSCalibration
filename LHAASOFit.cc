/*  2016-04-18: reconstruction code for LHAASO fastMC
 *  2017-11-08: modified for G4KM2A geant4 simulation
 *  2019-09-24: modified for KM2A data 
 *  2020-02-04: modified for trigger and spacetimefilter
 *  2020-06-05: add detector status result
 *  2021-05-20: add LHCalibraiton class
 *  If you find any bug please send email to
 *      chensz@ihep.ac.cn  !!!
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TNtupleD.h"
#include "KM2AEvent.h"
//#include "EOSopen.h"
#include "LHCalibration.h"
#include "G4KM2A_Reconstruction.h"
#include "EOSopen.h"
int main(int argc, char *argv[]) {
    int i,j,k,n,m,id,nentries,nTotalHits,mode,tag,nflag,nFile,Nout,arrayFlag,Num[4];
    double t,pe,x,y,z,dt,mjd,pd,peda,pedd,Tresolutioni,ntu[30],Dnum[10];	 
    float dtt,daa,drr;
    char sTemp[80],name[50];
    if(argc<4) {
        printf("Usage: %s arrayflag outfile in_file [ADCcalfile TDCcalfile DetectorStatus] \n", argv[0]);
        printf("arrayflag: 0:unit; 1:ybj array; 2:33ED; 3:71ED+10MD, 4:1/4 array; 5:1/2 array  6: LHAASO all\n");
        return 1;
    }
    arrayFlag=atoi(argv[1]);
    printf("Flag=%d\n",arrayFlag);

    printf("***********************************************\n");
    printf("***  LHAASO reconstruction v2.1 @2021-06-16 ***\n");
    printf("***  Please email to chensz@ihep.ac.cn for  ***\n");
    printf("***  any problem or bug !                   ***\n");
    printf("***********************************************\n");
    //init the detector position firstly !!!
    G4KM2A_Geometry *geom= G4KM2A_Geometry::GetInstance(arrayFlag);
    //Calibration class 
    LHCalibration *cal= LHCalibration::GetInstance();
    if(argc>4){if(cal->SetADCcalibrationKM2A(argv[4])) return 1;     } 
    if(argc>5){if(cal->SetTDCcalibrationKM2A(argv[5])) return 1;     }
    if(argc>7){if(cal->SetMDcalibration(argv[7])) return 1;     }
    if(argc>6){if(cal->SetStatusKM2A(argv[6],1)) return 1;     }
    if(argc>4)cal->SetFinalStatusKM2A();
    char Recversion[300],Adccalname[300]="none",Tdccalname[300]="none",Statusname[300]="none",MDname[300]="none";
    strcpy(Recversion,"LHAASO reconstruction v2.1 @2021-06-16");
    if(argc>4)strcpy(Adccalname,argv[4]);
    if(argc>5)strcpy(Tdccalname,argv[5]);
    if(argc>6)strcpy(Statusname,argv[6]);
    if(argc>7)strcpy(MDname,argv[7]);
    //output file and tree
    LHRecEvent *rec = new LHRecEvent();
    
    G4KM2A_Reconstruction *reconstruct=G4KM2A_Reconstruction::GetInstance(arrayFlag);
    reconstruct->SetDetetorStatus(cal);
    //decode data event
    KM2AEvent *event = new KM2AEvent();
    //rec event
    LHEvent *lhevent =new LHEvent();
    //output root file
    TFile *newfile = EOSopen(argv[2],1); //0 open, 1: created
    TNtupleD *nt1 = new TNtupleD("nt1", "LHAASO Reconstructed Data",
    "Nevent:mjd:NhitE:NfiltE:NhitM:NfiltM:NpE1:NpE2:NpE3:NuM1:NuM2:NuM3:NuM4:NuM5:r_Corex:r_Corey:r_Theta:r_Phi:r_C0:a:chi:NtrigE:NtrigE2:size:age:dr:nMD50");
   int nED,nMD,nEDgood,nMDgood;
    nEDgood=cal->GetNED();
    nMDgood=cal->GetNMD();
   TTree *info=new TTree("info","rec information");
      info->Branch("nED",&nED,"nED/I");
      info->Branch("nMD",&nMD,"nMD/I");
      info->Branch("nEDgood",&nEDgood,"nEDgood/I");
      info->Branch("nMDgood",&nMDgood,"nEDgood/I");
      info->Branch("version",Recversion,"version/C");
      info->Branch("ADCcal",Adccalname,"ADCcal/C");
      info->Branch("TDCcal",Tdccalname,"TDCcal/C");
      info->Branch("MDcal",MDname,"MDcal/C"); 
      info->Branch("Status",Statusname,"Status/C");

  //import the files and reconstructed them
    for(nFile=3;nFile<4;nFile++){
        TFile *hfile=EOSopen(argv[nFile],0);
        printf("inputing file -> %s\n",argv[nFile]);
        if(!hfile||hfile->IsZombie()||hfile->GetEND()<50){
            printf("file error %s!!\n",argv[nFile]);
            continue;
        }
        //input file and tree
        TTree *EventTree = (TTree *)gDirectory->Get("event");
        EventTree->SetBranchAddress("Event", &event);
        nentries = Int_t(EventTree->GetEntriesFast());
        printf("Total event %d\n",nentries);
   
        Nout=int(nentries/10.);
        //Nout=10000;
        for(i=0; i<1000; i++) {
//        for(i=0; i<nentries; i++) {
            EventTree->GetEntry(i);
            nflag=cal->ApplyKM2AEvent(event,lhevent);

            //error
            if (nflag < 2)continue;

            reconstruct->eventrecline(lhevent,rec,nflag); 
            ntu[0]=event->EvN();
            ntu[1]=event->Mjd();
            ntu[2]=event->NHit();
            ntu[3]=rec->NfiltE;
            ntu[4]=event->Nmd();
            ntu[5]=rec->NfiltM;
            ntu[6]=rec->NpE1;
            ntu[7]=rec->NpE2;
            ntu[8]=rec->NpE3;
            ntu[9]=rec->NuM1;
            ntu[10]=rec->NuM2;
            ntu[11]=rec->NuM3;
            ntu[12]=rec->NuM4;
            ntu[13]=rec->NuW1;
            ntu[14]=rec->rec_x;
            ntu[15]=rec->rec_y;
            ntu[16]=rec->rec_theta;
            ntu[17]=rec->rec_phi;
            ntu[18]=rec->rec_c0;
            ntu[19]=rec->rec_a;
            ntu[20]=rec->rec_sigma;
            ntu[21]=rec->NtrigE;
            ntu[22]=rec->NtrigW;
            ntu[23]=rec->rec_Esize;
            ntu[24]=rec->rec_Eage;
            ntu[25]=rec->NpW;
            ntu[26]=rec->NhitW; 
            nt1->Fill(ntu); 
            if(i%Nout==0||i==nentries-1)printf("---------->Reconstructed event %d %d\n",i,k);
        
        }
        nED=cal->GetEDcount();
        nMD=cal->GetMDcount();
        nEDgood=cal->GetEDcountGood();
        nMDgood=cal->GetMDcountGood();
        printf("live nED=%d nMD=%d \t good nED=%d  nMD=%d\n",nED,nMD,nEDgood,nMDgood);
        info->Fill(); 
        hfile->Close();
    }
    printf("----------End of Globle Reconstruction--------------\n");
    newfile->Write();
    newfile->Close();
    return 0;
}

