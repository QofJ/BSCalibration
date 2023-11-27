/*****************************************
 * Songzhan Chen    -- 2021-06-07 created
 * for EOS file open
 *
 ****************************************/
#ifndef __EOSOPEN_HH__
#define __EOSOPEN_HH__

#include "TFile.h"

TFile *EOSopen(const char filename[100],int Flag);
/*{
    char Name1[300]="root://eos01.ihep.ac.cn/";
    char Name2[300];
    char name[100];
    TFile *hfile;
    sscanf(filename,"%4s",name);
    if(strcmp(name,"/eos")==0){
        strcpy(Name2,Name1);
        strcat(Name2,filename);
    }
    else {
        strcpy(Name2,filename);
    }
    if(Flag==1){
        hfile = TFile::Open(Name2, "recreate");
        hfile->SetCompressionSettings(101);
    }
    else hfile = TFile::Open(Name2);
 
    if(!hfile){
        printf("file open error %s !!\n",filename);
    }
    if(hfile->IsZombie()){
        printf("file error %s !!\n",filename);
        hfile->Close();
    }
    return hfile;
}*/
#endif 
