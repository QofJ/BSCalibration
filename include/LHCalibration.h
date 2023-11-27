// -----------------------------------------------------------------------
//                                LHAASO
//                       --- LHCalibration ---
//                           class header file
// -----------------------------------------------------------------------
//  Description:
//
//    Manager of the LHAASO detector calibration and status
//
// -----------------------------------------------------------------------
// History:
//    Songzhan Chen    -- 2021.4 created
//           
//========================================================================
//
#ifndef __LHCALIBRATION_HH__
#define __LHCALIBRATION_HH__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "KM2AEvent.h"
#include "LHEvent.h"
#include "G4KM2A_Geometry.h"
#include "TNtupleD.h"
#include "TClonesArray.h"
#include "EOSopen.h"
#define NNED  5500  //maximum value to define vector   
#define NNMD  1500    
#define NNWCDA 3500   

using namespace std;

class LHCalibration
{
  public:
     static LHCalibration * GetInstance()
      {
       if (m_myself == 0) m_myself = new LHCalibration();
       return m_myself;
      }
     ~LHCalibration();

  private:
     LHCalibration();
     static LHCalibration * m_myself;
     static G4KM2A_Geometry *geom;
  public:
     int SetTDCcalibrationKM2A(const char filename[100]);
     int SetADCcalibrationKM2A(const char filename[100]);
     int SetMDcalibration(const char filename[100]);
     int SetStatusKM2A(const char filename[100],int cflag);
     int SetTDCcalibrationWCDA(const char filename[100]);
     int SetADCcalibrationWCDA(const char filename[100]);
     int SetFinalStatusKM2A();
     int SetStatusWCDA(const char filename[100]);
     int ApplyKM2AEvent(KM2AEvent *ev,LHEvent *lhev);
     int ApplyWCDAEvent();
     //return number of good detector 
     int GetNED(){return fNED;};
     int GetNMD(){return fNMD;};
     int GetNWCDA(){return fNWCDA;};

     float GetEDdt(int n){return fEDdt[n];};
     float GetMDdt(int n){return fMDdt[n];}; 
     float GetWCDAdt(int n){return fWCDAdt[n];}; 
     float GetWCDA2dt(int n){return fWCDA2dt[n];};

     float GetEDratio(int n){return fEDratio[n];};
     float GetMDratio(int n){return fMDratio[n];};
     float GetWCDAratio(int n){return fWCDAratio[n];};
     float GetWCDA2ratio(int n){return fWCDA2ratio[n];};

     float GetEDda(int n){return fEDda[n];};
     float GetMDda(int n){return fMDda[n];};
     float GetWCDAda(int n){return fWCDAda[n];};
     float GetWCDA2da(int n){return fWCDA2da[n];};

     int GetEDstatus(int n){return fEDstatus[n];};
     int GetMDstatus(int n){return fMDstatus[n];};
     int GetWCDAstatus(int n){return fWCDAstatus[n];};
     int GetWCDA2status(int n){return fWCDA2status[n];};
     int GetEDcount();
     int GetMDcount();
     int GetEDcountGood();
     int GetMDcountGood();

   private:
     int    fNED,fNMD,fNWCDA,fNWCDA2;
     string KM2ATDCfilename,KM2AADCfilename,MDfilename,KM2AStatusfilename;
     string WCDATDCfilename,WCDAADCfilename,WCDAStatusfilename;
     std::vector<float>    fEDdt;
     std::vector<float>    fEDda;
     std::vector<float>    fEDratio;
     std::vector<int>      fEDstatus;
     std::vector<int>      fEDcount;

     std::vector<float>    fMDdt;
     std::vector<float>    fMDda;
     std::vector<float>    fMDratio;
     std::vector<int>      fMDstatus;
     std::vector<int>      fMDcount;

     std::vector<float>    fWCDAdt;
     std::vector<float>    fWCDAda;
     std::vector<float>    fWCDAratio;
     std::vector<int>      fWCDAstatus;

     std::vector<float>    fWCDA2dt;
     std::vector<float>    fWCDA2da;
     std::vector<float>    fWCDA2ratio;
     std::vector<int>      fWCDA2status;

};

#endif  
