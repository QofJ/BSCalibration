// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIsrcdILHRecEventDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "LHRecEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LHRecEvent(void *p = 0);
   static void *newArray_LHRecEvent(Long_t size, void *p);
   static void delete_LHRecEvent(void *p);
   static void deleteArray_LHRecEvent(void *p);
   static void destruct_LHRecEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LHRecEvent*)
   {
      ::LHRecEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LHRecEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LHRecEvent", ::LHRecEvent::Class_Version(), "LHRecEvent.h", 7,
                  typeid(::LHRecEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LHRecEvent::Dictionary, isa_proxy, 4,
                  sizeof(::LHRecEvent) );
      instance.SetNew(&new_LHRecEvent);
      instance.SetNewArray(&newArray_LHRecEvent);
      instance.SetDelete(&delete_LHRecEvent);
      instance.SetDeleteArray(&deleteArray_LHRecEvent);
      instance.SetDestructor(&destruct_LHRecEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LHRecEvent*)
   {
      return GenerateInitInstanceLocal((::LHRecEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LHRecEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LHRecEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LHRecEvent::Class_Name()
{
   return "LHRecEvent";
}

//______________________________________________________________________________
const char *LHRecEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHRecEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LHRecEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHRecEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LHRecEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHRecEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LHRecEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHRecEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LHRecEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class LHRecEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LHRecEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(LHRecEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LHRecEvent(void *p) {
      return  p ? new(p) ::LHRecEvent : new ::LHRecEvent;
   }
   static void *newArray_LHRecEvent(Long_t nElements, void *p) {
      return p ? new(p) ::LHRecEvent[nElements] : new ::LHRecEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_LHRecEvent(void *p) {
      delete ((::LHRecEvent*)p);
   }
   static void deleteArray_LHRecEvent(void *p) {
      delete [] ((::LHRecEvent*)p);
   }
   static void destruct_LHRecEvent(void *p) {
      typedef ::LHRecEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LHRecEvent

namespace {
  void TriggerDictionaryInitialization_LHRecEventDict_Impl() {
    static const char* headers[] = {
"LHRecEvent.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/lhaaso.ihep.ac.cn/anysw/slc5_ia64_gcc73/external/root/6.22.00/include/",
"/workfs2/ybj/chensz/software/DataRec/KM2Arec_3quarter_V2/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LHRecEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$LHRecEvent.h")))  LHRecEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LHRecEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "LHRecEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"LHRecEvent", payloadCode, "@",
"StrDup", payloadCode, "@",
"operator+", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LHRecEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LHRecEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LHRecEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LHRecEventDict() {
  TriggerDictionaryInitialization_LHRecEventDict_Impl();
}
