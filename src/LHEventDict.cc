// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIsrcdILHEventDict
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
#include "LHEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LHEvent(void *p = 0);
   static void *newArray_LHEvent(Long_t size, void *p);
   static void delete_LHEvent(void *p);
   static void deleteArray_LHEvent(void *p);
   static void destruct_LHEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LHEvent*)
   {
      ::LHEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LHEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LHEvent", ::LHEvent::Class_Version(), "LHEvent.h", 12,
                  typeid(::LHEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LHEvent::Dictionary, isa_proxy, 4,
                  sizeof(::LHEvent) );
      instance.SetNew(&new_LHEvent);
      instance.SetNewArray(&newArray_LHEvent);
      instance.SetDelete(&delete_LHEvent);
      instance.SetDeleteArray(&deleteArray_LHEvent);
      instance.SetDestructor(&destruct_LHEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LHEvent*)
   {
      return GenerateInitInstanceLocal((::LHEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LHEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LHHit(void *p = 0);
   static void *newArray_LHHit(Long_t size, void *p);
   static void delete_LHHit(void *p);
   static void deleteArray_LHHit(void *p);
   static void destruct_LHHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LHHit*)
   {
      ::LHHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LHHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LHHit", ::LHHit::Class_Version(), "LHEvent.h", 153,
                  typeid(::LHHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LHHit::Dictionary, isa_proxy, 4,
                  sizeof(::LHHit) );
      instance.SetNew(&new_LHHit);
      instance.SetNewArray(&newArray_LHHit);
      instance.SetDelete(&delete_LHHit);
      instance.SetDeleteArray(&deleteArray_LHHit);
      instance.SetDestructor(&destruct_LHHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LHHit*)
   {
      return GenerateInitInstanceLocal((::LHHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LHHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LHWave(void *p = 0);
   static void *newArray_LHWave(Long_t size, void *p);
   static void delete_LHWave(void *p);
   static void deleteArray_LHWave(void *p);
   static void destruct_LHWave(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LHWave*)
   {
      ::LHWave *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LHWave >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LHWave", ::LHWave::Class_Version(), "LHEvent.h", 221,
                  typeid(::LHWave), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LHWave::Dictionary, isa_proxy, 4,
                  sizeof(::LHWave) );
      instance.SetNew(&new_LHWave);
      instance.SetNewArray(&newArray_LHWave);
      instance.SetDelete(&delete_LHWave);
      instance.SetDeleteArray(&deleteArray_LHWave);
      instance.SetDestructor(&destruct_LHWave);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LHWave*)
   {
      return GenerateInitInstanceLocal((::LHWave*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LHWave*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LHEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LHEvent::Class_Name()
{
   return "LHEvent";
}

//______________________________________________________________________________
const char *LHEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LHEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LHEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LHEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LHHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LHHit::Class_Name()
{
   return "LHHit";
}

//______________________________________________________________________________
const char *LHHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LHHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LHHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LHHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LHWave::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LHWave::Class_Name()
{
   return "LHWave";
}

//______________________________________________________________________________
const char *LHWave::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHWave*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LHWave::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHWave*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LHWave::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHWave*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LHWave::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHWave*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LHEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class LHEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LHEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(LHEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LHEvent(void *p) {
      return  p ? new(p) ::LHEvent : new ::LHEvent;
   }
   static void *newArray_LHEvent(Long_t nElements, void *p) {
      return p ? new(p) ::LHEvent[nElements] : new ::LHEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_LHEvent(void *p) {
      delete ((::LHEvent*)p);
   }
   static void deleteArray_LHEvent(void *p) {
      delete [] ((::LHEvent*)p);
   }
   static void destruct_LHEvent(void *p) {
      typedef ::LHEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LHEvent

//______________________________________________________________________________
void LHHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class LHHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LHHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(LHHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LHHit(void *p) {
      return  p ? new(p) ::LHHit : new ::LHHit;
   }
   static void *newArray_LHHit(Long_t nElements, void *p) {
      return p ? new(p) ::LHHit[nElements] : new ::LHHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_LHHit(void *p) {
      delete ((::LHHit*)p);
   }
   static void deleteArray_LHHit(void *p) {
      delete [] ((::LHHit*)p);
   }
   static void destruct_LHHit(void *p) {
      typedef ::LHHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LHHit

//______________________________________________________________________________
void LHWave::Streamer(TBuffer &R__b)
{
   // Stream an object of class LHWave.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LHWave::Class(),this);
   } else {
      R__b.WriteClassBuffer(LHWave::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LHWave(void *p) {
      return  p ? new(p) ::LHWave : new ::LHWave;
   }
   static void *newArray_LHWave(Long_t nElements, void *p) {
      return p ? new(p) ::LHWave[nElements] : new ::LHWave[nElements];
   }
   // Wrapper around operator delete
   static void delete_LHWave(void *p) {
      delete ((::LHWave*)p);
   }
   static void deleteArray_LHWave(void *p) {
      delete [] ((::LHWave*)p);
   }
   static void destruct_LHWave(void *p) {
      typedef ::LHWave current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LHWave

namespace {
  void TriggerDictionaryInitialization_LHEventDict_Impl() {
    static const char* headers[] = {
"LHEvent.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/lhaaso.ihep.ac.cn/anysw/slc5_ia64_gcc73/external/root/6.22.00/include/",
"/workfs2/ybj/chensz/software/DataRec/KM2Arec_3quarter_V2/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LHEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$LHEvent.h")))  LHEvent;
class __attribute__((annotate("$clingAutoload$LHEvent.h")))  LHHit;
class __attribute__((annotate("$clingAutoload$LHEvent.h")))  LHWave;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LHEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "LHEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"LHEvent", payloadCode, "@",
"LHHit", payloadCode, "@",
"LHWave", payloadCode, "@",
"StrDup", payloadCode, "@",
"operator+", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LHEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LHEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LHEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LHEventDict() {
  TriggerDictionaryInitialization_LHEventDict_Impl();
}
