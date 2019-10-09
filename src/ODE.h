#pragma once

#include "pch.h"

#include "CRTPBase.h"

using namespace Eigen;

namespace WGL_DG {
    
    template<class Derived, int _XV, int _UV, int _PV>
    struct ODE : VectorFunction<Derived, _XV+1+_UV+_PV, _XV> {

        public:
            /// Typedefs


        //======================================================================
            /// Constructors
            ODE(){
                
            }


        //======================================================================
            /// Properties
            static const int XV = _XV;
            static const int UV = _UV;
            static const int PV = _PV;
            static const int XtUV = XV + 1 + UV;
            static const int XtUPV = XV + 1 + UV + PV;
            

    };

}