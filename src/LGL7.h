#pragma once

#include "pch.h"

#include "VectorFunction.h"

namespace WGL_DG {
    
    template<class ODE>
    struct LGL7 {

        public:
            /// Properties
            static const int inputDim = 2*ODE::XtUVars;


        //======================================================================
            /// Typedefs


        //======================================================================
            /// Constructors
            LGL7(){

            }


        //======================================================================
            /// Defect Constraint
            struct Defect : VectorFunction<Defect, 2*ODE::XtUVars, 3>

    };

}