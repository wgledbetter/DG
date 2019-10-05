#pragma once

#include "pch.h"

#include "ODE.h"

using namespace Eigen;

namespace WGL_DG {

    template<class SeparableGame>
    struct SemiDirect : ODE<SemiDirect, /* P_XV + E_XV + E_costate */, /* P_UV + E_UV */, /* P_PV + E_PV */> {

        public:
            /// Typedefs


        //======================================================================
            /// Constructors
            SemiDirect(){

            }


        //======================================================================
            /// Primary Methods
            template<class InType, class OutType>
            inline void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {

                using IScalar = typename InType::Scalar;



            }


        //======================================================================
            /// Properties



////////////////////////////////////////////////////////////////////////////////
        private:
            SeparableGame* game;

    };

}