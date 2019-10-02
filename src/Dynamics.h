#pragma once

#include <vector>
#include <Eigen/Core>

#include "CRTPBase.h"

using namespace Eigen;

namespace WGL_DG {
    
    template<class Derived, int _XV, int _UV, int _PV>
    struct Dynamics : CRTPBase<Derived> {

        public:
            /// Typedefs


        //======================================================================
            /// Constructors


        //======================================================================
            /// Properties
            static const int XV = _XV;
            static const int UV = _UV;
            static const int PV = _PV;


        //======================================================================
            /// Primary Methods
            template<class InType, class OutType>
            void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {
                this->derived().compute(x, fx_);
            }

        //______________________________________________________________________
            template<class InType, class JacType>
            void jacobian(const MatrixBase<InType> & x, MatrixBase<JacType> const & jx_) const {
                this->derived().jacobian(x, jx_);
            }

            template<class InType, class OutType, class JacType>
            void compute_jacobian(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_, MatrixBase<JacType> const & jx_) const {
                this->derived().compute_jacobian(x, fx_, jx_);
            }

        //______________________________________________________________________
            

    };

}