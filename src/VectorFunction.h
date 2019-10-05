#pragma once

#include "pch.h"

#include "CRTPBase.h"

namespace WGL_DG {

    template<class Derived, int _IR, int _OR>
    struct VectorFunction : CRTPBase<Derived> {

        /// Properties
            static const int IR = _IR;
            static const int OR = _OR;


    //==========================================================================
        /// Typedefs
            template<class Scalar>
            using InputVec = Eigen::Matrix<Scalar, IR, 1>;
            template<class Scalar>
            using OutputVec = Eigen::Matrix<Scalar, OR, 1>;


    //==========================================================================
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


    //==========================================================================
        /// Other Methods


    };

}