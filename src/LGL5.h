#pragma once

#include "pch.h"

#include "VectorFunction.h"

namespace WGL_DG {
    
    template<class ODE>
    struct LGL5 {

        public:
            /// Properties
            static const int inputDim = 2*ODE::XtUV;
            static const int nDfct = 2;


        //======================================================================
            /// Typedefs


        //======================================================================
            /// Constructors
            LGL5(){

            }


        //======================================================================
            /// Defect Constraint
            struct Defect : VectorFunction<Defect, inputDim, nDfct> {

                public:
                    /// Typedefs
                    using Base = VectorFunction<Defect, inputDim, nDfct>;
                    template<class Scalar>
                    using ODEInputVec = ODE::template InputVec<Scalar>;
                    template<class Scalar>
                    using ODEOutputVec = ODE::template OutputVec<Scalar>;


                //==============================================================
                    /// Constructors
                    Defect(const ODE & od){
                        this->ode = od;
                    }


                //==============================================================
                    /// Primary Methods
                    template<class InType, class OutType>
                    inline void compute(const Eigen::MatrixBase<InType> & x, Eigen::MatrixBase<OutType> const & fx_) const {
                        // Expect x = [x1 t1 u1 x2 t2 u2 p]

                        using InScalar = typename InType::Scalar;
                        using OutScalar = typename OutType::Scalar;

                        // Extract left X
                        ODEInputVec<InScalar> X_K(ODE::IR);
                        X_K.template head<ODE::XtUV>(ODE::XtUV) = x.template head<ODE::XtUV>(ODE::XtUV);
                        X_K.template tail<ODE::PV>(ODE::PV) = x.template tail<ODE::PV>(ODE::PV);

                        // Extract right X
                        ODEInputVec<InScalar> X_Kp1(ODE::IR);
                        X_Kp1.template head<ODE::XtUV>(ODE::XtUV) = x.template segment<ODE::XtUV>(ODE::XtUV, ODE::XtUV);
                        X_Kp1.template tail<ODE::PV>(ODE::PV) = x.template tail<ODE::PV>(ODE::PV);

                        // Calculate time difference
                        InScalar h = X_Kp1[ODE::XV + 1] - X_K[ODE::XV + 1];

                        // Compute F(x_K)
                        ODEOutputVec<InScalar> F_K(ODE::XV);
                        ode.compute(X_K, F_K);

                        // Compute F(x_Kp1)
                        ODEOutputVec<InScalar> F_Kp1(ODE::XV);
                        ode.compute(X_Kp1, F_Kp1)

                        // Calculate central point
                        ODEInputVec<InScalar> X_C(ODE::IR);
                        X_C = (X_K + X_Kp1)/2;
                        X_C.template head<ODE::XV>(ODE::XV) = 

                    }



            ////////////////////////////////////////////////////////////////////
                protected:
                    /// Properties
                    ODE ode;

            };

    };

}