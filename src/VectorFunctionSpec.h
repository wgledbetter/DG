#pragma once

#include "pch.h"

#include "VectorFunction.h"

namespace WGL_DG {

    template<int dim>
    struct ParabolicSinkVectorFunction : VectorFunction<ParabolicSinkVectorFunction<dim>, dim, dim> {

        /// Typedefs
        using Base = VectorFunction<ParabolicSinkVectorFunction<dim>, dim, dim>;

        template<class Scalar>
        using InputVec = typename Base::template InputVec<Scalar>;
        template<class Scalar>
        using OutputVec = typename Base::template OutputVec<Scalar>;

        /// Primary Methods
        inline OutputVec<double> compute(InputVec<double> x) {
            return -x * x.norm();
        }

    };

}