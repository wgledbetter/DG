#pragma once

#include "pch.h"

#include "VectorFunction.h"

namespace WGL_DG {

    template<int dim>
    struct ParabolicSinkVectorFunction : VectorFunction<ParabolicSinkVectorFunction<dim>, dim, dim> {

        /// Typedefs
        using Base = VectorFunction<ParabolicSinkVectorFunction<dim>, dim, dim>;

        template<class Scalar>
        using Vector = typename Base::template Vector<Scalar>;

        /// Primary Methods
        inline Vector<double> compute(Vector<double> x) {
            return -x * x.norm();
        }

    };

}