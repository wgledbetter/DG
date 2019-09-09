#pragma once

#include <Eigen/Core>

namespace WGL_DG {

    template<int dim>
    struct VectorFunction {

        /// Typedefs
        template<class Scalar>
        using Vector = Eigen::Matrix<Scalar, dim, 1>;

    };

////////////////////////////////////////////////////////////////////////////////

    template<int dim>
    struct ParabolicSinkVectorFunction : VectorFunction<dim> {

        /// Typedefs
        using Base = VectorFunction<dim>;

        template<class Scalar>
        using Vector = typename Base::template Vector<Scalar>;

        /// Primary Methods
        inline Vector<double> compute(Vector<double> x) {
            return -x * x.norm();
        }

    };

}