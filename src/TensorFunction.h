#pragma once

#include <Eigen/Core>

namespace WGL_DG {

    template<int dim>
    struct TensorFunction {
        
        /// Typedefs
        template<class Scalar>
        using Vector = Eigen::Matrix<Scalar, dim, 1>;
        
        template<class Scalar>
        using Matrix = Eigen::Matrix<Scalar, dim, dim>;

    };
    
////////////////////////////////////////////////////////////////////////////////
    
    template<int dim>
    struct IdentityTensorFunction : TensorFunction<dim> {
        
        /// Typedefs
        using Base = TensorFunction<dim>;
        
        template<class Scalar>
        using Matrix = typename Base::template Matrix<Scalar>;
        
        template<class Scalar>
        using Vector = typename Base::template Vector<Scalar>;
        
        /// Primary Methods
        inline Matrix<double> compute(Vector<double> x){
            return Matrix<double>::Identity();
        }
        
    };
    
////////////////////////////////////////////////////////////////////////////////
    
    template<int dim>
    struct FastEdgesTensorFunction : TensorFunction<dim> {
        
        /// Typedefs
        using Base = TensorFunction<dim>;
        
        template<class Scalar>
        using Matrix = typename Base::template Matrix<Scalar>;
        
        template<class Scalar>
        using Vector = typename Base::template Vector<Scalar>;
        
        /// Primary Methods
        inline Matrix<double> compute(Vector<double> x){
            return (std::pow(2*x.norm(), 4) + 1) * Matrix<double>::Identity();
        }
        
    };

}