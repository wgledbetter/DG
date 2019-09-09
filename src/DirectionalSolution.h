#pragma once

#include "Mesh.h"

namespace WGL_DG {
    
    template<int dim, class VectorFunction, class Mesh>
    struct DirectionalSolution : Mesh {

        public:
            /// Typedefs
            template<class Scalar>
            using Vector = Eigen::Matrix<Scalar, dim, 1>;

            template<class Scalar>
            using Matrix = Eigen::Matrix<Scalar, dim, dim>;


        //======================================================================
            /// Constructors
            DirectionalSolution(){
                fim.convTol = 1e-6;
                threads = thread::hardware_concurrency();

                
            }

    }

}