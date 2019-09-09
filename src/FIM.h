#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace WGL_DG {

    template<int dim, class Mesh, class ProblemType>
    struct FIM {

        public:
            /// Typedefs
            template<class Scalar>
            using Vector = Eigen::Matrix<Scalar, dim, 1>;

            template<class Scalar>
            using Matirx = Eigen::Matirx<Scalar, dim, dim>;


        //======================================================================
            /// Constructors
            FIM(){

            }


        //======================================================================
            /// Primary Methods
            inline void compute(){
                while(activeList.size()){
                    for(int i=0; i<activeList.size(); i++){
                        
                    }
                }
            }

    }

}