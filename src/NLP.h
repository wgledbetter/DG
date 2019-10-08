#pragma once

#include "pch.h"

#include "VectorFunction.h"
#include "VectorFunctionTypeErasure.h"
#include "ctpl.h"

using namespace Eigen;

namespace WGL_DG {

    template<int IR>
    struct NLP {

        // Objective function
        struct Objective : VectorFunction<Objective, IR, 1> {

            

        };

        //----------------------------------------------------------------------

        // Constraint function
        struct Constraints : VectorFunction<Constraints, IR, -1> {

            public:
                /// Typedefs
                
                
            //==================================================================
                /// Properties
                
                
            //==================================================================
                /// Constructors
                Constraints(){
                    nConstr = 0;
                    computeIdx.push_back(0);
                }


            //==================================================================
                /// Setup
                inline void addConstr(VectorFunction_EigenRefCall* pfun){
                    this->funx.push_back(pfun);
                    this->computeIdx.push_back( computeIdx.back() + pfun->OR() );
                }


            //==================================================================
                /// Primary Methods
                template<class InType, class OutType>
                inline void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {
                    
                    MatrixBase<OutType> & fx = fx_.const_cast_derived();
                    
                    ConstrThreads.restart();
                    for(int i=0; i<nConstr; i++){
                        ConstrThreads.push( [this, x, fx](int thr){ this->funx[i].compute(x, fx_.segment(computeIdx[i], funx[i]->getOR())); } );
                    }
                    ConstrThreads.stop();

                }
                
                
        ////////////////////////////////////////////////////////////////////////
            protected:
                /// Properties
                int nConstr;

                std::vector<int> computeIdx;  // RENAME THIS
                std::vector<VectorFunction_EigenRefCall*> funx;

                ctpl::thread_pool ConstrThreads;

        };
        
        Constraints equalityConstraints;
        Constraints inequalityConstraints;

    };

}