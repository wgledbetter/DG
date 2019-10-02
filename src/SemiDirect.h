#pragma once

#include <Eigen/Core>

#include "DynamicGame.h"

namespace WGL_DG {

    template<class SeparableGame, class PursuerTranscription, class EvaderTranscription, class CostateTranscription=LGL7>
    struct SemiDirect : SeparableGame, PursuerTranscription, EvaderTranscription, CostateTranscription {

        public:
            /// Typedefs


        //======================================================================
            /// Constructors
            SemiDirect(){

            }


        //======================================================================
            /// Primary Methods


////////////////////////////////////////////////////////////////////////////////
        private:
            SeparableGame* game;

    };

}