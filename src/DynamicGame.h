#pragma once

#include <vector>
#include <Eigen/Core>

namespace WGL_DG {

    enum {
        Coupled = true,
        Decoupled = false,
    }

    template<class PursuerDynamics, class EvaderDynamics, bool Coupling>
    struct SeparableDynamicGame : PursuerDynamics {

    };

//------------------------------------------------------------------------------

    template<class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGame<PursuerDynamics, EvaderDynamics, Coupled> :
        PursuerDynamics, EvaderDynamics {

        public:
            /// Typedefs
            using Pursuer = PursuerDynamics;
            using Evader = EvaderDynamics;


        //======================================================================
            /// Constructors


        //======================================================================
            /// Setup


        //======================================================================
            /// Primary Methods
            inline void compute()


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Methods

        //======================================================================
            /// Properties

    };

//------------------------------------------------------------------------------

    template<class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGame<PursuerDynamics, EvaderDynamics, Decoupled> :
        PursuerDynamics, EvaderDynamics {

        public:
            /// Typedefs
            using Pursuer = PursuerDynamics;
            using Evader = EvaderDynamics;


        //======================================================================
            /// Constructors


        //======================================================================
            /// Setup


        //======================================================================
            /// Primary Methods
            inline void compute()


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Methods

        //======================================================================
            /// Properties

    };

}