#pragma once

#include "pch.h"

#include "SeparableDynamicGameBase.h"

using namespace Eigen;

namespace WGL_DG {

    enum {
        Decoupled,
        FullyCoupled,
        PartiallyCoupled,
    };

    /// Template
    template<int _PV, class PursuerDynamics, class EvaderDynamics, int Coupling=Decoupled>
    struct SeparableDynamicGame : SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, Coupling>, -1, _PV, PursuerDynamics, EvaderDynamics> {

    };

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //// Specializations

    /// Decoupled
    template<int _PV, class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, Decoupled> : 
        SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, Decoupled>, (PursuerDynamics::XV + EvaderDynamics::XV), _PV, PursuerDynamics, EvaderDynamics> {

        public:
            /// Typedefs
            using GameBase = SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, Decoupled>, (PursuerDynamics::XV + EvaderDynamics::XV), _PV, PursuerDynamics, EvaderDynamics>;

            using typename GameBase::Pursuer;
            using typename GameBase::Evader;

            using GameBase::XV;
            using GameBase::UV;
            using GameBase::PV;
            using GameBase::P_XV;
            using GameBase::P_UV;
            using GameBase::P_PV;
            using GameBase::E_XV;
            using GameBase::E_UV;
            using GameBase::E_PV;


        //======================================================================
            /// Properties
            static const int nPCostate = P_XV;
            static const int nECostate = E_XV;

            Array<int, nPCostate, 1> pCostateIdx;
            Array<int, nECostate, 1> eCostateIdx;


        //======================================================================
            /// Constructors
            SeparableDynamicGame() : GameBase::SeparableDynamicGameBase() {
                constructorFunction();
            }

            SeparableDynamicGame(Pursuer* p_set, Evader* e_set) : GameBase::SeparableDynamicGameBase(p_set, e_set) {
                constructorFunction();
            }


        //======================================================================


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            inline void constructorFunction(){
                int pxv = P_XV;
                int puv = P_UV;
                int exv = E_XV;

                const Array<int, P_XV, 1> purVar = Array<int, P_XV, 1>::LinSpaced(P_XV, 0, P_XV-1);
                const Array<int, P_UV, 1> purCon = Array<int, P_UV, 1>::LinSpaced(P_UV, 0, P_UV-1);
                GameBase::pInStateIdx = purVar;
                GameBase::pOutStateIdx = purVar;
                GameBase::pInControlIdx = purCon;
                GameBase::pInStateControlIdx.template head<P_XV>() = purVar;
                GameBase::pInStateControlIdx.template tail<P_UV>() = purCon + pxv;
                pCostateIdx = purVar;

                const Array<int, E_XV, 1> evaVar = Array<int, E_XV, 1>::LinSpaced(E_XV, pxv, P_XV+E_XV-1);
                const Array<int, E_UV, 1> evaCon = Array<int, E_UV, 1>::LinSpaced(E_UV, puv, P_UV+E_UV-1);
                GameBase::eInStateIdx = evaVar;
                GameBase::eOutStateIdx = evaVar;
                GameBase::eInControlIdx = evaCon;
                GameBase::eInStateControlIdx.template head<E_XV>() = evaVar;
                GameBase::eInStateControlIdx.template tail<E_UV>() = evaCon + exv;
                eCostateIdx = evaVar;
            }


        //======================================================================
            /// Properties



    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties


    };
    
    
//******************************************************************************
//******************************************************************************
    /// Fully Coupled
    template<int _PV, class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, FullyCoupled> : 
        SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, FullyCoupled>, PursuerDynamics::XV, _PV, PursuerDynamics, EvaderDynamics> {

        public:
            /// Typedefs
            using GameBase = SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, FullyCoupled>, PursuerDynamics::XV, _PV, PursuerDynamics, EvaderDynamics>;

            using GameBase::XV;
            using GameBase::UV;
            using GameBase::PV;
            using GameBase::P_XV;
            using GameBase::P_UV;
            using GameBase::P_PV;
            using GameBase::E_XV;
            using GameBase::E_UV;
            using GameBase::E_PV;


        //======================================================================
            /// Properties
            static const int nPCostate = XV;
            static const int nECostate = XV;

            Array<int, nPCostate, 1> pCostateIdx;
            Array<int, nECostate, 1> eCostateIdx;


        //======================================================================
            /// Constructors
            SeparableDynamicGame() : GameBase::SeparableDynamicGameBase() {
                constructorFunction();
            }

            SeparableDynamicGame(typename GameBase::Pursuer* p_set, typename GameBase::Evader* e_set) : GameBase::SeparableDynamicGameBase(p_set, e_set) {
                constructorFunction();
            }
            

        //======================================================================


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            inline void constructorFunction(){
                if( (P_XV == XV) && (E_XV == XV) ){
                    const Array<int, XV, 1> fullVars = Array<int, XV, 1>::LinSpaced(XV, 0, XV-1);
                    GameBase::set_pursuer_state_input(fullVars);
                    GameBase::set_pursuer_state_input(fullVars);
                    GameBase::set_evader_state_input(fullVars);
                    GameBase::set_evader_state_input(fullVars);
                    pCostateIdx = fullVars;
                    eCostateIdx = fullVars;
                }else{
                    // BAD
                }
            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties


    };


//******************************************************************************
//******************************************************************************
    /// Partially Coupled
    template<int _PV, class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, PartiallyCoupled> : 
        SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, PartiallyCoupled>, -1, _PV, PursuerDynamics, EvaderDynamics> {

        public:
            /// Typedefs
            using GameBase = SeparableDynamicGameBase<SeparableDynamicGame<_PV, PursuerDynamics, EvaderDynamics, PartiallyCoupled>, -1, _PV, PursuerDynamics, EvaderDynamics>;

            using GameBase::XV;
            using GameBase::UV;
            using GameBase::PV;
            using GameBase::P_XV;
            using GameBase::P_UV;
            using GameBase::P_PV;
            using GameBase::E_XV;
            using GameBase::E_UV;
            using GameBase::E_PV;


        //======================================================================
            /// Properties
            static const int nPCostate = -1;
            static const int nECostate = -1;

            Array<int, nPCostate, 1> pCostateIdx;
            Array<int, nECostate, 1> eCostateIdx;


        //======================================================================
            /// Constructors
            SeparableDynamicGame() : GameBase::SeparableDynamicGameBase() {
                constructorFunction();
            }

            SeparableDynamicGame(typename GameBase::Pursuer* p_set, typename GameBase::Evader* e_set) : GameBase::SeparableDynamicGameBase(p_set, e_set) {
                constructorFunction();
            }


        //======================================================================
            /// Setup
            inline void set_XVars(int x){

            }



    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            inline void constructorFunction(){
                // Mixed/Partial Coupling
            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties


    };

}