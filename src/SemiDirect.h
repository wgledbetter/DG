#pragma once

#include "pch.h"

#include "ODE.h"

using namespace Eigen;

namespace WGL_DG {

    template<class SeparableGame>
    struct SemiDirect : ODE<SemiDirect<SeparableGame>, ((SeparableGame::nECostate == -1) ? -1 : SeparableGame::XV + SeparableGame::nECostate), SeparableGame::P_UV, SeparableGame::PV> {

        public:
            /// Typedefs


        //======================================================================
            /// Constructors
            SemiDirect(){

            }


        //======================================================================
            /// Primary Methods
            template<class InType, class OutType>
            inline void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {

                using IScalar = typename InType::Scalar;

                Matrix<typename OutType::Scalar, SeparableGame::E_UV, 1> Econtrol;
                Matrix<IScalar, SD_XtUPV, 1> state;
                Matrix<IScalar, G_XV, 1> stateDot;
                Matrix<IScalar, G_CV, 1> costate;
                Matrix<IScalar, G_CV, 1> costateDot;
                Matrix<IScalar, SeparableGame::XtUPV, 1> fullState = Matrix<IScalar, SeparableGame::XtUPV, 1>::Zero();

                process_input_state(x, state, costate);

                // Compute the game's dynamics
                    // Given the input [state, Ecostate], calculate the evader's control
                SeparableGame::evader_max_adjointtransposegradient_control(state, costate, Econtrol);
                fullState.template head<G_XV+1>() = state.template head<G_XV+1>();
                fullState(game->pInControlIdx) = state.segment<SeparableGame::P_UV>(G_XV+1);
                fullState(game->eInControlIdx) += Econtrol;
                fullState.template tail<SeparableGame::PV>() = state.template tail<SeparableGame::PV>();

                game->compute(fullState, stateDot);


                // Compute the game's evader's costate dynamics
                game->adjointtransposegradient(fullState, costate, costateDot);


                // Output
                MatrixBase<OutType> & fx = fx_.const_cast_derived();
                process_output_state(stateDot, costateDot, fx);

            }

        //______________________________________________________________________

            template<class InType, class JacType>
            inline void jacobian(const MatrixBase<InType> & x, MatrixBase<JacType> const & jx_) const {

            }


        //======================================================================
            /// Properties
            static const int G_XV = SeparableGame::XV;
            static const int G_CV = SeparableGame::nECostate;
            static const int SD_XtUV = G_XV + 1 + SeparableGame::P_UV;
            static const int SD_XtUPV = SD_XtUV + SeparableGame::PV;



////////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            template<class InType, class StateOutType, class CostateOutType>
            inline void process_input_state(const MatrixBase<InType> & x, MatrixBase<StateOutType> const & st_, MatrixBase<CostateOutType> const & cs_) const {
                // x = [state Ecostate time Pcontrol params]
                MatrixBase<StateOutType> & st = st_.const_cast_derived();
                MatrixBase<CostateOutType> & cs = cs_.const_cast_derived();

                cs = x.segment<G_CV>(G_XV);

                st.template head<G_XV>() = x.template head<G_XV>();
                st.template tail<SD_XtUPV - G_XV>() = x.template tail<SD_XtUPV - G_XV>();
            }

            template<class StateInType, class CostateInType, class OutType>
            inline void process_output_state(const MatrixBase<StateInType> & st, const MatrixBase<CostateInType> & cs, MatrixBase<OutType> const & x_) const {
                MatrixBase<OutType> & x = x_.const_cast_derived();

                x.template head<G_XV>() = st.template head<G_XV>();
                x.segment<G_CV>(G_XV) = cs;
            }



////////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            SeparableGame* game;

    };

}