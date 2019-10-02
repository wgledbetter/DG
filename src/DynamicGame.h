#pragma once

#include <vector>
#include <Eigen/Core>

using namespace Eigen;

namespace WGL_DG {

    template<int _XV, int _UV, class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGame {

        public:
            /// Typedefs
            using Pursuer = PursuerDynamics;
            using Evader = EvaderDynamics;

            template<class Scalar>
            using PxVector = Matrix<Scalar, Pursuer::XV, 1>;
            template<class Scalar>
            using PuVector = Matrix<Scalar, Pursuer::UV, 1>;
            template<class Scalar>
            using PxuVector = Matrix<Scalar, Pursuer::XV+Pursuer::UV, 1>;
            template<class Scalar>
            using PJacMatrix = Matrix<Scalar, Pursuer::XV, Pursuer::XV+Pursuer::UV>;

            template<class Scalar>
            using ExVector = Matrix<Scalar, Evader::XV, 1>;
            template<class Scalar>
            using EuVector = Matrix<Scalar, Evader::UV, 1>;
            template<class Scalar>
            using ExuVector = Matrix<Scalar, Evader::XV+Evader::UV, 1>;
            template<class Scalar>
            using EJacMatrix = Matrix<Scalar, Evader::XV, Evader::XV+Evader::UV>;


        //======================================================================
            /// Constructors
            SeparableDynamicGame(Pursuer* p_set, Evader* e_set){
                set_pursuer(p_set);
                set_evader(e_set);
            }

            SeparableDynamicGame(){
                p = NULL;
                e = NULL;
            }


        //======================================================================
            /// Setup
            inline void set_pursuer(const Pursuer* p_set){
                p = p_set;
            }

            inline void gen_pursuer(){
                p = new Pursuer;
            }

            inline void set_evader(const Evader* e_set){
                e = e_set;
            }

            inline void gen_evader(){
                e = new Evader;
            }

        //______________________________________________________________________

            inline void set_pursuer_state_input(const Array<int, Pursuer::XV, 1> purVar){
                // Which part of the full state affects the pursuer's dynamics?
                pInStateIdx = purVar;
                pInStateControlIdx.template head<P_XV>() = purVar;
            }

            inline void set_pursuer_state_output(const Array<int, Pursuer::XV, 1> purVar){
                // Which state variables is the pursuer affecting?
                pOutStateIdx = purVar;
            }

            inline void set_pursuer_control_input(Array<int, Pursuer::UV, 1> purCon){
                // What part of the control vector is for the pursuer?
                pInControlIdx = purCon;
                pInStateControlIdx.template tail<P_UV>() = purCon + XV;
            }

        //______________________________________________________________________

            inline void set_evader_state_input(const Array<int, Evader::XV, 1> evaVar){
                // Which part of the full state affects the evader's dynamics?
                eInStateIdx = evaVar;
                eInStateControlIdx.template head<E_XV>() = evaVar;
            }

            inline void set_evader_state_output(const Array<int, Evader::XV, 1> evaVar){
                // Which state variables is the evader affecting?
                eOutStateIdx = evaVar;
            }

            inline void set_evader_control_input(const Array<int, Evader::UV, 1> evaCon){
                // What part of the control vector is for the evader?
                eInControlIdx = evaCon;
                eInStateControlIdx.template tail<E_UV>() = evaCon + XV;
            }

        //______________________________________________________________________

            inline void assume_full_default_coupling(){
                if( (P_XV == XV) && (E_XV == XV) ){
                    const Array<int, XV, 1> fullVars = Array<int, XV, 1>::LinSpaced(XV, 0, XV-1);
                    set_pursuer_state_input(fullVars);
                    set_pursuer_state_input(fullVars);
                    set_evader_state_input(fullVars);
                    set_evader_state_input(fullVars);
                }else{
                    // BAD
                }
            }

            inline void assume_decoupled(){
                if( P_XV + E_XV == XV ){
                    const Array<int, P_XV, 1> purVar = Array<int, P_XV, 1>::LinSpaced(P_XV, 0, P_XV-1);
                    set_pursuer_state_input(purVar);
                    set_pursuer_state_input(purVar);
                    const Array<int, P_UV, 1> purCon = Array<int, P_UV, 1>::LinSpaced(P_UV, 0, P_UV-1);
                    set_pursuer_control_input(purCon);

                    const Array<int, E_XV, 1> evaVar = Array<int, E_XV, 1>::LinSpaced(E_XV, P_XV, P_XV+E_XV-1);
                    set_evader_state_input(evaVar);
                    set_evader_state_output(evaVar);
                    const Array<int, E_UV, 1> evaCon = Array<int, E_UV, 1>::LinSpaced(E_UV, P_UV, P_UV+E_UV-1);
                    set_evader_control_input(evaCon);
                }else{
                    // BAD
                }
            }


        //======================================================================
            /// Primary Methods
            template<class InType, class OutType>
            inline void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {

                using Scalar = typename InType::Scalar;

                MatrixBase<OutType> & fx = fx_.const_cast_derived();

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);

                PxVector<Scalar> fxP;
                ExVector<Scalar> fxE;

                p->compute(xP, fxP);
                e->compute(xE, fxE);

                process_output_state(fxP, fxE, fx);
            }

        //______________________________________________________________________

            template<class InType, class JacType>
            inline void jacobian(const MatrixBase<InType> & x, MatrixBase<JacType> const & jx_) const {

                using Scalar = typename InType::Scalar;

                MatrixBase<JacType> & jx = jx_.const_cast_derived();

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);

                PJacMatrix<Scalar> jP;
                EJacMatrix<Scalar> jE;

                p->jacobian(xP, jP);
                e->jacobian(xE, jE);

                process_output_jacobian(jP, jE, jx);
            }

            template<class InType, class OutType, class JacType>
            inline void compute_jacobian(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_, MatrixBase<JacType> const & jx_) const {

                using Scalar = typename InType::Scalar;

                MatrixBase<OutType> fx = fx_.const_cast_derived();
                MatrixBase<JacType> jx = jx_.const_cast_derived();

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);

                PxVector<Scalar> fxP;
                ExVector<Scalar> fxE;
                PJacMatrix<Scalar> jP;
                EJacMatrix<Scalar> jE;

                p->compute(xP, fxP);
                p->jacobian(xP, jP);
                e->compute(xE, fxE);
                e->jacobian(xE, jE);

                process_output_state(fxP, fxE, fx);
                process_output_jacobian(jP, jE, jx);
            }

        //______________________________________________________________________

            inline void hessian(){

            }


        //======================================================================
            /// Access
            Pursuer* pointer_to_pursuer() const {
                return p;
            }

            Evader* pointer_to_evader() const {
                return e;
            }


        //======================================================================
            /// Properties
            static const int XV = _XV;
            static const int UV = _UV;
            static const int P_XV = Pursuer::XV;
            static const int P_UV = Pursuer::UV;
            static const int E_XV = Evader::XV;
            static const int E_UV = Evader::UV;


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            template<class InType, class POutType, class EOutType>
            inline void process_input_state(const MatrixBase<InType> & x, MatrixBase<POutType> const & px_, MatrixBase<EOutType> const & ex_) const {
                MatrixBase<POutType> & px = px_.const_cast_derived();
                MatrixBase<EOutType> & ex = ex_.const_cast_derived();

                px = x(pInStateControlIdx);
                ex = x(eInStateControlIdx);
            }

        //______________________________________________________________________
            template<class PInType, class EInType, class OutType>
            inline void process_output_state(const MatrixBase<PInType> & px, const MatrixBase<EInType> & ex, MatrixBase<OutType> const & x_) const {
                MatrixBase<OutType> x = x_.const_cast_derived();

                x.fill(0);
                x(pOutStateIdx) = px;
                x(eOutStateIdx) += ex;
            }

        //______________________________________________________________________
            template<class PJacType, class EJacType, class JacType>
            inline void process_output_jacobian(const MatrixBase<PJacType> & pj, const MatrixBase<EJacType> & ej, MatrixBase<JacType> const & jx_) const {
                MatrixBase<JacType> jx = jx_.const_cast_derived();

                jx.fill(0);
                jx(pOutStateIdx, pInStateControlIdx) = pj;
                jx(eOutStateIdx, eInStateControlIdx) += ej;
            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Methods

        //======================================================================
            /// Properties
            Pursuer* p;
            Evader* e;

            Array<int, P_XV, 1> pInStateIdx, pOutStateIdx;
            Array<int, P_UV, 1> pInControlIdx;
            Array<int, P_XV+P_UV, 1> pInStateControlIdx;
            Array<int, E_XV, 1> eInStateIdx, eOutStateIdx;
            Array<int, E_UV, 1> eInControlIdx;
            Array<int, E_XV+E_UV, 1> eInStateControlIdx;

    };

}