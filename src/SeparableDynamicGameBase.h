#pragma once

#include "pch.h"

#include "ODE.h"
#include "CRTPBase.h"

using namespace Eigen;

namespace WGL_DG {

    template<class Derived, int _XV, int _PV, class PursuerDynamics, class EvaderDynamics>
    struct SeparableDynamicGameBase : ODE<Derived, _XV, PursuerDynamics::UV + EvaderDynamics::UV, _PV> {

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
            /// Properties
            static const int XV = _XV;
            static const int UV = Pursuer::UV + Evader::UV;
            static const int PV = _PV;
            static const int XtUV = XV + 1 + UV;
            static const int P_XV = Pursuer::XV;
            static const int P_UV = Pursuer::UV;
            static const int P_PV = Pursuer::PV;
            static const int E_XV = Evader::XV;
            static const int E_UV = Evader::UV;
            static const int E_PV = Evader::PV;

            Array<int, P_XV, 1> pInStateIdx, pOutStateIdx;
            Array<int, P_UV, 1> pInControlIdx;
            Array<int, P_XV+P_UV, 1> pInStateControlIdx;

            Array<int, E_XV, 1> eInStateIdx, eOutStateIdx;
            Array<int, E_UV, 1> eInControlIdx;
            Array<int, E_XV+E_UV, 1> eInStateControlIdx;


        //======================================================================
            /// Constructors
            SeparableDynamicGameBase(Pursuer* p_set, Evader* e_set){
                set_pursuer(p_set);
                set_evader(e_set);
            }

            SeparableDynamicGameBase(){
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
            
            inline void set_pursuer_state_input(const Array<int, P_XV, 1> purVar){
                // Which part of the full state affects the pursuer's dynamics?
                pInStateIdx = purVar;
                pInStateControlIdx.template head<P_XV>() = purVar;
            }

            inline void set_pursuer_state_output(const Array<int, P_XV, 1> purVar){
                // Which state variables is the pursuer affecting?
                pOutStateIdx = purVar;
            }

            inline void set_pursuer_control_input(Array<int, P_UV, 1> purCon){
                // What part of the control vector is for the pursuer?
                pInControlIdx = purCon;
                pInStateControlIdx.template tail<P_UV>() = purCon + P_XV;
            }

        //______________________________________________________________________

            inline void set_evader_state_input(const Array<int, E_XV, 1> evaVar){
                // Which part of the full state affects the evader's dynamics?
                eInStateIdx = evaVar;
                eInStateControlIdx.template head<E_XV>() = evaVar;
            }

            inline void set_evader_state_output(const Array<int, E_XV, 1> evaVar){
                // Which state variables is the evader affecting?
                eOutStateIdx = evaVar;
            }

            inline void set_evader_control_input(const Array<int, E_UV, 1> evaCon){
                // What part of the control vector is for the evader?
                eInControlIdx = evaCon;
                eInStateControlIdx.template tail<E_UV>() = evaCon + E_XV;
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
            template<class InType, class AdjGradType, class AdjVarType>
            inline void adjointgradient(const MatrixBase<InType> & x, MatrixBase<AdjGradType> const & adjgrad_, const MatrixBase<AdjVarType> & adjvars) const {

                using Scalar = typename AdjGradType::Scalar;

                Matrix<Scalar, XV, XV+UV> jac;
                this->jacobian(x, jac);

                MatrixBase<AdjGradType> & adjgrad = adjgrad_.const_cast_derived();
                adjgrad = jac*adjvars;

            }
            
        //______________________________________________________________________
            template<class InType, class AdjGradType, class AdjVarType>
            inline void adjointtransposegradient(const MatrixBase<InType> & x, MatrixBase<AdjGradType> const & adjgrad_, const MatrixBase<AdjVarType> & adjvars) const {

                using Scalar = typename AdjGradType::Scalar;

                Matrix<Scalar, XV, XV+UV> jac;
                this->jacobian(x, jac);

                MatrixBase<AdjGradType> & adjgrad = adjgrad_.const_cast_derived();
                adjgrad = jac.transpose()*adjvars;

            }


        //======================================================================
            /// Secondary Methods
            template<class InType, class AdjVarType, class OutType>
            inline void pursuer_max_adjointtransposegradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

                using Scalar = typename InType::Scalar;

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);
                
                auto adjP = adj(this->derived().pCostateIdx);
                Matrix<typename OutType::Scalar, P_UV, 1> uP;
                
                p->max_adjointtransposegradient_control(xP, adjP, u_);

            }
            
            template<class InType, class AdjVarType, class OutType>
            inline void pursuer_min_adjointtransposegradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

                using Scalar = typename InType::Scalar;

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);

                auto adjP = adj(this->derived().pCostateIdx);
                Matrix<typename OutType::Scalar, P_UV, 1> uP;

                p->min_adjointtransposegradient_control(xP, adjP, u_);

            }
            
        //______________________________________________________________________
            
            template<class InType, class AdjVarType, class OutType>
            inline void evader_max_adjointtransposegradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

                using Scalar = typename InType::Scalar;

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);

                auto adjE = adj(this->derived().eCostateIdx);
                Matrix<typename OutType::Scalar, E_UV, 1> uE;

                e->max_adjointtransposegradient_control(xE, adjE, u_);

            }
            
            template<class InType, class AdjVarType, class OutType>
            inline void evader_min_adjointtransposegradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

                using Scalar = typename InType::Scalar;

                PxuVector<Scalar> xP;
                ExuVector<Scalar> xE;

                process_input_state(x, xP, xE);

                auto adjE = adj(this->derived().eCostateIdx);
                Matrix<typename OutType::Scalar, E_UV, 1> uE;

                e->min_adjointtransposegradient_control(xE, adjE, u_);

            }
            
            
        //======================================================================
            /// Access
            Pursuer* pointer_to_pursuer() const {
                return p;
            }

            Evader* pointer_to_evader() const {
                return e;
            }
            
        //______________________________________________________________________


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            template<class InType, class POutType, class EOutType>
            inline void process_input_state(const MatrixBase<InType> & x, MatrixBase<POutType> const & px_, MatrixBase<EOutType> const & ex_) const {
                MatrixBase<POutType> & px = px_.const_cast_derived();
                MatrixBase<EOutType> & ex = ex_.const_cast_derived();

                // This type of array-indexing requires Eigen>3.3
                px = x(pInStateControlIdx);
                ex = x(eInStateControlIdx);
            }

        //______________________________________________________________________
            template<class PInType, class EInType, class OutType>
            inline void process_output_state(const MatrixBase<PInType> & px, const MatrixBase<EInType> & ex, MatrixBase<OutType> const & x_) const {
                MatrixBase<OutType> x = x_.const_cast_derived();

                x.fill(0);
                // This type of array-indexing requires Eigen>3.3
                x(pOutStateIdx) = px;
                x(eOutStateIdx) += ex;
            }

        //______________________________________________________________________
            template<class PJacType, class EJacType, class JacType>
            inline void process_output_jacobian(const MatrixBase<PJacType> & pj, const MatrixBase<EJacType> & ej, MatrixBase<JacType> const & jx_) const {
                MatrixBase<JacType> jx = jx_.const_cast_derived();

                jx.fill(0);
                // This type of array-indexing requires Eigen>3.3
                jx(pOutStateIdx, pInStateControlIdx) = pj;
                jx(eOutStateIdx, eInStateControlIdx) += ej;
            }
            
            
        //======================================================================
            /// Properties
            Pursuer* p;
            Evader* e;


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Methods

        //======================================================================
            /// Properties

    };

}