#pragma once

#include "pch.h"

#include "TypeErasure.h"

using namespace Eigen;

namespace WGL_DG {

    // Type Erasure of a VectorFunction that is called with only Eigen::Ref to Dynamic Vectors and Matrices
    struct VectorFunction_EigenRefCall_Spec {
        struct Concept {
            virtual ~Concept() = default;
            // VectorFunction's interface:
            virtual void compute(const Ref<const VectorXd> & x, Ref<VectorXd> fx_) const = 0;
            virtual void jacobian(const Ref<const VectorXd> & x, Ref<MatrixXd> jx_) const = 0;
            virtual void compute_jacobian(const Ref<const VectorXd> & x, Ref<VectorXd> fx_, Ref<MatrixXd> jx_) const = 0;
            virtual void adjointgradient(const Ref<const VectorXd> & x, Ref<VectorXd> adjgrad_, const Ref<const VectorXd> & adjvars) const = 0;
        };

        template<class Holder>
        struct Model : public Holder, public virtual Concept {
            using Holder::Holder;
            virtual void compute(const Ref<const VectorXd> & x, Ref<VectorXd> fx_) const override {
                return rubber_types::model_get(this).compute(x, fx_);
            }
            virtual void jacobian(const Ref<const VectorXd> & x, Ref<MatrixXd> jx_) const override {
                return rubber_types::model_get(this).jacobian(x, jx_);
            }
            virtual void compute_jacobian(const Ref<const VectorXd> & x, Ref<VectorXd> fx_, Ref<MatrixXd> jx_) const override {
                return rubber_types::model_get(this).compute_jacobian(x, fx_, jx_);
            }
            virtual void adjointgradient(const Ref<const VectorXd> & x, Ref<VectorXd> adjgrad_, const Ref<const VectorXd> & adjvars) const override {
                return rubber_types::model_get(this).adjointgradient(x, adjgrad_, adjvars);
            }
        };

        template<class Container>
        struct ExternalInterface : public Container {
            using Container::Container;
            virtual void compute(const Ref<const VectorXd> & x, Ref<VectorXd> fx_) const {
                return rubber_types::interface_get(this).compute(x, fx_);
            }
            virtual void jacobian(const Ref<const VectorXd> & x, Ref<MatrixXd> jx_) const {
                return rubber_types::interface_get(this).jacobian(x, jx_);
            }
            virtual void compute_jacobian(const Ref<const VectorXd> & x, Ref<VectorXd> fx_, Ref<MatrixXd> jx_) const {
                return rubber_types::interface_get(this).compute_jacobian(x, fx_, jx_);
            }
            virtual void adjointgradient(const Ref<const VectorXd> & x, Ref<VectorXd> adjgrad_, const Ref<const VectorXd> & adjvars) const {
                return rubber_types::interface_get(this).adjointgradient(x, adjgrad_, adjvars);
            }
        };
    };

    using VectorFunction_EigenRefCall = rubber_types::TypeErasure<VectorFunction_EigenRefCall_Spec>;

////////////////////////////////////////////////////////////////////////////////

    // Type Erasure of a VectorFunction with some other calling convention

}