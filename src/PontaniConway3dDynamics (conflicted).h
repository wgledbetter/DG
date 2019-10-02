#pragma once

#include <vector>
#include <cmath>
#include <Eigen/Core>

#include "Dynamics.h"

using namespace Eigen;

namespace WGL_DG {

    struct PontaniConway3dDynamics : Dynamics<PontaniConway3dDynamics, 6, 2, 0> {

        public:
            /// Typedefs
            using Base = Dynamics;


        //======================================================================
            /// Constructors
            PontaniConway3dDynamics() {

            }


        //======================================================================
            /// Setup
            void set_mu(double m) {
                mu = m;
            }

            void set_thrust(double Th) {
                T = Th;
            }

            void set_mass(double ma) {
                m = ma;
            }


        //======================================================================
            /// Primary Methods
            template<class InType, class OutType>
            inline void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {

                using Scalar = typename InType::Scalar;
                using std::cos;
                using std::sin;

                Scalar r, v, gam, xi, phi, zeta, alph, beta;

                r = x[0];
                v = x[1];
                gam = x[2];
                xi = x[3];
                phi = x[4];
                zeta = x[5];
                alph = x[6];
                beta = x[7];

                Scalar sg, cg, sx, cx, sp, cp, sz, cz, sa, ca, sb, cb;
                sg = sin(gam);
                cg = cos(gam);
                sx = sin(xi);
                cx = cos(xi);
                sp = sin(phi);
                cp = cos(phi);
                sz = sin(zeta);
                cz = cos(zeta);
                sa = sin(alph);
                ca = cos(alph);
                sb = sin(beta);
                cb = cos(beta);

                MatrixBase<OutType> fx = fx_.const_cast_derived();

                fx[0] = v*sin(gam);
                fx[1] = T*cos(alph)*cos(beta)/m - mu*sin(gam)/(r*r);
                fx[2] = v*cos(gam)/r + T*sin(alph)*cos(beta)/(m*v) - mu*cos(gam)/(r*r*v);
                fx[3] = v*cos(gam)*cos(zeta)/(r*cos(phi));
                fx[4] = v*cos(gam)*sin(zeta)/r;
                fx[5] = T*sin(beta)/(m*v*cos(gam)) - v*cos(gam)*tan(phi)*cos(zeta)/r;

            }

        //______________________________________________________________________
            template<class InType, class JacType>
            inline void jacobian(const MatrixBase<InType> & x, MatrixBase<JacType> const & jx_) const {

                using Scalar = typename InType::Scalar;
                using std::cos;
                using std::sin;
                using std::tan;

                Scalar r, v, gam, xi, phi, zeta, alph, beta;

                r = x[0];
                v = x[1];
                gam = x[2];
                xi = x[3];
                phi = x[4];
                zeta = x[5];
                alph = x[6];
                beta = x[7];

                Scalar sg, cg, sx, cx, sp, cp, sz, cz, sa, ca, sb, cb;
                sg = sin(gam);
                cg = cos(gam);
                sx = sin(xi);
                cx = cos(xi);
                sp = sin(phi);
                cp = cos(phi);
                sz = sin(zeta);
                cz = cos(zeta);
                sa = sin(alph);
                ca = cos(alph);
                sb = sin(beta);
                cb = cos(beta);

                MatrixBase<JacType> jx = jx_.const_cast_derived();

                jx.fill(0);
                jx(0,1) = sg;
                jx(0,2) = v*cg;
                jx(1,0) = 2*mu*sg/(r*r*r);
                jx(1,2) = -mu*cg/(r*r);
                jx(1,6) = -T*sin(alph)*cos(beta)/m;
                jx(1,7) = -T*cos(alph)*sin(beta)/m;
                jx(2,0) = cg*(2*mu-r*v*v)/(r*r*r*v);
                jx(2,1) = -T*sin(alph)*cos(beta)/(m*v*v) + mu*cg/(r*r*v*v) + cg/r;
                jx(2,2) = sg*(mu-r*v*v)/(r*r*v);
                jx(2,6) = T*cos(alph)*cos(beta)/(m*v);
                jx(2,7) = -T*sin(alph)*sin(beta)/(m*v);
                jx(3,0) = -v*cg*cos(zeta)/(r*r*cos(phi));
                jx(3,1) = cg*cos(zeta)/(r*cos(phi));
                jx(3,2) = -v*sg*cos(zeta)/(r*cos(phi));
                jx(3,4) = v*cg*cos(zeta)*tan(phi)/(r*cos(phi));
                jx(3,5) = -v*cg*sin(zeta)/(r*cos(phi));
                jx(4,0) = -v*cg*sin(zeta)/(r*r);
                jx(4,1) = cg*sin(zeta)/r;
                jx(4,2) = -v*sg*sin(zeta)/r;
                jx(4,5) = v*cg*cos(zeta)/r;
                jx(5,0) = v*cg*tan(phi)*cos(zeta)/(r*r);
                jx(5,1) = -T*sin(beta)/(m*v*v*cg) - cg*tan(phi)*cos(zeta)/r;
                jx(5,2) = T*sin(beta)*tan(gam)/(m*v*cg) + v*sg*tan(phi)*cos(zeta)/r;
                jx(5,4) = -v*cg*cos(zeta)/(r*cos(phi)*cos(phi));
                jx(5,5) = v*cg*tan(phi)*sin(zeta)/r;
                jx(5,7) = T*cos(beta)/(m*v*cg);

            }

            template<class InType, class OutType, class JacType>
            inline void compute_jacobian(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_, MatrixBase<JacType> const & jx_) const {
                MatrixBase<OutType> fx = fx_.const_cast_derived();
                MatrixBase<JacType> jx = jx_.const_cast_derived();

                compute(x, fx);
                jacobian(x, jx);
            }

        //______________________________________________________________________
            template<class InType, class AdjVarType, class AdjGradType>
            inline void adjointgradient(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<AdjGradType> const & adjgrad_) const {



            }

        //______________________________________________________________________
            template<class InType, class AdjVarType, class OutType>
            inline void max_adjointgradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

            }

            template<class InType, class AdjVarType, class OutType>
            inline void min_adjointgradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            double mu;  // Gravitational Parameter (G*m)
            double T;  // Thrust
            double m;  // Spacecraft Mass


    };

}