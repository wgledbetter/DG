#pragma once

#include "pch.h"

#include "ODE.h"

using namespace Eigen;

namespace WGL_DG {

    struct PontaniConway3dDynamics : ODE<PontaniConway3dDynamics, 6, 2, 0> {

        public:
            /// Typedefs
            using Base = ODE<PontaniConway3dDynamics, 6, 2, 0>;


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
                alph = x[7];
                beta = x[8];

                Scalar sg, cg, cp, sz, cz, sa, ca, sb, cb;
                sg = sin(gam);
                cg = cos(gam);
                cp = cos(phi);
                sz = sin(zeta);
                cz = cos(zeta);
                sa = sin(alph);
                ca = cos(alph);
                sb = sin(beta);
                cb = cos(beta);

                MatrixBase<OutType> & fx = fx_.const_cast_derived();

                fx[0] = v*sg;
                fx[1] = T*ca*cb/m - mu*sg/(r*r);
                fx[2] = v*cg/r + T*sa*cb/(m*v) - mu*cg/(r*r*v);
                fx[3] = v*cg*cz/(r*cp);
                fx[4] = v*cg*sz/r;
                fx[5] = T*sb/(m*v*cg) - v*cg*tan(phi)*cz/r;

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
                alph = x[7];
                beta = x[8];

                Scalar sg, cg, cp, sz, cz, sa, ca, sb, cb;
                sg = sin(gam);
                cg = cos(gam);
                cp = cos(phi);
                sz = sin(zeta);
                cz = cos(zeta);
                sa = sin(alph);
                ca = cos(alph);
                sb = sin(beta);
                cb = cos(beta);

                MatrixBase<JacType> & jx = jx_.const_cast_derived();

                jx.fill(0);
                jx(0,1) = sg;
                jx(0,2) = v*cg;
                jx(1,0) = 2*mu*sg/(r*r*r);
                jx(1,2) = -mu*cg/(r*r);
                jx(1,6) = -T*sa*cb/m;
                jx(1,7) = -T*ca*sb/m;
                jx(2,0) = cg*(2*mu-r*v*v)/(r*r*r*v);
                jx(2,1) = -T*sa*cb/(m*v*v) + mu*cg/(r*r*v*v) + cg/r;
                jx(2,2) = sg*(mu-r*v*v)/(r*r*v);
                jx(2,6) = T*ca*cb/(m*v);
                jx(2,7) = -T*sa*sb/(m*v);
                jx(3,0) = -v*cg*cz/(r*r*cp);
                jx(3,1) = cg*cz/(r*cp);
                jx(3,2) = -v*sg*cz/(r*cp);
                jx(3,4) = v*cg*cz*tan(phi)/(r*cp);
                jx(3,5) = -v*cg*sz/(r*cp);
                jx(4,0) = -v*cg*sz/(r*r);
                jx(4,1) = cg*sz/r;
                jx(4,2) = -v*sg*sz/r;
                jx(4,5) = v*cg*cz/r;
                jx(5,0) = v*cg*tan(phi)*cz/(r*r);
                jx(5,1) = -T*sb/(m*v*v*cg) - cg*tan(phi)*cz/r;
                jx(5,2) = T*sb*tan(gam)/(m*v*cg) + v*sg*tan(phi)*cz/r;
                jx(5,4) = -v*cg*cz/(r*cp*cp);
                jx(5,5) = v*cg*tan(phi)*sz/r;
                jx(5,7) = T*cb/(m*v*cg);

            }

            template<class InType, class OutType, class JacType>
            inline void compute_jacobian(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_, MatrixBase<JacType> const & jx_) const {
                compute(x, fx_);
                jacobian(x, jx_);
            }

        //______________________________________________________________________
            template<class InType, class AdjVarType, class AdjGradType>
            inline void adjointgradient(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<AdjGradType> const & adjgrad_) const {

                using Scalar = typename InType::Scalar;

                Matrix<Scalar, 6, 8> jac;
                jacobian(x,jac);

                MatrixBase<AdjGradType> & adjgrad = adjgrad_.const_cast_derived();
                adjgrad[0] = (jac(0,1)*adj[1]) + (jac(0,2)*adj[2]);
                adjgrad[1] = (jac(1,0)*adj[0]) + (jac(1,2)*adj[2]);
                adjgrad[2] = (jac(2,0)*adj[0]) + (jac(2,1)*adj[1]) + (jac(2,2)*adj[2]);
                adjgrad[3] = (jac(3,0)*adj[0]) + (jac(3,1)*adj[1]) + (jac(3,2)*adj[2]) + (jac(3,4)*adj[4]) + (jac(3,5)*adj[5]);
                adjgrad[4] = (jac(4,0)*adj[0]) + (jac(4,1)*adj[1]) + (jac(4,2)*adj[2]) + (jac(4,5)*adj[5]);
                adjgrad[5] = (jac(5,0)*adj[0]) + (jac(5,1)*adj[1]) + (jac(5,2)*adj[2]) + (jac(5,4)*adj[4]) + (jac(5,5)*adj[5]);

            }

            template<class InType, class AdjVarType, class AdjGradType>
            inline void adjointtransposegradient(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<AdjGradType> const & adjgrad_) const {

                using Scalar = typename InType::Scalar;
                using std::cos;
                using std::sin;

                Scalar x1, x2, x3, x4, x5, x6, u1, u2;
                x1 = x[0];
                x2 = x[1];
                x3 = x[2];
                x4 = x[3];
                x5 = x[4];
                x6 = x[5];
                u1 = x[7];
                u2 = x[8];

                Scalar sx3, cx3, sx5, cx5, sx6, cx6, su1, cu1, su2, cu2;
                sx3 = sin(x3);
                cx3 = cos(x3);
                sx5 = sin(x5);
                cx5 = cos(x5);
                sx6 = sin(x6);
                su1 = sin(x1);
                cu1 = cos(u1);
                su2 = sin(u2);
                cu2 = cos(u2);

                typename AdjVarType::Scalar l1, l2, l3, l4, l5, l6;
                l1 = adj[0];
                l2 = adj[1];
                l3 = adj[2];
                l4 = adj[3];
                l5 = adj[4];
                l6 = adj[5];

                MatrixBase<AdjGradType> & adjgrad = adjgrad_.const_cast_derived();
                adjgrad[0] = (x2*l3*cx3 + l4*x2*cx3*cx6/cx5 + x2*l5*cx3*sx6 - x2*l6*cx3*sx5*cx6/cx5)/(x1*x1) - (2*mu*l2*sx3 + 2*mu*l3*cx3/x2)/(x1*x1*x1);
                adjgrad[1] = -l1*sx3 = l3*(cx3*(1/x1 + mu/(x1*x1*x2*x2)) - T*su1*cu2/(m*x2*x2)) - l4*cx3*cx6/(x1*cx5) - l5*cx3*sx6/x1 + l6*(T*su2/(m*x2*x2*cx3) + cx3*sx5*cx6/(x1*cx5));
                adjgrad[2] = -x2*l1*cx3 + mu*l2*cx3/(x1*x1) + l3*sx3*(x2/x1 - mu/(x1*x1*x2)) + x2*l4*sx3*cx6/(x1*cx5) + x2*l5*sx3*sx6/x1 - l6*(T*su2*sx3/(m*x2*cx3*cx3) + x2*sx3*sx5*cx6/(x1*cx5));
                adjgrad[3] = 0;
                adjgrad[4] = x2*cx3*cx6*(-l4*sx5 + l6)/(x1*cx5*cx5);
                adjgrad[5] = x2*cx3*(l4*sx6 - l5*cx5*cx6 - l6*sx5*sx6)/(x1*cx5);
            }

        //______________________________________________________________________
            template<class InType, class AdjVarType, class OutType>
            inline void max_adjointtransposegradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

                using Scalar = typename InType::Scalar;
                using OScalar = typename OutType::Scalar;
                using UVec = Matrix<OScalar, 2, 1>;
                using std::cos;
                using std::sin;
                using std::atan;
                
                Scalar x2, x3;
                x2 = x[1];
                x3 = x[2];
                
                typename AdjVarType::Scalar l2, l3, l6;
                l2 = adj[1];
                l3 = adj[2];
                l6 = adj[5];

                MatrixBase<OutType> & u = u_.const_cast_derived();

                std::vector<UVec> uOpts = ext_adjTransGrad_subroutine(x, adj);

                // Use 2nd adjoint derivative to find min or max solution
                u = UVec::Zero();
                for(int i=0; i<4; i++){
                    OScalar H11, H12, H22;
                    H11 = -T*cos(uOpts[i][1])*(x2*l2*cos(uOpts[i][0]) + l3*sin(uOpts[i][0]))/(m*x2);
                    H12 = T*sin(uOpts[i][1])*(x2*l2*sin(uOpts[i][0]) - l3*cos(uOpts[i][0]))/(m*x2);
                    H22 = -T*((l3*sin(uOpts[i][0])*cos(uOpts[i][1]) + l6*sin(uOpts[i][1])/cos(x3))/x2 + l2*cos(uOpts[i][0])*cos(uOpts[i][1]))/m;
                    if( (H11 + H22 <= 0) && (H11*H22 - H12*H12 <=0) ){  // Second-order maximum condition
                        u = uOpts[i];
                        return;
                    }
                }
            }

            template<class InType, class AdjVarType, class OutType>
            inline void min_adjointtransposegradient_control(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj, MatrixBase<OutType> const & u_) const {

                using Scalar = typename InType::Scalar;
                using OScalar = typename OutType::Scalar;
                using UVec = Matrix<OScalar, 2, 1>;
                using std::cos;
                using std::sin;
                using std::atan;
                
                Scalar x2, x3;
                x2 = x[1];
                x3 = x[2];

                typename AdjVarType::Scalar l2, l3, l6;
                l2 = adj[1];
                l3 = adj[2];
                l6 = adj[5];

                MatrixBase<OutType> & u = u_.const_cast_derived();

                std::vector<UVec> uOpts = ext_adjTransGrad_subroutine(x, adj);

                // Use 2nd adjoint derivative to find min or max solution
                u = UVec::Zero();
                for(int i=0; i<4; i++){
                    OScalar H11, H12, H22;
                    H11 = -T*cos(uOpts[i][1])*(x2*l2*cos(uOpts[i][0]) + l3*sin(uOpts[i][0]))/(m*x2);
                    H12 = T*sin(uOpts[i][1])*(x2*l2*sin(uOpts[i][0]) - l3*cos(uOpts[i][0]))/(m*x2);
                    H22 = -T*((l3*sin(uOpts[i][0])*cos(uOpts[i][1]) + l6*sin(uOpts[i][1])/cos(x3))/x2 + l2*cos(uOpts[i][0])*cos(uOpts[i][1]))/m;
                    if( (H11 + H22 >= 0) && (H11*H22 - H12*H12 >=0) ){  // Second-order minimum condition
                        u = uOpts[i];
                        return;
                    }
                }

            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Methods
            template<class InType, class AdjVarType, class OutType>
            inline std::vector<Matrix<typename OutType::Scalar, 2, 1> > ext_adjTransGrad_subroutine(const MatrixBase<InType> & x, const MatrixBase<AdjVarType> & adj) const {
                using Scalar = typename InType::Scalar;
                using OScalar = typename OutType::Scalar;
                using UVec = Matrix<OScalar, 2, 1>;
                using std::cos;
                using std::sin;
                using std::atan;

                Scalar x1, x2, x3, x4, x5, x6, u1, u2;
                x1 = x[0];
                x2 = x[1];
                x3 = x[2];
                x4 = x[3];
                x5 = x[4];
                x6 = x[5];
                u1 = x[7];
                u2 = x[8];

                typename AdjVarType::Scalar l1, l2, l3, l4, l5, l6;
                l1 = adj[0];
                l2 = adj[1];
                l3 = adj[2];
                l4 = adj[3];
                l5 = adj[4];
                l6 = adj[5];

                OScalar u1a, u1b, u2aa, u2ab, u2ba, u2bb;
                u1a = atan(l3/(x2*l2));
                u1b = u1a + M_PI;
                u2aa = atan((l6 - l3*cos(x3)*sin(u1a))/(x2*l2*cos(x3)*cos(u1a)));
                u2ab = u2aa + M_PI;
                u2ba = atan((l6 - l3*cos(x3)*sin(u1b))/(x2*l2*cos(x3)*cos(u1b)));
                u2bb = u2ba + M_PI;

                std::vector<UVec> uOpts;
                uOpts.resize(4);
                uOpts[0][0] = u1a;
                uOpts[0][1] = u2aa;
                uOpts[1][0] = u1a;
                uOpts[1][1] = u2ab;
                uOpts[2][0] = u1b;
                uOpts[2][1] = u2ba;
                uOpts[3][0] = u1b;
                uOpts[3][1] = u2bb;

                return uOpts;
            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            double mu;  // Gravitational Parameter (G*m)
            double T;  // Thrust
            double m;  // Spacecraft Mass


    };

}