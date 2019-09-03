#pragma once

#include <vector>
#include <numeric>
#include <Eigen/Core>

namespace WGL_DG {

    enum MeshType {
        Regular,
    };

    /// Template
    template<int dim, int meshType>
    struct Mesh {

    };

////////////////////////////////////////////////////////////////////////////////
    /// Specializations
    template<int dim>
    struct Mesh<dim, MeshType::Regular> {

        public:
            /// Typedefs
            template<class Scalar>
            using Vector = Eigen::Matrix<Scalar, dim, 1>;


        //======================================================================
            /// Constructors
            Mesh(){
                nDisc.resize(dim);
                dx.resize(dim);
                loBounds.resize(dim);
                hiBounds.resize(dim);
            }


        //======================================================================
            /// Setup
            inline void set_bounds(int d, double lo, double hi){
                loBounds[d] = lo;
                hiBounds[d] = hi;
            }
            inline void set_bounds(std::vector< std::vector<double> > bds){
                for(int i=0; i<bds.size(); i++){
                    set_bounds(i, bds[i][0], bds[i][1]);
                }
            }
            inline void set_bounds(std::vector<double> lows, std::vector<double> highs){
                for(int i=0; i<lows.size(); i++){
                    set_bounds(i, lows[i], highs[i]);
                }
            }

        //______________________________________________________________________
            inline void set_nDisc(int d, int n){
                nDisc[d] = n;
            }
            inline void set_nDisc(std::vector<int> nD){
                for(int i=0; i<nD.size(); i++){
                    set_nDisc(i, nD[i]);
                }
            }


        //======================================================================
            /// Primary Methods
            inline void gen_mesh(){
                // Calc dt for each dimension
                for(int i=0; i<dim; i++){
                    dx[i] = (hiBounds[i] - loBounds[i])/(nDisc[i]-1);
                }

                // Calc total verts
                nVert = nDisc.prod();
                verts.resize(nVert);
                neighbors.resize(nVert);

                // Calculate vertex locations
                Vector<int> idxVec = Vector<int>::Zeros();
                for(int i=0; i<nVert; i++){
                    verts[i] = calcVert(idxVec);
                    neighbors[i] = calcNeighbors(idxVec);
                    incrementNdIndex(idxVec, nDisc);
                }
                
            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline void index_one2n(int idx, Vector<int> &vec, const Vector<int> &discretization) const {
                int denom = discretization.prod();
                for(int i=dim-1; i>-1; i--){
                    denom /= discretization[i];
                    vec[i] = idx % denom;
                    idx -= vec[i]*discretization[i];
                }
            }

            inline void index_n2one(Vector<int> vec, int &idx, const Vector<int> &discretization) const {
                int factor = discretization.prod();
                idx = 0;
                for(int i=dim-1; i>-1; i--){
                    factor /= discretization[i];
                    idx += factor*vec[i];
                }
            }
            
            inline void index_n2one(Vector<int> vec, int &idx, const int disc){
                int factor = pow(disc, dim);
                idx = 0;
                for(int i=dim-1; i>-1; i--){
                    factor /= disc;
                    idx += factor*vec[i];
                }
            }

            inline void incrementNdIndex(Vector<int> &vec, const Vector<int> &discretization) const {
                for(int i=0; i<dim; i++){
                    if(vec[i]+1 < discretization[i]){
                        vec[i]++;
                        return;
                    }else{
                        vec[i] = 0;
                    }
                }
            }
            
            inline void incrementNdIndex(Vector<int> &vec, const int disc){
                for(int i=0; i<dim; i++){
                    if(vec[i]+1 < disc){
                        vec[i]++;
                        return;
                    }else{
                        vec[i] = 0;
                    }
                }
            }
            

        //______________________________________________________________________
            inline Vector<double> calcVert(Vector<int> idxVec) const {
                return loBounds + dx.cwiseProduct(idxVec);
            }

        //______________________________________________________________________
            inline std::vector<int> calcNeighbors(const Vector<int> idxVec) const {

                std::vector<int> out;

                const Vector<int> lowVec = -Vector<int>::Ones();
                const Vector<int> dxVec = Vector<int>::Ones();

                Vector<int> nhbIdxVec = Vector<int>::Zeros();
                Vector<int> delta;
                Vector<int> nhbVec;
                int nb;

                int nNhb = pow(3, dim);
                for(int i=0; i<nNhb; i++){
                    delta = lowVec + dxVec.cwiseProduct(nhbIdxVec);
                    nhbVec = idxVec + delta;
                    // If in bounds AND not central
                    bool cond1 = (nhbVec == nhbVec.cwiseAbs());  // No element is subzero
                    bool cond2 = ( (nDisc-nhbVec) == (nDisc-nhbVec).cwiseAbs() );  // No element is above nDisc
                    bool cond3 = (delta != Vector<int>::Zero());  // This is not the central point
                    if(cond1 && cond2 && cond3){
                        index_n2one(nhbVec, nb, nDisc);
                        out.push_back(nb);
                    }
                    incrementNdIndex(nhbIdxVec, 3);
                }
                
                return out;
            }

        //======================================================================
            /// Properties
            std::vector< Vector<double> > verts;
            std::vector< std::vector<int> > neighbors;
            int nVert;


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            Vector<int> nDisc;
            Vector<double> dx;
            Vector<double> loBounds;
            Vector<double> hiBounds;

    };

}