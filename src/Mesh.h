#pragma once

#include <vector>
#include <numeric>
#include <Eigen/Core>

namespace WGL_DG {

    enum MeshType {
        Regular,
    }

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
                for(int i=0; i<nD.size(), i++){
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
                nVert = std::accumulate(nDisc.begin(), nDisc.end(), 1, std::multiplies<int>());
                verts.resize(nVert);
                neighbors.resize(nVert);

                // Calculate vertex locations
                Eigen::Matrix<int, dim, 1> idxVec = Eigen::Matrix<int, dim, 1>::Zeros();
                for(int i=0; i<nVert; i++){
                    verts[i] = calcVert(idxVec);
                    incrementNdIndex(idxVec);
                }
            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline void index_one2n(int idx, Eigen::Matrix<int, dim, 1> &vec) const {
                int denom = std::accumulate(nDisc.begin(), nDisc.end(), 1, std::multiplies<int>());
                for(int i=dim-1; i>-1; i--){
                    denom /= nDisc[i];
                    vec[i] = idx % denom;
                    idx -= vec[i]*nDisc[i];
                }
            }

            inline void index_n2one(Eigen::Matrix<int, dim, 1> vec, int &idx) const {
                int factor = std::accumulate(nDisc.begin(), nDisc.end(), 1, std::multiplies<int>());
                idx = 0;
                for(int i=dim-1, i>-1; i--){
                    factor /= nDisc[i];
                    idx += factor*vec[i];
                }
            }

            inline void incrementNdIndex(Eigen::Matrix<int, dim, 1> &vec) const {
                for(int i=0; i<dim; i++){
                    if(vec[i]+1 < nDisc[i]){
                        vec[i]++;
                    }else{
                        vec[i] = 0;
                    }
                }
            }

        //______________________________________________________________________
            inline Eigen::Matrix<double, dim, 1> calcVert(Eigen::Matrix<int, dim, 1> idxVec) const {
                return loBounds + dx.cwiseProduct(idxVec);
            }

        //______________________________________________________________________
            inline std::vector<int> calcNeighbors(Eigen::Matrix<int, dim, 1> idxVec) const {
                // All permutations of dim instances of {-1, 0, 1} except all zeros
                int nNhb = pow(3, dim) - 1;  // Max possible neighbors

            }

        //======================================================================
            /// Properties
            std::vector< Eigen::Matrix<double, dim, 1> > verts;
            std::vector< std::vector<int> > neighbors;
            int nVert;


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            std::vector<int> nDisc;
            Eigen::Matrix<double, dim, 1> dx;
            Eigen::Matrix<double, dim, 1> loBounds;
            Eigen::Matrix<double, dim, 1> hiBounds;

    }

}