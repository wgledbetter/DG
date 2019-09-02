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


        //======================================================================
            /// Setup
            inline void set_bounds(int d, double lo, double hi){

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

            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline void index_one2n(int idx, std::vector<int> &vec) const {
                int denom = std::accumulate(nDisc.begin(), nDisc.end(), 1, std::multiplies<int>());
                for(int i=dim-1; i>-1; i--){
                    denom /= nDisc[i];
                    vec[i] = idx % denom;
                    idx -= vec[i]*nDisc[i];
                }
            }

            inline void index_n2one(std::vector<int> vec, int &idx) const {
                int factor = std::accumulate(nDisc.begin(), nDisc.end(), 1, std::multiplies<int>());
                idx = 0;
                for(int i=dim-1, i>-1; i--){
                    factor /= nDisc[i];
                    idx += factor*vec[i];
                }
            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            std::vector<int> nDisc(dim);


    }

}