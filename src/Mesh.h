#pragma once

#include <vector>
#include <numeric>
#include <Eigen/Core>

namespace WGL_DG {

    enum MeshType {
        Regular,
        RecTriangular,
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
                Vector<int> idxVec = Vector<int>::Zero();
                for(int i=0; i<nVert; i++){
                    verts[i] = calcVert(idxVec);
                    neighbors[i] = calcNeighbors(idxVec);
                    incrementNdIndex(idxVec, nDisc);
                }
                
            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline void index_one2n(int idx, Vector<int> &vec, Vector<int> discretization) const {
                int denom = discretization.prod();
                for(int i=dim-1; i>-1; i--){
                    denom /= discretization[i];
                    vec[i] = idx % denom;
                    idx -= vec[i]*discretization[i];
                }
            }

            inline void index_n2one(Vector<int> vec, int &idx, Vector<int> discretization) const {
                int factor = discretization.prod();
                idx = 0;
                for(int i=dim-1; i>-1; i--){
                    if(vec[i] >= discretization[i]){
                        idx = std::numeric_limits<int>::infinity();
                        return;
                    }
                    factor /= discretization[i];
                    idx += factor*vec[i];
                }
            }
            
            inline void index_n2one(Vector<int> vec, int &idx, int disc) const {
                int factor = pow(disc, dim);
                idx = 0;
                for(int i=dim-1; i>-1; i--){
                    if(vec[i] >= disc){
                        idx = std::numeric_limits<int>::infinity();
                        return;
                    }
                    factor /= disc;
                    idx += factor*vec[i];
                }
            }

            inline void incrementNdIndex(Vector<int> &vec, Vector<int> discretization) const {
                for(int i=0; i<dim; i++){
                    if(vec[i]+1 < discretization[i]){
                        vec[i]++;
                        return;
                    }else{
                        vec[i] = 0;
                    }
                }
            }
            
            inline void incrementNdIndex(Vector<int> &vec, int disc) const {
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
                return ( loBounds + dx.cwiseProduct(idxVec.template cast<double>()) );
            }

        //______________________________________________________________________
            inline std::vector<int> calcNeighbors(Vector<int> idxVec) const {

                std::vector<int> out = {};

                const Vector<int> lowVec = -Vector<int>::Ones();
                const Vector<int> dxVec = Vector<int>::Ones();

                Vector<int> nhbIdxVec = Vector<int>::Zero();
                Vector<int> delta;
                Vector<int> nhbVec;
                int nb;

                int nNhb = pow(3, dim);
                for(int i=0; i<nNhb; i++){
                    delta = lowVec + dxVec.cwiseProduct(nhbIdxVec);
                    nhbVec = idxVec + delta;
                    // If in bounds AND not central
                    bool cond1 = (nhbVec == nhbVec.cwiseAbs());  // No element is subzero
                    Vector<int> vec2 = (nDisc-nhbVec).array() - 1;
                    bool cond2 = ( vec2 == vec2.cwiseAbs() );  // No element is above nDisc
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

//==============================================================================
    template<>
    struct Mesh<2, MeshType::RecTriangular> {

        public:
            /// Typedefs
            template<class Scalar>
            using Vector = Eigen::Matrix<Scalar, 2, 1>;


        //======================================================================
            /// Constructors
            Mesh(){
                nDisc.resize(2);
                dx.resize(2);
                loBounds.resize(2);
                hiBounds.resize(2);
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
                for(int i=0; i<2; i++){
                    dx[i] = (hiBounds[i] - loBounds[i])/(nDisc[i]-1);
                }

                nVert = nDisc.prod();
                verts.resize(nVert);
                neighbors.resize(nVert);

                Vector<int> idxVec = Vector<int>::Zero();
                for(int i=0; i<nVert; i++){
                    verts[i] = calcVert(idxVec);
                    neighbors[i] = calcNeighbors(idxVec);
                    incrementNdIndex(idxVec, nDisc);
                }

                Vector<int> triDisc;
                triDisc[0] = nDisc[0]-1;
                triDisc[1] = 2*(nDisc[1]-1);
                nTri = triDisc.prod();
                triangles.resize(nTri);
                Vector<int> triIdx = Vector<int>::Zero();
                for(int i=0; i<nTri; i++){
                    triangles[i] = calcTri(triIdx);
                    incrementNdIndex(triIdx, triDisc);
                }

            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline Vector<double> calcVert(Vector<int> idxVec) const {
                return ( loBounds + dx.cwiseProduct(idxVec.template cast<double>()) );
            }

        //______________________________________________________________________
            inline std::vector<int> calcNeighbors(Vector<int> idxVec) const {

                std::vector<int> out = {};

                // Make "positive slope" triangles -> /, not \.
                std::vector< Vector<int> > triNhbs;
                triNhbs.resize(6);
                triNhbs[0] << -1, -1;
                triNhbs[1] << 0, -1;
                triNhbs[2] << 1, 0;
                triNhbs[3] << 1, 1;
                triNhbs[4] << 0, 1;
                triNhbs[5] << -1, 0;

                Vector<int> nhbVec;
                int nb;

                for(int i=0; i<6; i++){
                    nhbVec = idxVec + triNhbs[i];

                    bool cond1 = (nhbVec == nhbVec.cwiseAbs());
                    Vector<int> vec2 = (nDisc-nhbVec).array() -1;
                    bool cond2 = ( vec2 == vec2.cwiseAbs() );
                    if(cond1 && cond2){
                        index_n2one(nhbVec, nb, nDisc);
                        out.push_back(nb);
                    }
                }

                return out;

            }

        //______________________________________________________________________
            inline Eigen::Matrix<int, 3, 1> calcTri(Vector<int> idxVec){

                Eigen::Matrix<int, 3, 1> tri;
                bool odd = idxVec[0] % 2;
                int x, y;
                Vector<int> vec;
                int idx;

                if(odd){
                    x = idxVec[0];
                    y = idxVec[1]/2 + 1;

                    vec[0] = x;
                    vec[1] = y;
                    index_n2one(vec, idx, nDisc);
                    tri[0] = idx;

                    vec[0] = x;
                    vec[1] = y-1;
                    index_n2one(vec, idx, nDisc);
                    tri[1] = idx;

                    vec[0] = x+1;
                    vec[1] = y;
                    index_n2one(vec, idx, nDisc);
                    tri[2] = idx;
                }else{
                    x = idxVec[0];
                    y = idxVec[1]/2;

                    vec[0] = x;
                    vec[1] = y;
                    index_n2one(vec, idx, nDisc);
                    tri[0] = idx;

                    vec[0] = x+1;
                    vec[1] = y;
                    index_n2one(vec, idx, nDisc);
                    tri[1] = idx;

                    vec[0] = x+1;
                    vec[1] = y+1;
                    index_n2one(vec, idx, nDisc);
                    tri[2] = idx;
                }

                return tri;

            }

        //______________________________________________________________________
            inline void incrementNdIndex(Vector<int> &vec, Vector<int> discretization) const {
                for(int i=0; i<2; i++){
                    if(vec[i]+1 < discretization[i]){
                        vec[i]++;
                        return;
                    }else{
                        vec[i] = 0;
                    }
                }
            }

            inline void incrementNdIndex(Vector<int> &vec, int disc) const {
                for(int i=0; i<2; i++){
                    if(vec[i]+1 < disc){
                        vec[i]++;
                        return;
                    }else{
                        vec[i] = 0;
                    }
                }
            }

            inline void index_n2one(Vector<int> vec, int &idx, Vector<int> discretization) const {
                int factor = discretization.prod();
                idx = 0;
                if(vec[1] >= discretization[1]){
                    idx = std::numeric_limits<int>::infinity();
                    return;
                }
                factor /= discretization[1];
                idx += factor*vec[1];
                if(vec[0] >= discretization[0]){
                    idx = std::numeric_limits<int>::infinity();
                    return;
                }
                factor /= discretization[0];
                idx += factor*vec[0];
            }

            inline void index_n2one(Vector<int> vec, int &idx, int disc) const {
                int factor = disc*disc;
                idx = 0;
                if(vec[1] >= disc){
                    idx = std::numeric_limits<int>::infinity();
                    return;
                }
                factor /= disc;
                idx += factor*vec[1];
                if(vec[0] >= disc){
                    idx = std::numeric_limits<int>::infinity();
                    return;
                }
                factor /= disc;
                idx += factor*vec[0];
            }


        //======================================================================
            /// Properties
            std::vector< Vector<double> > verts;
            std::vector< std::vector<int> > neighbors;
            int nVert;
            int nTri;
            std::vector< Eigen::Matrix<int, 3, 1> > triangles;


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            Vector<int> nDisc;
            Vector<double> dx;
            Vector<double> loBounds;
            Vector<double> hiBounds;

    };

}