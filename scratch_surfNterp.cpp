// Make an outline for your 'surface interpolation with HJB constraints' idea

#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

template<int dim, class VectorFunction, class Mesh>
struct SurfNTerp : Mesh {

};

template<class VectorFunction>
struct SurfNTerp<2, VectorFunction, Mesh<2, MeshType::RecTriangular> > : Mesh<2, MeshType::RecTriangular> {

    template<class Scalar>
    using Vector = Eigen::Matrix<Scalar, 2, 1>;

    void set_seed(vector<int> verts, vector<double> vals){
                seedVert = verts;
                seedVal = vals;
            }

    void init() {
        weights.resize(Mesh::nVerts);
        backTri.resize(Mesh::nVerts);
        val.resize(Mesh::nVerts);

        Vector<double> vert;
        Vector<double> dVert;
        Vector<double> back;
        int t;
        std::vector<int> backprop_active;
        Eigen::Matrix<double, 3, 3> weightMat;
        Eigen::Matrix<double, 3, 1> weightVec;
        for(int i=0; i<Mesh::nVerts; i++){
            // Calc VectorFunction at all meshpoints
            vert = Mesh::verts[i];
            dVert = vFunc->compute(vert);

            // Backprop by dt
            back = vert - dt*dVert;

            // Calc insideTriangle for all backprops
            t = Mesh::insideTriangle(back);

            // If backprop is in bounds, calc interpolation weights for each backprop inside its respective triangle
            if(t >= 0){
                backprop_active.push_back(i);
                backTri[i] = t;

                weightMat.row(2) = Eigen::Matrix<double, 1, 3>::Ones();
                for(int j=0; j<3; j++){
                    weightMat.block<2,1>(0, j) = Mesh::verts[Mesh::triangles[t][j]];
                }
                weightVec[0] = back[0];
                weightVec[1] = back[1];
                weightVec[2] = 1;
                weights[i] = weightMat.ldlt().solve(weightVec);
            }
        }
        // Calc total number of constraints
        int nBckp = backprop_active.size();
        int nSeed = seedVert.size();
        int nCon = nBckp + nSeed;
        A.resize(nCon, Mesh::nVerts);
        b.resize(nCon);

        // Store backprops in correct location in A
        for(int i=0; i<nBckp; i++){
            int v = backprop_active[i];
            A(i,v) = -1;
            for(int j=0; j<3; j++){
                A(i, Mesh::triangles[backTri[v]][j]) = weights[v][j];
            }

            // Calc HJB constraints
            // MORE CODE GOES HERE
        }
        // Add seed constraints
        for(int i=0; i<nSeed; i++){
            int idx = nBckp + i;
            A(idx, seedVert[i]) = 1;

            // Add seeds to 'b'
            b[idx] = seedVal[i];
        }

    }

    void solve() {
        Eigen::PardisoLDLT< Eigen::SparseMatrix<double> > solver;
        solver.compute(A);
        if(solver.info() == Eigen::Success){
            val = solver.solve(b);
        }
        // val = A.ldlt().solve(b);
    }

    void initVert(int i) {
        // Calc VectorFunction
        Vector<double> vert = Mesh::verts[i];
        Vector<double> dVert = fVunc->compute(vert);

        // Backprop by dt
        Vector<double> back = vert - dt*dVert;

        // Which triangle is it in?
        int t = Mesh::insideTriangle(back);
        backTri[i] = t;

        // Solve interpolation weights
        Eigen::Matrix<double, 3, 3> weightMat;
        weightMat.row(2) = Eigen::Matrix<double, 1, 3>::Ones();
        for(int j=0; j<3; j++){
            weightMat.block<2,1>(0, j) = Mesh::verts[Mesh::triangles[t][j]];
        }
        Eigen::Matrix<double, 3, 1> weightVec;
        weightVec[0] = back[0];
        weightVec[1] = back[1];
        weightVec[2] = 1;
        weights[i] = weightMat.ldlt().solve(weightVec);

        // Store weights in A
        A(i,i) = -1;
        for(int j=0; j<3; j++){
            A(i, Mesh::triangles[t][j]) = weights[i][j];
        }
    }

    std::vector< Eigen::Matrix<double, 3, 1> > weights;
    std::vector<int> backTri;
    Eigen::VectorXd val;

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;

};