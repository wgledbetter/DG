#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>

#include <Eigen/Sparse>
#include <Eigen/SPQRSupport>

#include "Mesh.h"

using namespace std;
using namespace Eigen;

namespace WGL_DG {
    
    /// Template
    template<int dim, class VectorFunction, class Mesh>
    struct SurfNTerp : Mesh {
        
    };
    
////////////////////////////////////////////////////////////////////////////////
    /// Specializations
    template<class VectorFunction>
    struct SurfNTerp<2, VectorFunction, Mesh<2, MeshType::RecTriangular> > : Mesh<2, MeshType::RecTriangular> {
        
        public:
            /// Typedefs
            template<class Scalar>
            using Vector = Matrix<Scalar, 2, 1>;
            
            
        //======================================================================
            /// Constructors
            SurfNTerp(){
                
            }
            
            
        //======================================================================
            /// Setup
            inline void set_seed(vector<int> verts, vector<double> vals){
                seedVert = verts;
                seedVal = vals;
            }
            
        //______________________________________________________________________
            inline void set_dt(double in){
                dt = in;
            }
            
        //______________________________________________________________________
            inline void set_vec_func(VectorFunction *vf){
                vFunc = vf;
            }
            
            
        //======================================================================
            /// Primary Methods
            inline void init(){
                weights.resize(Mesh::nVert);
                backTri.resize(Mesh::nVert);
                val.resize(Mesh::nVert);
                
                Vector<double> vert;
                Vector<double> dVert;
                Vector<double> back;
                int t;
                std::vector<int> backprop_active;
                Eigen::Matrix<double, 3, 3> weightMat;
                Eigen::Matrix<double, 3, 1> weightVec;
                for(int i=0; i<Mesh::nVert; i++){
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
                A.resize(nCon, Mesh::nVert);
                b.resize(nCon);
                typedef Triplet<double> Trip;
                vector<Trip> tripVec;

                // Store backprops in correct location in A
                for(int i=0; i<nBckp; i++){
                    int v = backprop_active[i];
                    
                    bool source_constraint = false;
                    for(int j=0; j<3; j++){
                        // A(i, Mesh::triangles[backTri[v]][j]) = weights[v][j];
                        if(Mesh::triangles[backTri[v]][j] == v){
                            tripVec.push_back(Trip(i, v, weights[v][j]-1));
                            source_constraint = true;
                        }else{
                            tripVec.push_back(Trip(i, Mesh::triangles[backTri[v]][j], weights[v][j]));
                        }
                    }
                    if(!source_constraint){
                        // A(i,v) = -1;
                        tripVec.push_back(Trip(i, v, -1));
                    }

                    // Calc HJB constraints
                    // Need integral payoff function
                    // Just use dt for now - will give travel time value function
                    b[i] = dt;
                }
                // Add seed constraints
                for(int i=0; i<nSeed; i++){
                    int idx = nBckp + i;
                    // A(idx, seedVert[i]) = 1;
                    tripVec.push_back(Trip(idx, seedVert[i], 1));

                    // Add seeds to 'b'
                    b[idx] = seedVal[i];
                }
                
                // Fill A from Triplets
                A.setFromTriplets(tripVec.begin(), tripVec.end());
                A.makeCompressed();
            }
            
        //______________________________________________________________________
            inline void solve(){
                SPQR< SparseMatrix<double> > solver;
                solver.compute(A);
                if(solver.info() == Eigen::Success){
                    val = solver.solve(b);
                    if(solver.info() == Eigen::Success){
                        cout << "  Solved" << endl;
                    }
                }
            }
            
            
        //======================================================================
            inline void textFileOutput(){
                std::time_t t = time(0);
                struct tm * now = localtime( & t );
                char currentTime_buffer [80];
                std::strftime (currentTime_buffer, 80, "%F_%T", now);
                
                Eigen::IOFormat Hector(Eigen::FullPrecision, 0, "", ", ", "", "", "", "");
                
                std::ostringstream foldName;
                foldName << "out/" << currentTime_buffer;
                if( mkdir(foldName.str().c_str(), 0777) != -1 ){
                    // Save Value and Vertex
                    std::ostringstream valFname;
                    valFname << foldName.str().c_str() << "/SurfNTerp_Value.txt";
                    std::ofstream save_val;
                    save_val.open(valFname.str());
                    
                    std::ostringstream vertFname;
                    vertFname << foldName.str().c_str() << "/SurfNTerp_Vertex.txt";
                    std::ofstream save_vert;
                    save_vert.open(vertFname.str());
                    
                    if( save_val.is_open() && save_vert.is_open()){
                        for(int i=0; i<Mesh::nVert; i++){
                            save_val << val[i] << '\n';
                            save_vert << Mesh::verts[i].format(Hector) << '\n';
                        }
                    }
                    save_val.close();
                    save_vert.close();
                
                }
            }
            
            
    ////////////////////////////////////////////////////////////////////////////
        protected:
            inline void something(){
                
            }
            
            
    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            VectorFunction* vFunc;

            vector< Matrix<double, 3, 1> > weights;
            vector<int> backTri;
            VectorXd val;

            SparseMatrix<double> A;
            VectorXd b;

            vector<int> seedVert;
            vector<double> seedVal;

            double dt;
        
    };
    
}