#pragma once

#include <vector>
#include <limits>
#include <algorithm>
#include <Eigen/Core>

#include "ctpl.h"
#include "Mesh.h"
#include "TensorFunction.h"

using namespace std;

namespace WGL_DG {

    template<int dim, class TensorFunction, class Mesh>
    struct EikonalSolution : Mesh {

        public:
            /// Typedefs
            //using MeshBase = Mesh;


        //======================================================================
            /// Constructors
            EikonalSolution(){
                vector<double> MeshData::val;
                vector< Eigen::Matrix<double, dim, dim> > MeshData::spd;
                vector<bool> MeshData::is_seed;
                vector<bool> MeshData::is_active;
                vector< vector<double> > MeshData::costToGo;

                vector<int> FIM::activeList;
                vector<int> FIM::removeThese;
                vector<int> FIM::addThese;
                double FIM::convTol;
            }


        //======================================================================
            /// Setup
            inline void set_spd_func(TensorFunction *tf){
                this->tFunc = tf;
            }

        //______________________________________________________________________
            inline void set_seed(vector<int> verts, vector<double> vals){
                seedVert = verts;
                seedVal = vals;
            }


        //======================================================================
            /// Primary Methods
            inline void init(){
                // Init mesh
                MeshData::val.resize(Mesh::nVert);
                MeshData::spd.resize(Mesh::nVert);
                MeshData::is_seed.resize(Mesh::nVert);
                MeshData::is_active.resize(Mesh::nVert);
                fill(MeshData::val.begin(), MeshData::val.end(), numeric_limits<double>::infinity());
                fill(MeshData::is_seed.begin(), MeshData::is_seed.end(), false);
                fill(MeshData::is_active.begin(), MeshData::is_active.end(), false);

                // Init seed
                fill(MeshData::val.begin(), MeshData::val.end(), numeric_limits<double>::infinity());
                fill(MeshData::is_seed.begin(), MeshData::is_seed.end(), false);
                for(int i=0; i<seedVert.size(); i++){
                    MeshData::val[seedVert[i]] = seedVal[i];
                    MeshData::is_seed[seedVert[i]] = true;
                }
                for(int i=0; i<seedVert.size(); i++){
                    for(int j=0; j<Mesh::neighbors[seedVert[i]]; j++){
                        int nb = Mesh::neighbors[seedVert[i]][j];
                        if(!MeshData::is_seed[nb] && !MeshData::is_active[nb]){
                            MeshData::is_active[nb] = true;
                            FIM::activeList.push_back(nb);
                        }
                    }
                }

                // Calc (1/sqrt((x-y).trps()*M*(x-y))) for all connections
                MeshData::costToGo.resize(Mesh::nVert);
                for(int i=0; i<Mesh::nVert; i++){
                    int nNhb = Mesh::neighbors[i].size();
                    MeshData::costToGo[i].resize(nNhb);
                    Eigen::Matrix<double, dim, 1> x = Mesh::verts[i];
                    Eigen::Matrix<double, dim, dim> M = tFunc->compute(x);
                    for(int j=0; j<nNhb; j++){
                        int nb = Mesh::neighbors[i][j];
                        Eigen::Matrix<double, dim, 1> y = Mesh::verts[nb];
                        Eigen::Matrix<double, dim, 1> invDist = (x-y).cwiseInverse();
                        double radical = invDist.transpose()*M*invDist;
                        costToGo[i][j] = 1/sqrt(rad);
                    }
                }

                // Set up memory structures for high speed & conflict prevention
                setupMemory();
            }

        //______________________________________________________________________
            inline void compute(){
                while(FIM::activeList.size()){
                    for(int i=0; i<FIM::activeList.size(); i++){
                        activeLoop(v);
                    }
                    updateActiveList();
                }
            }


        //======================================================================
            /// Properties


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline double solvePDE(int v) const {

            }

        //______________________________________________________________________
            inline void activeLoop(int v){
                double p = MeshData::val[v];
                double q = solvePDE(v);
                if(p > q){
                    MeshData::val[v] = q;
                }
                if(abs(p-q) < FIM::convTol){
                    for(int i=0; i<Mesh::neighbors[v].size(); i++){
                        int nb = Mesh::neighbors[v][i];
                        if(!FIM::is_active[nb]){
                            double ip = MeshData::val[nb];
                            double iq = solvePDE(nb);
                            if(ip > iq){
                                MeshData::val[nb] = iq;
                                FIM::addThese.push_back(nb);
                            }
                        }
                    }
                    FIM::removeThese.push_back(v);
                }
            }

        //______________________________________________________________________
            inline void updateActiveList(){
                vector<int>::iterator findInVec;
                for(int i=0; i<FIM::removeThese.size(); i++){
                    int rm = FIM::removeThese[i];
                    findInVec = find(FIM::activeList.begin(), FIM::activeList.end(), rm);
                    if(findInVec != FIM::activeList.end()){  // If it was found, remove it
                        FIM::activeList.erase(findInVec);
                        MeshData::is_active[rm] = false;
                    }
                }
                FIM::removeThese.clear();
                for(int i=0; i<FIM::addThese.size(); i++){
                    int ad = FIM::addThese[i];
                    findInVec = find(FIM::activeList.begin(), FIM::activeList.end(), ad);
                    if(findInVec == FIM::activeList.end()){  // If it was not found, add it
                        FIM::activeList.push_back(ad);
                        MeshData::is_active[ad] = true;
                    }
                }
                FIM::addThese.clear();
            }

        //______________________________________________________________________
            inline void setupMemory(){

            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            TensorFunction* tFunc;

            vector<int> seedVert;
            vector<double> seedVal;

            struct MeshData {
                static vector<double> val;
                static vector< Eigen::Matrix<double, dim, dim> > spd;
                static vector<bool> is_seed;
                static vector<bool> is_active;
                static vector< vector<double> > costToGo;  // costToGo[i][j] is cost to go from vertex i to its jth neighbor
            };

            struct FIM {
                static vector<int> activeList;
                static vector<int> removeThese;
                static vector<int> addThese;
                static double convTol;
            }


    };

}