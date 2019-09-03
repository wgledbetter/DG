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
                meshData.val.resize(Mesh::nVert);
                meshData.spd.resize(Mesh::nVert);
                meshData.is_seed.resize(Mesh::nVert);
                meshData.is_active.resize(Mesh::nVert);
                fill(meshData.val.begin(), meshData.val.end(), numeric_limits<double>::infinity());
                fill(meshData.is_seed.begin(), meshData.is_seed.end(), false);
                fill(meshData.is_active.begin(), meshData.is_active.end(), false);

                // Init seed
                fill(meshData.val.begin(), meshData.val.end(), numeric_limits<double>::infinity());
                fill(meshData.is_seed.begin(), meshData.is_seed.end(), false);
                for(int i=0; i<seedVert.size(); i++){
                    meshData.val[seedVert[i]] = seedVal[i];
                    meshData.is_seed[seedVert[i]] = true;
                }
                for(int i=0; i<seedVert.size(); i++){
                    for(int j=0; j<Mesh::neighbors[seedVert[i]]; j++){
                        int nb = Mesh::neighbors[seedVert[i]][j];
                        if(!meshData.is_seed[nb] && !meshData.is_active[nb]){
                            meshData.is_active[nb] = true;
                            fim.activeList.push_back(nb);
                        }
                    }
                }

                // Calc (1/sqrt((x-y).trps()*M*(x-y))) for all connections
                meshData.costToGo.resize(Mesh::nVert);
                for(int i=0; i<Mesh::nVert; i++){
                    int nNhb = Mesh::neighbors[i].size();
                    meshData.costToGo[i].resize(nNhb);
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
                while(fim.activeList.size()){
                    for(int i=0; i<fim.activeList.size(); i++){
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
                vector<double> sols;
                for(int i=0; i<Mesh::neighbors[v].size(); i++){
                    int nb = Mesh::neighbors[v][i];
                    if(meshData.val[nb] == numeric_limits<double>::infinity()){
                        // Skip
                    }else{
                        sols.push_back(costToGo[v][i] + meshData.val[nb]);
                    }
                }
                return *min_element(sols.begin(), sols.end());
            }

        //______________________________________________________________________
            inline void activeLoop(int v){
                double p = meshData.val[v];
                double q = solvePDE(v);
                if(p > q){
                    meshData.val[v] = q;
                }
                if(abs(p-q) < fim.convTol){
                    for(int i=0; i<Mesh::neighbors[v].size(); i++){
                        int nb = Mesh::neighbors[v][i];
                        if(!fim.is_active[nb]){
                            double ip = meshData.val[nb];
                            double iq = solvePDE(nb);
                            if(ip > iq){
                                meshData.val[nb] = iq;
                                fim.addThese.push_back(nb);
                            }
                        }
                    }
                    fim.removeThese.push_back(v);
                }
            }

        //______________________________________________________________________
            inline void updateActiveList(){
                vector<int>::iterator findInVec;
                for(int i=0; i<fim.removeThese.size(); i++){
                    int rm = fim.removeThese[i];
                    findInVec = find(fim.activeList.begin(), fim.activeList.end(), rm);
                    if(findInVec != fim.activeList.end()){  // If it was found, remove it
                        fim.activeList.erase(findInVec);
                        meshData.is_active[rm] = false;
                    }
                }
                fim.removeThese.clear();
                for(int i=0; i<fim.addThese.size(); i++){
                    int ad = fim.addThese[i];
                    findInVec = find(fim.activeList.begin(), fim.activeList.end(), ad);
                    if(findInVec == fim.activeList.end()){  // If it was not found, add it
                        fim.activeList.push_back(ad);
                        meshData.is_active[ad] = true;
                    }
                }
                fim.addThese.clear();
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
                vector<double> val;
                vector< Eigen::Matrix<double, dim, dim> > spd;
                vector<bool> is_seed;
                vector<bool> is_active;
                vector< vector<double> > costToGo;  // costToGo[i][j] is cost to go from vertex i to its jth neighbor
            };
            MeshData meshData;

            struct FIM {
                vector<int> activeList;
                vector<int> removeThese;
                vector<int> addThese;
                double convTol;
            }
            FIM fim;


    };

}