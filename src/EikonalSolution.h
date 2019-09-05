#pragma once

#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <thread>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ctpl.h"
#include "Mesh.h"

using namespace std;

namespace WGL_DG {

    template<int dim, class TensorFunction, class Mesh>
    struct EikonalSolution : Mesh {

        public:
            /// Typedefs
            //using MeshBase = Mesh;
            template<class Scalar>
            using Vector = Eigen::Matrix<Scalar, dim, 1>;
            
            template<class Scalar>
            using Matrix = Eigen::Matrix<Scalar, dim, dim>;


        //======================================================================
            /// Constructors
            EikonalSolution(){
                fim.convTol = 1e-6;
                threads = thread::hardware_concurrency();

                fim.addTheseThread.resize(threads);
                fim.removeTheseThread.resize(threads);
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
                    for(int j=0; j < Mesh::neighbors[seedVert[i]].size(); j++){
                        int nb = Mesh::neighbors[seedVert[i]][j];
                        if(!meshData.is_seed[nb] && !meshData.is_active[nb]){
                            meshData.is_active[nb] = true;
                            fim.activeList.push_back(nb);
                        }
                    }
                }

                // Calc (1/sqrt((x-y).trps()*M*(x-y))) for all connections
                meshData.costToGo.resize(Mesh::nVert);
                vector< Matrix<double> > M;
                M.resize(Mesh::nVert);
                ctpl::thread_pool initThreads(threads);
                Vector<double> x;
                for(int i=0; i<Mesh::nVert; i++){
                    x = Mesh::verts[i];
                    initThreads.push( [this, x, M, i, &tFunc](int thr){ M[i] = tFunc->compute(x); } );
                }
                initThreads.stop(true);
                initThreads.restart(threads);
                
                for(int i=0; i<Mesh::nVert; i++){
                    int nNhb = Mesh::neighbors[i].size();
                    meshData.costToGo[i].resize(nNhb);
                    Vector<double> x = Mesh::verts[i];
                    Matrix<double> M = tFunc->compute(x);
                    for(int j=0; j<nNhb; j++){
                        int nb = Mesh::neighbors[i][j];
                        Vector<double> y = Mesh::verts[nb];
                        Vector<double> delta = x-y;
                        Vector<double> temp = M.ldlt().solve(delta);
                        double radical = delta.dot(temp);
                        meshData.costToGo[i][j] = sqrt(radical);
                    }
                }

                // Set up memory structures for high speed & conflict prevention
                setupMemory();
            }

        //______________________________________________________________________
            inline void compute(){
                ctpl::thread_pool fimThreads(threads);
                int vert;
                
                while(fim.activeList.size()){
                    fimThreads.restart(threads);
                    for(int i=0; i<fim.activeList.size(); i++){
                        vert = fim.activeList[i];
                        fimThreads.push( [this, vert](int thr){ this->activeLoop(vert, thr); } );
                    }
                    fimThreads.stop(true);
                    updateActiveList();
                }
            }


        //======================================================================
            /// Output
            inline void textFileOutput(){
                std::time_t t = time(0);
                struct tm * now = localtime( & t );
                char currentTime_buffer [80];
                std::strftime (currentTime_buffer, 80, "%F_%T", now);
                
                Eigen::IOFormat Hector(Eigen::FullPrecision, 0, "", ", ", "", "", "", "");
                
                std::ostringstream foldName;
                foldName << "out/" << currentTime_buffer;
                if( mkdir(foldName.str().c_str(), 0777) != -1 ){
                    {  // Save Value and Vertex
                        std::ostringstream valFname;
                        valFname << foldName.str().c_str() << "/Value.txt";
                        std::ofstream save_val;
                        save_val.open(valFname.str());
                        
                        std::ostringstream vertFname;
                        vertFname << foldName.str().c_str() << "/Vertex.txt";
                        std::ofstream save_vert;
                        save_vert.open(vertFname.str());
                        
                        if( save_val.is_open() && save_vert.is_open()){
                            for(int i=0; i<Mesh::nVert; i++){
                                save_val << meshData.val[i] << '\n';
                                save_vert << Mesh::verts[i].format(Hector) << '\n';
                            }
                        }
                        save_val.close();
                        save_vert.close();
                    }
                }
            }


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods
            inline double solvePDE(int v) const {
                vector<double> sols;
                int nNhb = Mesh::neighbors[v].size();
                for(int i=0; i<nNhb; i++){
                    int nb = Mesh::neighbors[v][i];
                    if(meshData.val[nb] == numeric_limits<double>::infinity()){
                        // Skip
                    }else{
                        sols.push_back(meshData.costToGo[v][i] + meshData.val[nb]);
                    }
                }
                return *min_element(sols.begin(), sols.end());
            }

        //______________________________________________________________________
            inline void activeLoop(const int v, const int thr){
                double p = meshData.val[v];
                double q = solvePDE(v);
                if(p > q){
                    meshData.val[v] = q;
                }
                if(abs(p-q) < fim.convTol){
                    for(int i=0; i<Mesh::neighbors[v].size(); i++){
                        int nb = Mesh::neighbors[v][i];
                        if(!meshData.is_active[nb]){
                            double ip = meshData.val[nb];
                            double iq = solvePDE(nb);
                            if(ip > iq){
                                meshData.val[nb] = iq;
                                fim.addTheseThread[thr].push_back(nb);
                            }
                        }
                    }
                    fim.removeTheseThread[thr].push_back(v);
                }
            }

        //______________________________________________________________________
            inline void updateActiveList(){
                // Amalgamate *TheseThreads into *These
                fim.removeThese.clear();
                fim.addThese.clear();
                for(int i=0; i<threads; i++){
                    for(int j=0; j<fim.addTheseThread[i].size(); j++){
                        fim.addThese.push_back(fim.addTheseThread[i][j]);
                    }
                    for(int j=0; j<fim.removeTheseThread[i].size(); j++){
                        fim.removeThese.push_back(fim.removeTheseThread[i][j]);
                    }
                    fim.addTheseThread[i].clear();
                    fim.removeTheseThread[i].clear();
                }

                vector<int>::iterator findInVec;
                sort(fim.removeThese.begin(), fim.removeThese.end());
                fim.removeThese.erase( unique(fim.removeThese.begin(), fim.removeThese.end()), fim.removeThese.end() );
                for(int i=0; i<fim.removeThese.size(); i++){
                    int rm = fim.removeThese[i];
                    findInVec = find(fim.activeList.begin(), fim.activeList.end(), rm);
                    if(findInVec != fim.activeList.end()){  // If it was found, remove it
                        fim.activeList.erase(findInVec);
                        meshData.is_active[rm] = false;
                    }
                }
                fim.removeThese.clear();
                
                sort(fim.addThese.begin(), fim.addThese.end());
                fim.addThese.erase( unique(fim.addThese.begin(), fim.addThese.end()), fim.addThese.end() );
                for(int i=0; i<fim.addThese.size(); i++){
                    int ad = fim.addThese[i];
                    findInVec = find(fim.activeList.begin(), fim.activeList.end(), ad);
                    if(findInVec == fim.activeList.end()){  // If it was not found, add it
                        fim.activeList.push_back(ad);
                        meshData.is_active[ad] = true;
                    }
                }
                fim.addThese.clear();
                
                sort(fim.activeList.begin(), fim.activeList.end());
                fim.activeList.erase( unique(fim.activeList.begin(), fim.activeList.end()), fim.activeList.end() );
            }

        //______________________________________________________________________
            inline void setupMemory(){

            }
            
        //______________________________________________________________________
            inline void costToGoThread(const int i, const int j){
                
            }


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            TensorFunction* tFunc;

            vector<int> seedVert;
            vector<double> seedVal;

            int threads;

            struct MeshData {
                vector<double> val;
                vector< Matrix<double> > spd;
                vector<bool> is_seed;
                vector<bool> is_active;
                vector< vector<double> > costToGo;  // costToGo[i][j] is cost to go from vertex i to its jth neighbor
            };
            MeshData meshData;

            struct FIM {
                vector<int> activeList;
                vector<int> removeThese;
                vector< vector<int> > removeTheseThread;
                vector<int> addThese;
                vector< vector<int> > addTheseThread;
                double convTol;
            };
            FIM fim;


    };

}