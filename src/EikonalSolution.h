#pragma once

#include <vector>
#include <limits>
#include <Eigen/Core>

#include "ctpl.h"
#include "Mesh.h"
#include "TensorFunction.h"

using namespace std;

namespace WGL_DG {

    template<int dim, class TensorFunction, class Mesh>
    struct EikonalSolution {

        public:
            /// Typedefs


        //======================================================================
            /// Constructors
            EikonalSolution(){
                vector<double> MeshData::val;
                vector< Eigen::Matrix<double, dim, dim> > MeshData::spd;
                vector<bool> MeshData::is_seed;

                vector<int> FIM::activeList;
                vector<int> FIM::removeThese;
                vector<int> FIM::addThese;
            }


        //======================================================================
            /// Setup
            inline void set_mesh(Mesh *m){
                this->mesh = m;

                MeshData::val.resize(mesh->nVert);
                MeshData::spd.resize(mesh->nVert);
                MeshData::is_seed.resize(mesh->nVert);
            }

        //______________________________________________________________________
            inline void set_spd_func(TensorFunction *tf){
                this->tFunc = tf;
            }

        //______________________________________________________________________
            inline void set_seed(vector<int> verts, vector<double> vals){
                fill(MeshData::val.begin(), MeshData::val.end(), std::numeric_limits<double>::infinity());
                fill(MeshData::is_seed.begin(), MeshData::is_seed.end(), false);
                for(int i=0; i<verts.size(); i++){
                    MeshData::val[verts[i]] = vals[i];
                    MeshData::is_seed[verts[i]] = true;
                }
            }


        //======================================================================
            /// Primary Methods
            inline void init(){

            }

        //______________________________________________________________________
            inline void compute(){

            }


        //======================================================================
            /// Properties


    ////////////////////////////////////////////////////////////////////////////
        protected:
            /// Intermediate Methods


    ////////////////////////////////////////////////////////////////////////////
        private:
            /// Properties
            Mesh* mesh;
            TensorFunction* tFunc;

            struct MeshData {
                static vector<double> val;
                static vector< Eigen::Matrix<double, dim, dim> > spd;
                static vector<bool> is_seed;
            };

            struct FIM {
                static vector<int> activeList;
                static vector<int> removeThese;
                static vector<int> addThese;
            }


    };

}