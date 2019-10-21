#include <cstdlib>
#include <vector>
#include <iostream>

#include "src/pch.h"

#include "src/EikonalSolution.h"
#include "src/Mesh.h"
#include "src/TensorFunction.h"
#include "src/SurfNTerp.h"
#include "src/VectorFunctionSpec.h"

#include "src/PontaniConway3dDynamics.h"
#include "src/SeparableDynamicGameBase.h"
#include "src/SeparableDynamicGame.h"
#include "src/SemiDirect.h"

#include "src/TypeErasure.h"
#include "src/VectorFunctionTypeErasure.h"

using namespace std;
using namespace WGL_DG;

template<int N>
struct MyStruct1 {
    template<typename T>
    void theFunction(T & x) const {
        x = x*N;
    }
};

template<int N>
struct MyStruct2 {
    template<typename T>
    void theFunction(T & x) const {
        x = N;
    }
    
    void thisOtherFunction() const {
        int wow = 12;
    }
};

struct MySpec {
    struct Concept {
        virtual ~Concept() = default;
        virtual void theFunction(int & x) const = 0;
        virtual void theFunction(double & x) const = 0;
    };
    template<class Holder>
    struct Model : public Holder, public virtual Concept {
        using Holder::Holder;
        virtual void theFunction(int & x) const override {
            return rubber_types::model_get(this).theFunction(x);
        }
        virtual void theFunction(double & x) const override {
            return rubber_types::model_get(this).theFunction(x);
        }
    };
    template<class Container>
    struct ExternalInterface : public Container {
        using Container::Container;
        void theFunction(int & x) const {
            return rubber_types::interface_get(this).theFunction(x);
        }
        void theFunction(double & x) const {
            return rubber_types::interface_get(this).theFunction(x);
        }
    };
};

using MyConcept = rubber_types::TypeErasure<MySpec>;



template<class Scalar>
using DgInVec = SeparableDynamicGame<0, PontaniConway3dDynamics, PontaniConway3dDynamics, Decoupled>::InputVec<Scalar>;
template<class Scalar>
using DgOutVec = SeparableDynamicGame<0, PontaniConway3dDynamics, PontaniConway3dDynamics, Decoupled>::OutputVec<Scalar>;



int main() {
    
    bool breakp;
    
    
    
    if(1){
        PontaniConway3dDynamics pc3d_1, pc3d_2;
        pc3d_1.set_mass(10);
        pc3d_1.set_mu(9.81);
        pc3d_1.set_thrust(0.001);
        pc3d_2.set_mass(20);
        pc3d_2.set_mu(9.81);
        pc3d_2.set_thrust(0.0015);
        
        PontaniConway3dDynamics::InputVec<double> x;
        x = PontaniConway3dDynamics::InputVec<double>::Random();
        
        cout << PontaniConway3dDynamics::UV << endl;
        
        cout << SeparableDynamicGame<0, PontaniConway3dDynamics, PontaniConway3dDynamics, Decoupled>::GameBase::P_XV << endl;
        SeparableDynamicGame<0, PontaniConway3dDynamics, PontaniConway3dDynamics, Decoupled> DG;
        
        DgInVec<double> dgIn = DgInVec<double>::Random();
        DgOutVec<double> dgOut;
        
        DG.compute(dgIn,dgOut);

    }
    
    if(0){
        vector<MyConcept> vec;
        MyStruct1<3> ms13;
        MyStruct1<7> ms17;
        MyStruct2<5> ms25;
        vec.push_back(ms13);
        vec.push_back(ms17);
        vec.push_back(ms25);
        
        int x1 = 12;
        double x2 = 3.14;
        vec[0].theFunction(x1);
        vec[1].theFunction(x2);
        vec[2].theFunction(x2);
        
        std::vector<VectorFunction_EigenRefCall> vf_vec;
        PontaniConway3dDynamics pc3d_1, pc3d_2;
        pc3d_1.set_mass(10);
        pc3d_1.set_mu(9.81);
        pc3d_1.set_thrust(0.001);
        pc3d_2.set_mass(20);
        pc3d_2.set_mu(9.81);
        pc3d_2.set_thrust(0.0015);
        
        vf_vec.push_back(pc3d_1);
        vf_vec.push_back(pc3d_2);
        
        PontaniConway3dDynamics::InputVec<double> x;
        x = PontaniConway3dDynamics::InputVec<double>::Random();
        Matrix<double, 8, 1> fx;
        Array<int, 6, 1> idx;
        idx.LinSpaced(6, 0, 5);
        idx[0] = 0;
        idx[1] = 2;
        idx[2] = 4;
        idx[3] = 1;
        idx[4] = 3;
        idx[5] = 5;
        cout << "Before Calculation fx:" << endl << fx << endl;
        Matrix<double, 6, 1> tempVec = fx(idx);
        cout << "Before Calculation tempVec:" << endl << tempVec << endl << endl;
        vf_vec[1].compute(x, tempVec);
        cout << "After Calculation tempVec:" << endl << tempVec << endl;
        fx(idx) = tempVec;
        cout << "After Calculation fx:" << endl << fx << endl;
    }
    
    if(0){
        
        PontaniConway3dDynamics testPursuer, testEvader;
        SeparableDynamicGame<0, PontaniConway3dDynamics, PontaniConway3dDynamics> DG;
        
        DG.gen_pursuer();
        DG.gen_evader();
        
        DG.pointer_to_evader()->set_mass(12);
        
    }
    
    if(0) {
        const int dim = 2;
    
        FastEdgesTensorFunction<dim> itf;
    
        EikonalSolution< dim, FastEdgesTensorFunction<dim>, Mesh<dim, MeshType::RecTriangular> > eikSol;
    
        for(int i=0; i<dim; i++){
            eikSol.set_bounds(i, -1, 1);
            eikSol.set_nDisc(i, 250);
        }
    
        eikSol.gen_mesh();
    
        eikSol.set_spd_func(&itf);
    
        std::vector<int> seedVert;
        std::vector<double> seedVal;
    
        seedVert.push_back(0);
        seedVal.push_back(0.0);

        seedVert.push_back(12240);
        seedVal.push_back(0.0);
        
        seedVert.push_back(55200);
        seedVal.push_back(0.0);
        
        seedVert.push_back(31375);
        seedVal.push_back(0.0);
    
        eikSol.set_seed(seedVert, seedVal);
    
        eikSol.init();
    
        cout << "Initialized..." << endl;
    
        breakp = true;
    
        eikSol.compute();
        eikSol.textFileOutput();
    }
    
    if(0){
        SurfNTerp<2, ParabolicSinkVectorFunction<2>, Mesh<2, MeshType::RecTriangular> > snt;
        
        ParabolicSinkVectorFunction<2> psvf;
        
        snt.set_vec_func(&psvf);
        
        snt.set_bounds(0, -2, 2);
        snt.set_bounds(1, -2, 2);
        snt.set_nDisc(0, 50);
        snt.set_nDisc(1, 50);
        
        snt.gen_mesh();
        
        std::vector<int> seedVert;
        std::vector<double> seedVal;
        
        seedVert.push_back(0);
        seedVal.push_back(0.0);
        
        snt.set_seed(seedVert, seedVal);
        
        snt.set_dt(4.0/100);
        
        snt.init();
        cout << "Initialized..." << endl;
        
        snt.solve();
        breakp = false;
        
        snt.textFileOutput();
    }

    return 0;
}

