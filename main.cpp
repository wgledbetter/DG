#include <cstdlib>
#include <vector>
#include <iostream>

#include "src/EikonalSolution.h"
#include "src/Mesh.h"
#include "src/TensorFunction.h"
#include "src/SurfNTerp.h"
#include "src/VectorFunction.h"

using namespace std;
using namespace WGL_DG;

int main() {
    
    bool breakp;
    
    if(0) {
        const int dim = 2;
    
        FastEdgesTensorFunction<dim> itf;
    
        EikonalSolution< dim, FastEdgesTensorFunction<dim>, Mesh<dim, MeshType::Regular> > eikSol;
    
        for(int i=0; i<dim; i++){
            eikSol.set_bounds(i, -1, 1);
            eikSol.set_nDisc(i, 50);
        }
    
        eikSol.gen_mesh();
    
        eikSol.set_spd_func(&itf);
    
        std::vector<int> seedVert;
        std::vector<double> seedVal;
    
        seedVert.push_back(0);
        seedVal.push_back(0.0);

        seedVert.push_back(624);
        seedVal.push_back(0.0);
    
        eikSol.set_seed(seedVert, seedVal);
    
        eikSol.init();
    
        cout << "Initialized..." << endl;
    
        breakp = true;
    
        eikSol.compute();
        eikSol.textFileOutput();
    }
    
    if(1){
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

