#include <cstdlib>
#include <vector>
#include <iostream>

#include "src/EikonalSolution.h"
#include "src/Mesh.h"
#include "src/TensorFunction.h"
#include "src/SurfNTerp.h"

using namespace std;
using namespace WGL_DG;

int main() {
    
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
    
    bool breakp = true;
    
    eikSol.compute();
    eikSol.textFileOutput();

    return 0;
}

