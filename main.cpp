#include <cstdlib>
#include <vector>
#include <iostream>

#include "src/EikonalSolution.h"
#include "src/Mesh.h"
#include "src/TensorFunction.h"

using namespace std;
using namespace WGL_DG;

int main() {
    
    vector<int> addThese = {1,2,3,2,4,1,3,5,7,1,6};
    sort(addThese.begin(), addThese.end());
    addThese.erase( unique(addThese.begin(), addThese.end()), addThese.end() );

    for(int i=0; i<addThese.size(); i++){
      cout << addThese[i] << '\n';
    }
    
    return -1;
    
    const int dim = 2;
    
    FastEdgesTensorFunction<dim> itf;
    
    EikonalSolution< dim, FastEdgesTensorFunction<dim>, Mesh<dim, MeshType::Regular> > eikSol;
    
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

    seedVert.push_back(62400);
    seedVal.push_back(0.0);
    
    eikSol.set_seed(seedVert, seedVal);
    
    eikSol.init();
    
    cout << "Initialized..." << endl;
    
    bool breakp = true;
    
    eikSol.compute();
    eikSol.textFileOutput();

    return 0;
}

