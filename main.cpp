#include <cstdlib>

#include "src/EikonalSolution.h"
#include "src/Mesh.h"

using namespace std;
using namespace WGL_DG;

int main() {
    
    Mesh<4, MeshType::Regular> m;
    EikonalSolution< 4, TensorFunction<4>, Mesh<4, MeshType::Regular> > eikSol;
    
    m.set_bounds(0, -1, 1);
    m.set_bounds(1, -1, 1);
    m.set_bounds(2, -1, 1);
    m.set_bounds(3, -1, 1);
    
    m.set_nDisc(0, 10);
    m.set_nDisc(1, 10);
    m.set_nDisc(2, 10);
    m.set_nDisc(3, 10);
    
    m.gen_mesh();

    return 0;
}

