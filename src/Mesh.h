#pragma once

#include <vector>
#include <Eigen/Core>

namespace WGL_DG {

    enum MeshType {
        Regular,
    }

    /// Template
    template<int dim, int meshType>
    struct Mesh {
        
    };

////////////////////////////////////////////////////////////////////////////////
    /// Specializations
    template<int dim>
    struct Mesh<dim, MeshType::Regular> {

    }

}