template<int dim>
class EikSol {

    vector< Eigen::Matrix<double,dim,1> > meshVerts;  // List of all vertices in mesh
    vector< vector<int> > meshConn;  // 

    struct MeshData {
        static vector< vector<int> > neigh;
        static vector< vector<double> > val;
        static vector< Eigen::Matrix<double, dim, dim> > spd;
    }

    vector<int> actList;
    vector<int> removeThese;
    vector<int> addThese;

    void mainAlg(){
        while(actList.size()){
            for(int i=0; i<actList.size(); i++){  // This is parallelizable
                double p = MeshData::val[actList[i]];
                double q = localApprox(actList[i]);  // Godunov or whatever
                if(p > q){
                    MeshData::val[actList[i]] = q;  // 
                }
                if(abs(p-q) < tol){
                    for(int j=0; j<MeshData::neigh[actList[i]].size(); j++){  // This is also parallelizable
                        int nb = MeshData::neigh[actList[i]][j];
                        if(nb IS NOT IN actList){  // The jth neighbor of vertex actList[i]
                            double jp = MeshData::val[nb];  // Potential memory access conflict if vertex nb is the neigbor of two points in actList
                            double jq = localApprox(nb);
                            if(jp > jq){
                                MeshData::val[nb] = jq;
                                // Queued until the next while pass
                                addThese.push_back(nb);
                            }
                        }
                    }
                    // remove actList[i] from actList. Issues with parallel vector resizing? Queue index for later deletion?
                    removeThese.push_back(actList[i]);
                }
            }
            // remove vertices listed in removeThese
            // add vertices listed in addThese
            updateActiveList();
        }
    }

    inline double localApprox(int v, vector< Eigen::Matrix<double, dim, 1> > nbPos, vector<double> nbVals){

    }

    inline double solvePDE(int v, ){
        Eigen::Matrix<double, dim, dim> M = MeshData::spd[v];
        Eigen::Matrix<double, dim, 1> x = meshVerts[v];
        vector<double> sols;
        for(int n=0; n<MeshData::neigh[v].size(); n++){  // Could do this in parallel, but probably not necessary?
            int nb = MeshData::neigh[v][n];
            if(MeshData::val[nb] == std::numeric_limits<double>::infinity()){
                wow fake;
            }else{
                Eigen::Matrix<double, dim, 1> y = meshVerts[nb];
                double uy = MeshData::val[nb];

                Eigen::Matrix<double, dim, 1> del = (x-y).cwiseInverse();
                double rad = del.transpose()*M*del;
                sols.push_back(1/sqrt(rad) + uy);
            }
        }
        return std::min_element(sols.begin(), sols.end());
    }

};