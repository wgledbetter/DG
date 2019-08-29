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
            for(int i=0; i<actList.size(); i++){  // This is parallellizable
                double p = MeshData::val[actList[i]];
                double q = localApprox(actList[i]);  // Godunov or whatever
                if(p > q){
                    MeshData::val[actList[i]] = q;  // 
                }
                if(abs(p-q) < tol){
                    for(int j=0; j<MeshData::neigh[actList[i]].size(); j++){  // For each neighbor of actList[i]
                        if(MeshData::neigh[actList[i]][j] IS NOT IN actList){  // The jth neighbor of vertex actList[i]
                            // This also maybe should be queued until the next while pass
                            addThese.push_back(MeshData::neigh[actList[i]][j])
                        }
                    }
                    // remove actList[i] from actList. Issues with parallel vector resizing? Queue index for later deletion?
                    removeThese.push_back(actList[i]);
                }
            }
            // remove vertices listed in removeThese
            // add vertices listed in addThese
        }
    }

};