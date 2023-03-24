#include "Smoothing.hh"

double constant_edge_weight(TriMesh& mesh_, OpenMesh::EPropHandleT<double> edge_weight_, TriMesh::EdgeHandle eh) {
    /*
     Function for the constant weight.
    */
    return 1.0;
}

double cotan_edge_weight(TriMesh& mesh_, OpenMesh::EPropHandleT<double> edge_weight_, TriMesh::EdgeHandle eh) {
    /*
     Function for the constant weight.
    */
    return mesh_.property(edge_weight_, eh);
}

double area_vertex_weight(TriMesh& mesh_, OpenMesh::VPropHandleT<double> vertex_weight_, TriMesh::VertexHandle vh) {
    /*
     Function for the constant weight.
    */
    return mesh_.property(vertex_weight_, vh);;
}

Eigen::SparseMatrix<double> Smoothing::get_M_matrix(double (*edge_weight)(TriMesh& mesh_, OpenMesh::EPropHandleT<double> edge_weight_, TriMesh::EdgeHandle eh)){
   Eigen::SparseMatrix<double> M(mesh_.n_vertices(), mesh_.n_vertices());

   std::vector<Eigen::Triplet<double>> tripletList2;
   tripletList2.reserve(7*mesh_.n_vertices());

   for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
        double total = 0.0;

        TriMesh::VertexFaceIter vf_it = mesh_.vf_iter(*v_it);
        TriMesh::VertexHandle next;

        while (vf_it.is_valid()) {
            TriMesh::FaceVertexIter fv_it = mesh_.fv_iter(*vf_it);
            TriMesh::VertexHandle v0 = *fv_it;
            TriMesh::VertexHandle v1 = *(++fv_it);
            TriMesh::VertexHandle v2 = *(++fv_it);
            
            if(v0 == *v_it) {
                next = v1;
            } else if (v1 == *v_it) {
                next = v2;
            } else {
                next = v0;
            }
            
            TriMesh::HalfedgeHandle heh = mesh_.find_halfedge(*v_it, next);
            if (heh.is_valid())
            {
                TriMesh::EdgeHandle eh = mesh_.edge_handle(heh);
                double weight = edge_weight(mesh_, edge_weight_, eh);
                total += weight;
                tripletList2.push_back(Eigen::Triplet<double>((*v_it).idx(), next.idx(), weight));
            }
            
            ++vf_it;
        }
        tripletList2.push_back(Eigen::Triplet<double>((*v_it).idx(), (*v_it).idx(), -total));
    }
    M.setFromTriplets(tripletList2.begin(), tripletList2.end());

    return M;
}

Eigen::SparseMatrix<double> Smoothing::get_D_matrix(double (*vertex_weight)(TriMesh& mesh_, OpenMesh::VPropHandleT<double> vertex_weight_, TriMesh::VertexHandle vh)){
    Eigen::SparseMatrix<double> D(mesh_.n_vertices(), mesh_.n_vertices());

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mesh_.n_vertices());
    for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
        tripletList.push_back(Eigen::Triplet<double>((*v_it).idx(), (*v_it).idx(), vertex_weight(mesh_, vertex_weight_, *v_it)));
    }
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    return D;
}

Eigen::SparseMatrix<double> Smoothing::get_inverse_D_matrix(double (*vertex_weight)(TriMesh& mesh_, OpenMesh::VPropHandleT<double> vertex_weight_, TriMesh::VertexHandle vh)){
   Eigen::SparseMatrix<double> D(mesh_.n_vertices(), mesh_.n_vertices());

   std::vector<Eigen::Triplet<double>> tripletList(mesh_.n_vertices());

    for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
        tripletList.push_back(Eigen::Triplet<double>((*v_it).idx(), (*v_it).idx(), 1/vertex_weight(mesh_, vertex_weight_, *v_it)));
    }
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    return D;
}

Eigen::MatrixXd Smoothing::get_vertices_matrix(){
    Eigen::MatrixXd vertices(mesh_.n_vertices(), 3);

    // Loop over all vertices and extract their (x, y, z) coordinates
    for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it)
    {
        Smoothing::Point p = mesh_.point(*v_it);
        vertices.row(v_it->idx()) << p[0], p[1], p[2];
    }

    return vertices;
}

void Smoothing::set_vertices(Eigen::MatrixXd vertices){
    for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it)
    {
        if(!mesh_.is_boundary(*v_it)){
            mesh_.set_point(*v_it, Smoothing::Point(vertices(v_it->idx(), 0), vertices(v_it->idx(), 1), vertices(v_it->idx(), 2)));
        }
    }
}

Eigen::MatrixXd Smoothing::explicit_smoothing(int _iterations, double constant, double (*edge_weight)(TriMesh&, OpenMesh::EPropHandleT<double>, TriMesh::EdgeHandle)) {
    Eigen::MatrixXd vertices = get_vertices_matrix();
    for (int iter = 0; iter < _iterations; ++iter) {
        calc_edges_weights();
        
        // TODO parallel
        for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it){
            TriMesh::VertexHandle v_i = *v_it;
            if(!mesh_.is_boundary(v_i)) {
                double total_weight = 0;
                Smoothing::Point p = Smoothing::Smoothing::Point(0);
                for (TriMesh::VertexVertexIter vv_it=mesh_.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
                    TriMesh::VertexHandle v_j = *vv_it; 
                    TriMesh::HalfedgeHandle heh = mesh_.find_halfedge(v_i, v_j);
                    if (heh.is_valid())
                    {
                        TriMesh::EdgeHandle eh = mesh_.edge_handle(heh);
                        double weight = edge_weight(mesh_, edge_weight_, eh);
                        total_weight += weight;
                        
                        p += weight*(mesh_.point(v_j) - mesh_.point(v_i));
                    }
                }
                if (total_weight != 0) {
                    p *= constant/total_weight;
                    p += mesh_.point(v_i);  
                    vertices.row(v_it->idx()) << p[0], p[1], p[2];
                }
            }
        }
        set_vertices(vertices);
    }
    mesh_.update_normals();
    return vertices;
}

void Smoothing::implicit_smoothing(double timestep, double diffuse, double (*edge_weight)(TriMesh& mesh_, OpenMesh::EPropHandleT<double> edge_weight_, TriMesh::EdgeHandle eh), double (*vertex_weight)(TriMesh& mesh_, OpenMesh::VPropHandleT<double> vertex_weight_, TriMesh::VertexHandle vh)){
    calc_weights();
    Eigen::SparseMatrix<double> M = get_M_matrix(edge_weight);
    Eigen::SparseMatrix<double> D_inv = get_inverse_D_matrix(vertex_weight);
    Eigen::MatrixXd vertices = get_vertices_matrix();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(D_inv - timestep*diffuse*M);
    vertices = solver.solve(D_inv*vertices);

    set_vertices(vertices);

    mesh_.update_normals();

}

void Smoothing::feature_enhancement(int iterations, int coeff, double diffuse, double (*edge_weight)(TriMesh& mesh_, OpenMesh::EPropHandleT<double> edge_weight_, TriMesh::EdgeHandle eh)) {
    Eigen::MatrixXd p_in = get_vertices_matrix();
    Eigen::MatrixXd p_out = explicit_smoothing(iterations, diffuse, edge_weight);
    p_out += coeff*(p_in - p_out);
    set_vertices(p_out);
    mesh_.update_normals();
}


// ======================================================================
// EXERCISE 1.1
// ========================================================================
void Smoothing::uniform_smooth(const int _iterations) {
    
    explicit_smoothing(_iterations, 0.5, &constant_edge_weight);
}

// ======================================================================
// EXERCISE 1.2
// ========================================================================
void Smoothing::cotan_laplacian_smooth(const int _iterations) {

    explicit_smoothing(_iterations, 0.5, &cotan_edge_weight);
}

// ======================================================================
// EXERCISE 2
// ========================================================================
void Smoothing::implicit_smooth(const double _timestep) {
    implicit_smoothing(_timestep, 0.5, &cotan_edge_weight, &area_vertex_weight);
}

// ======================================================================
// EXERCISE 3.1
// ========================================================================
void Smoothing::uniform_laplacian_enhance_feature(const int _iterations, const int _coefficient) {
    feature_enhancement(_iterations, _coefficient, 0.5, &constant_edge_weight);
}

// ======================================================================
// EXERCISE 3.2
// ========================================================================
void Smoothing::cotan_laplacian_enhance_feature(const int _iterations, const int _coefficient) {
    feature_enhancement(_iterations, _coefficient, 0.5, &cotan_edge_weight);
}