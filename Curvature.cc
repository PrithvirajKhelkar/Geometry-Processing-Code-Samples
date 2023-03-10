#include <math.h>
#include <algorithm>
// #include "Curvature.hh"


void Curvature::show_normal(const int _normal_type) {
    switch(_normal_type) {
        case 0: compute_normals_with_constant_weights();
            break;
        case 1: compute_normals_by_area_weights();
            break;
        default: compute_normals_with_angle_weights();
    }
}

void Curvature::show_curvature(const int _curvature_type) {
    switch(_curvature_type) {
        case 0: calc_uniform_laplacian();
            break;
        case 1: calc_mean_curvature();
            break;
        default: calc_gauss_curvature();
    }

//    auto values = mesh_.property(vertex_curvature_).data_vector();
    color_coding(vertex_curvature_);
}

void Curvature::show_valence() {
    for(auto vh : mesh_.vertices())
        mesh_.property(vertex_valence_, vh) = mesh_.valence(vh);

//    auto values = mesh_.property(vertex_valence_).data_vector();
    color_coding(vertex_valence_, 3, 8);
}

Curvature::Point normal(Curvature::Point p1, Curvature::Point p2) {
    /*
     Function for the normal.
    */
    return (p1.cross(p2)).normalize();
}

float constant_weight(Curvature::Point p1, Curvature::Point p2) {
    /*
     Function for the constant weight.
    */
    return 1.0;
}

float area_weight(Curvature::Point p1, Curvature::Point p2) {
    /*
     Function for the area weight of the triangle between the two vectors.
    */
    return abs((p1.cross(p2)).norm())/2;
}

float angle_weight(Curvature::Point p1, Curvature::Point p2) {
    /*
     Function for the angle weight between the two vectors.
    */
    return acos(p1.dot(p2)/p1.norm()*p2.norm());
}

Curvature::Operator<Curvature::Point> uniform_laplace {
    /*
        Operator to apply the uniform laplace approximation
    */
    [](TriMesh& mesh, TriMesh::VertexHandle vh_0, TriMesh::VertexHandle vh_1, TriMesh::VertexHandle vh_2, double weight) -> Curvature::Point {return mesh.point(vh_1) - mesh.point(vh_0);},
    [](TriMesh& mesh, Curvature::Point vertex_normal, TriMesh::VertexHandle vh, double weight, Curvature::Point p) -> double {return vertex_normal.dot(p)/(mesh.valence(vh)*2);}
};

Curvature::Operator<Curvature::Point> laplace_beltrami {
    /*
        Laplace Beltrami Operator
    */
    [](TriMesh& mesh, TriMesh::VertexHandle vh_0, TriMesh::VertexHandle vh_1, TriMesh::VertexHandle vh_2, double weight) -> Curvature::Point {return weight*(mesh.point(vh_1) - mesh.point(vh_0));},
    [](TriMesh& mesh, Curvature::Point vertex_normal, TriMesh::VertexHandle vh, double weight, Curvature::Point p) -> double {return vertex_normal.dot(p)*weight/2;}
};

Curvature::Operator<double> gaussian {
    /*
        Gaussian Curvature
    */
    [](TriMesh& mesh, TriMesh::VertexHandle vh_0, TriMesh::VertexHandle vh_1, TriMesh::VertexHandle vh_2, double weight) -> double {return acos(((mesh.point(vh_1) - mesh.point(vh_0)).normalize()).dot((mesh.point(vh_2) - mesh.point(vh_0)).normalize()));},
    [](TriMesh& mesh, Curvature::Point vertex_normal, TriMesh::VertexHandle vh, double weight, double value) -> double {return ((2*M_PI) - value)*weight*2;}
};

template <typename T>
T Curvature::aggregate_over_neighborhood(TriMesh::VertexHandle vh, std::function<T(TriMesh::VertexHandle, TriMesh::VertexHandle)> func) {
    /*
        Function to interate over the ring neighborhood of a vertex and apply a computation func
    */
        
    T value = T(0);

    Mesh::VertexFaceIter vf_it = mesh_.vf_iter(vh);

    while (vf_it.is_valid()) {
        TriMesh::FaceVertexIter fv_it = mesh_.fv_iter(*vf_it);
        TriMesh::VertexHandle v0 = *fv_it;
        TriMesh::VertexHandle v1 = *(++fv_it);
        TriMesh::VertexHandle v2 = *(++fv_it);
        
        if(v0 == vh) {
            value += func(v1, v2);
        } else if (v1 == vh) {
            value += func(v2, v0);
        } else {
            value += func(v0, v1);
        }
        ++vf_it;
    }

    return value;
}


std::function<void(TriMesh::VertexHandle vh)> Curvature::calculate_vertex_normal_by_weight_function(float (*weight)(Point, Point)){
    /*
        Function to return a function that calculates the summation of the weighted normals of all the faces around a vertex, and stores its normalized value as the vertex normal.
    */
    return [this, weight](TriMesh::VertexHandle vh){
        Point p0 = mesh_.point(vh);
        Point vertex_normal = aggregate_over_neighborhood<Point>(vh, [this, p0, weight](TriMesh::VertexHandle vh_1, TriMesh::VertexHandle vh_2) -> Curvature::Point {return weight(mesh_.point(vh_1) - p0, mesh_.point(vh_2) - p0) * normal(mesh_.point(vh_1) - p0, mesh_.point(vh_2) - p0);});
        mesh_.property(vertex_normal_, vh) = vertex_normal.normalize();
    };
}

template <typename T>
std::function<void(TriMesh::VertexHandle vh)> Curvature::calculate_vertex_curvature_by_curvature_operator(Operator<T> op){
    /*
    A template function for calculating curvature using any operator.
    */
    return [this, op](TriMesh::VertexHandle vh){
        if(!mesh_.is_boundary(vh)){
            Point p0 = mesh_.point(vh);
            Point vertex_normal = mesh_.property(vertex_normal_, vh);
            double vertex_weight = mesh_.property(vertex_weight_, vh);

            double vertex_curvature = op.outer(mesh_, vertex_normal, vh, vertex_weight, aggregate_over_neighborhood<T>(vh, [this, vh, op](TriMesh::VertexHandle vh_1, TriMesh::VertexHandle vh_2) -> T {
                    double edge_weight = mesh_.property(edge_weight_, mesh_.edge_handle(mesh_.find_halfedge(vh_1, vh)));
                    return op.inner(mesh_, vh, vh_1, vh_2,  edge_weight);
                }));
            mesh_.property(vertex_curvature_, vh) = vertex_curvature;

            if ( vertex_curvature > max_curvature_){
                max_curvature_ = vertex_curvature;
            }

            if ( vertex_curvature < min_curvature_){
                min_curvature_ = vertex_curvature;
            }
        }
    };
}

void Curvature::iterate_over_vertices(std::function<void(TriMesh::VertexHandle)> func){
    /*
        Function to iterate over all vertices in the mesh and perform the specified function on each vertex.
    */
    for (TriMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
        func(*v_it);
    }
}



//====================================================================================================================//
// ========================================================================
// EXERCISE 1.1
// ========================================================================
void Curvature::compute_normals_with_constant_weights() {
    iterate_over_vertices(calculate_vertex_normal_by_weight_function(&constant_weight));
    update_normals();
}
// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Curvature::compute_normals_by_area_weights() {
    iterate_over_vertices(calculate_vertex_normal_by_weight_function(&area_weight));
    update_normals();
}
// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Curvature::compute_normals_with_angle_weights() {
    iterate_over_vertices(calculate_vertex_normal_by_weight_function(&angle_weight));
    update_normals();
}
// ========================================================================
// EXERCISE 2.1
// ========================================================================
void Curvature::calc_uniform_laplacian() {
    min_curvature_ = DBL_MAX;
    max_curvature_ = -DBL_MAX;
    compute_normals_with_constant_weights();
    iterate_over_vertices(calculate_vertex_curvature_by_curvature_operator(uniform_laplace));

    std::cout<<"Min Uniform Laplace value is: " << min_curvature_ << std::endl;
    std::cout<<"Max Uniform Laplace value is: " << max_curvature_ << std::endl;
}
// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Curvature::calc_mean_curvature() {
    min_curvature_ = DBL_MAX;
    max_curvature_ = -DBL_MAX;

    calc_weights();
    compute_normals_with_constant_weights();
    iterate_over_vertices(calculate_vertex_curvature_by_curvature_operator(laplace_beltrami));

    std::cout<<"Min Laplace-Beltrami value is: " << min_curvature_ << std::endl;
    std::cout<<"Max Laplace-Beltrami value is: " << max_curvature_ << std::endl;
}
// ========================================================================
// EXERCISE 2.3
// ========================================================================
void Curvature::calc_gauss_curvature() {
    min_curvature_ = DBL_MAX;
    max_curvature_ = -DBL_MAX;

    calc_weights();
    compute_normals_with_constant_weights();
    iterate_over_vertices(calculate_vertex_curvature_by_curvature_operator(gaussian));

    std::cout<<"Min Gauss Curvature value is: " << min_curvature_ << std::endl;
    std::cout<<"Max Gauss Curvature value is: " << max_curvature_ << std::endl;
}

//====================================================================================================================//
void Curvature::calc_weights() {
    calc_vertices_weights();
    calc_edges_weights();
}

void Curvature::calc_vertices_weights() {
    double area = 0.;
    for (auto vh: mesh_.vertices()) {
        area = 0.0;

        for(auto vih_it = mesh_.vih_iter(vh); vih_it.is_valid(); ++vih_it) {
            if(mesh_.is_boundary(*vih_it))
                continue;

            area += mesh_.calc_sector_area(*vih_it) * 0.3333f;
        }

        mesh_.property(vertex_weight_, vh) = 0.5 / area;
    }
}

void Curvature::calc_edges_weights() {
    OpenMesh::HalfedgeHandle heh0, heh1, heh2;
    Point p0(0.), p1(0.), p2(0.), d0(0.), d1(0.);
    double cross;
    double w;
    for (auto eh: mesh_.edges()) {
        w = 0.0;

        heh0 = mesh_.halfedge_handle(eh, 0);
        p0 = mesh_.point(mesh_.to_vertex_handle(heh0));

        heh1 = mesh_.halfedge_handle(eh, 1);
        p1 = mesh_.point(mesh_.to_vertex_handle(heh1));

        if (!mesh_.is_boundary(heh0))
        {
            heh2 = mesh_.next_halfedge_handle(heh0);
            p2 = mesh_.point(mesh_.to_vertex_handle(heh2));
            d0 = p0 - p2;
            d1 = p1 - p2;
            cross = std::max(1e-16, (d0 % d1).norm());
            w += (d0|d1) / cross;
        }

        if (!mesh_.is_boundary(heh1))
        {
            heh2 = mesh_.next_halfedge_handle(heh1);
            p2 = mesh_.point(mesh_.to_vertex_handle(heh2));
            d0 = p0 - p2;
            d1 = p1 - p2;
            cross = std::max(1e-16, (d0 % d1).norm());
            w += (d0|d1) / cross;
        }

        w = std::max(0., w);
        mesh_.property(edge_weight_, eh) = w;
    }
}

void Curvature::update_normals() {
    for(auto vh : mesh_.vertices())
        mesh_.set_normal(vh, 0.5*avr_e_length_*mesh_.property(vertex_normal_, vh));
}

template <typename T>
void Curvature::color_coding(const OpenMesh::VPropHandleT<T>& _vprop, const double _min_value, const double _max_value, const int _bound) {
    auto values = mesh_.property(_vprop).data_vector();
    auto min_value = _min_value;
    auto max_value = _max_value;

    if(min_value == 0 && max_value == 0) {
        // discard upper and lower bound
        auto n = values.size()-1;
        auto i = n / _bound;

        std::sort(values.begin(), values.end());
        min_value = values[i];
        max_value = values[n-1-i];
    }
//    std::cerr<<"\nmax1 "<<max_value<<" "<<min_value<<" ";

    const auto range = max_value - min_value;
    color_.set_range(0, 1.0, false);

    for(auto vh : mesh_.vertices()) {
        auto t = (mesh_.property(_vprop, vh) - min_value)/range;
        mesh_.set_color(vh, color_.color_float4(t));
    }
}
