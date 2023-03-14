

void insert_point(TriMeshObject * _tri_obj, const bool _on_edge, const TriMesh::FaceHandle& _fh,
    /*
        Hidden code to add the vertex in the mesh
    */

   
    std::vector<TriMesh::EdgeHandle> edges;
    for (TriMesh::VertexFaceIter vf_it = mesh->vf_iter(vh); vf_it.is_valid(); ++vf_it) {
    	for (TriMesh::FaceHalfedgeIter fh_it = mesh->fh_iter(*vf_it); fh_it.is_valid(); ++ fh_it) {
		if (mesh->to_vertex_handle(*fh_it) != vh) {
			TriMesh::EdgeHandle eh = mesh->edge_handle(*fh_it);
			edges.push_back(eh);
			}
	}
    }

    for(TriMesh::EdgeHandle edge : edges) {
    	if(!is_delaunay(mesh, edge) && mesh->is_flip_ok(edge)){
		mesh->flip(edge);
	}
    }

    mesh->garbage_collection();
}

TriMesh::Point get_point_on_paraboloid(TriMesh::Point point) {
	double x = point[0];
	double y = point[1];
	double z = x*x + y*y;
	return TriMesh::Point(x, y, z);
}


bool is_same_side_of_plane(TriMesh::Point p1, TriMesh::Point p2, TriMesh::Point p3, TriMesh::Point p4){
	TriMesh::Point v1 = p2 - p1;
	TriMesh::Point v2 = p3 - p1;
	TriMesh::Point normal = v2.cross(v1);

	double signed_distance = normal.dot(p4-p1);

	return signed_distance > 0;
}

bool is_delaunay(TriMesh * _mesh, TriMesh::EdgeHandle _eh) {

    if (_mesh == 0) {
        return false;
    }
    bool result = true;
    // Get half edge corresponding to the edge _eh
    TriMesh::HalfedgeHandle heh = _mesh->halfedge_handle(_eh, 0);
	
    // Get the vertices for the adjacent triangle
    TriMesh::VertexHandle v1 = _mesh->to_vertex_handle(heh);
    TriMesh::VertexHandle v2 = _mesh->to_vertex_handle(_mesh->next_halfedge_handle(heh));
    TriMesh::VertexHandle v3 = _mesh->to_vertex_handle(_mesh->next_halfedge_handle(_mesh->opposite_halfedge_handle(heh)));

    // Get the other vertex
    heh = _mesh->opposite_halfedge_handle(heh);
    TriMesh::VertexHandle v4 = _mesh->to_vertex_handle(heh);

    // Get points
    TriMesh::Point p1 = _mesh->point(v1);
    TriMesh::Point p2 = _mesh->point(v2);
    TriMesh::Point p3 = _mesh->point(v3);
    TriMesh::Point p4 = _mesh->point(v4);
	
    std::cout<<"p1: "<<p1<<", p2: "<<p2<<", p3: "<<p3<<", p4: "<<p4<<std::endl;
    // Get points on paraboloid
    TriMesh::Point pp1 = get_point_on_paraboloid(p1);
    TriMesh::Point pp2 = get_point_on_paraboloid(p2);
    TriMesh::Point pp3 = get_point_on_paraboloid(p3);
    TriMesh::Point pp4 = get_point_on_paraboloid(p4);

    std::cout<<"pp1: "<<pp1<<", pp2: "<<pp2<<", pp3: "<<pp3<<", pp4: "<<pp4<<std::endl;

    return is_same_side_of_plane(pp1, pp2, pp3, pp4) && is_same_side_of_plane(pp2, pp1, pp4, pp3);
}
