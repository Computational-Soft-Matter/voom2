#include "LoopShellMesh.h"
#include "LoopShellElement.h"
#include "HalfEdgeMesh.h"

namespace voom {

  // Constructor from input file
  LoopShellMesh::LoopShellMesh(const string Nodes, const string ConnTable): Mesh(Nodes, ConnTable) {
    ifstream inp(ConnTable.c_str());
    uint NumFaces = 0, temp = 0;
    string ElType;
    // First line
    // First number is Faces, Second is type of element
    inp >> NumFaces >> ElType;
    _elements.resize(NumFaces);

    uint NumNodesEl = 0;
    // Full integration quadratic triangular element
    _quadrature = new TriQuadrature(2); 
    vector<VectorXd > QuadPoints = _quadrature->getQuadPoints();
      
    // Compute connectivity table
    vector<TriangleConnectivity> connectivities;
    for (uint e = 0; e < NumFaces; e++) {
      Vector3i triangleConn;
      int temp;
      inp >> temp >> triangleConn(0) >> triangleConn(1) >> triangleConn(2);
      //inp >> temp >> triangleConn(0) >> triangleConn(2) >> triangleConn(1);
      connectivities.push_back(triangleConn);
      //cout << triangleConn <<endl;
    }
    initializeMesh(connectivities);
    cout << "Shapes' size: "<<_shapes.size()<<endl;
    for (std::map< CornerValences,vector<LoopShellShape*> >::iterator it=_shapes.begin();
	 it!=_shapes.end(); ++it)
      {
	cout << it->first[0] << "," << it->first[1]<< "," << it->first[2] <<endl;
      }
  } 
  void LoopShellMesh::initializeMesh(vector<TriangleConnectivity> & connectivities){

    //------------------------------------------------------------------------
    // build elements using simple data structure
    std::cout << "Building elements..." << std::endl;
    clock_t t1=clock();

    HalfEdgeMesh * mesh = new HalfEdgeMesh(connectivities, _X.size());
    
    // add a boundary layer of ghost faces and vertices

    // step 1: assuming a single boundary, walk around it starting
    // from one edge

    // find a boundary edge to start
    HalfEdge * Hstart=0;
    for(int h=0; h<mesh->halfEdges.size(); h++) {
      HalfEdge * H = mesh->halfEdges[h];
      if( H->opposite == 0 ) {
	Hstart = H; 
	break;
      }
    }

    // found one, now walk around boundary, storing boundary edges in
    // order
    std::vector< HalfEdge * > boundaryEdges;
    if( Hstart !=0 ) {
      HalfEdge * H = Hstart; 
      do {
	boundaryEdges.push_back( H );
	// find next boundary edge (CCW around the boundary) by
	// walking around H's vertex CW
	H = H->next;
	while ( H->opposite != 0 ) {
	  H = H->opposite->next;
	}
      } while ( H != Hstart );
    }

    std::cout << "Identified " << boundaryEdges.size() << " boundary edges." 
	      << std::endl;

    // SD: Something is wrong, for a cylinder it detects only one boundary.
    // SD: For my periodic BC I manually added ghost nodes, so this step is not needed. 
    // SD: commented this -> if( boundaryEdges.size() > 0 ) {
    if( boundaryEdges.size() > 0 ) {
      // step 2: add a ghost node and face for each boundary edge,
      // also create ghostBC for each boundary edge
      for(int h=0; h<boundaryEdges.size(); h++) {
	HalfEdge * H = boundaryEdges[h];
	
	// add new ghost node
	int V = _X.size();
	//NodeBase::DofIndexMap idx(3,-1);
	Eigen::VectorXd gsnode(3); gsnode << 0,0,0;
	//FeNode_t * N = new FeNode_t(V, idx, X, X);
	
	
	// Create ghost triangle V0 V2 V and ghostBC for adjacent faces
	//      
	//          V
	//         / \
	//       /     \
	//     V0-------V2
	//       \  H  /
	// H->next \ /  H->prev
	//          V1
	//
	int V0 = H->vertex->id;
	int V1 = H->next->vertex->id;
	int V2 = H->prev->vertex->id;

	//Next line is needed to get the actual nodal position (see LoopGhostBC)
	gsnode = _X[V0]+_X[V2]-_X[V1];
	TriangleConnectivity c;
	c << V0, V2, V;
	_X.push_back(gsnode);
	connectivities.push_back(c);

	//LoopGhostBC * bc = new LoopGhostBC(N0,N1,N2,N);
	//_constraints.push_back( bc );
      
      }

      // step 3: add ghost faces
      // connecting neighboring ghost vertices
      for(int h=0; h<boundaryEdges.size(); h++) {
	//            
	//      VHH-------VH
	//       / \ Add / \  
	//     /     \ /     \
	//    *-------V-------*
	//        HH      H
	//
	HalfEdge * H = boundaryEdges[h];
	HalfEdge * HH = boundaryEdges[(h+1) % boundaryEdges.size()];

	int V = H->vertex->id;
	int VH = mesh->vertices.size() + h;
	int VHH = mesh->vertices.size() + (h+1) % boundaryEdges.size();

	TriangleConnectivity c;
	c << V, VH, VHH;
	connectivities.push_back(c);
      }

      // rebuild halfedge mesh
      delete mesh;
      mesh = new HalfEdgeMesh(connectivities, _X.size());
    }
       
    for(int f=0; f<mesh->faces.size(); f++) {
      Face * F = mesh->faces[f];
      // get valences of corners and find one-ring of
      // next-nearest-neighbor vertices
      vector<VectorXd> nds;
      vector<int> ndIDs;

      CornerValences v(3,0);
      //don't create element for this face if any of the vertices is on the boundary
      if( F->halfEdges[0]->vertex->boundary ) continue;
      if( F->halfEdges[1]->vertex->boundary ) continue;
      if( F->halfEdges[2]->vertex->boundary ) continue;

      // make sure that irregular node occurs at v[0]
      HalfEdge * H0 = F->halfEdges[0];
      HalfEdge * H1 = F->halfEdges[1];
      HalfEdge * H2 = F->halfEdges[2];
      int vv0 = H0->vertex->halfEdges.size();
      int vv1 = H1->vertex->halfEdges.size();
      int vv2 = H2->vertex->halfEdges.size();  
      int myidx = (vv0!=6)*0 + (vv1!=6)*1 + (vv2!=6)*2; //one of it will be nonzero for irregular patch and it will pick out that indx, esle it will be 0
      
      //int  myidx = 0; //this would be the naive way, without making sure v[0] is irregular
      int idv[3]={myidx,((myidx+1) % 3),((myidx+2) % 3)};
      //cout << idv[0] << " " << idv[1] << " " << idv[2] << endl;
      for(int a=0; a<3; a++) {
	
	  // add corner vertices
	  HalfEdge * H = F->halfEdges[idv[a]];
	  int Id = H->vertex->id;
	  ndIDs.push_back( Id );
	  //cout << Id << ", "; 
	  nds.push_back( _X[Id] );
	  // valence is the number of incident halfedges
	  v[a] = H->vertex->halfEdges.size();
	}
      //cout << v[0] << " " << v[1] << " " << v[2] << endl;
      //cout << endl;
      // walk around the one-ring adding vertices
      //
      // Start by walking CCW through the vertices incident to the
      // second corner vertex.  Stop at first vertex incident to the
      // third corner vertex.
      //cout<< F->halfEdges[1]->vertex->id<<": \n    ";
      int v1=3;
      for(HalfEdge * H = F->halfEdges[idv[1]]->opposite->next; 
    	  H->vertex->id != F->halfEdges[idv[2]]->opposite->next->vertex->id; 
    	  H = H->next->opposite->next  ) {
	int Id = H->vertex->id;
	ndIDs.push_back( Id );
	//cout << Id << ", "; 
    	nds.push_back( _X[Id] );
    	v1++;
      }
      //cout << endl;
      // Now continue CCW through the vertices incident to the third
      // corner vertex. Stop at first vertex incident to the first
      // corner vertex.
      int v2=3;
      //cout<< F->halfEdges[2]->vertex->id<<": \n    ";
      for(HalfEdge * H = F->halfEdges[idv[2]]->opposite->next; 
    	  H->vertex->id != F->halfEdges[idv[0]]->opposite->next->vertex->id; 
    	  H = H->next->opposite->next  ) {
	int Id = H->vertex->id;
	ndIDs.push_back( Id );
	//cout << Id << ", "; 
    	nds.push_back( _X[Id] );
    	v2++;
      }
      //cout << endl;
      // Now continue CCW through the vertices incident to the first
      // corner vertex. Stop at first vertex incident to the second
      // corner vertex.
      //cout<< F->halfEdges[0]->vertex->id<<": \n    ";
      int v0=3;
      for(HalfEdge * H = F->halfEdges[idv[0]]->opposite->next; 
    	  H->vertex->id != F->halfEdges[idv[1]]->opposite->next->vertex->id; 
    	  H = H->next->opposite->next  ) {
	int Id = H->vertex->id;
	ndIDs.push_back( Id );
	//cout << Id << ", "; 
    	nds.push_back( _X[Id] );
    	v0++;
      }
      //cout << endl << "--------------------" << endl;
      assert(v0==v[0]);
      assert(v1==v[1]);
      assert(v2==v[2]);

      //const unsigned npn = v(0) + v(1) + v(2) - 6;
      
      vector<VectorXd> QuadPoints = _quadrature->getQuadPoints();
      map< CornerValences, vector<LoopShellShape*> >::iterator it = _shapes.find(v);
      if ( it != _shapes.end()) { //If _shapes contains key 'V'
      	_elements[f] = new LoopShellElement(f, ndIDs, nds, _shapes[v] , _quadrature, v);
      }
      else{
      	for (uint q = 0; q < QuadPoints.size(); q++ ){
      	  LoopShellShape* shape = new LoopShellShape(nds.size(), v, QuadPoints[q]);
      	  _shapes[v].push_back( shape );
      	}
      	_elements[f] = new LoopShellElement(f, ndIDs, nds, _shapes[v] , _quadrature, v);
      }
      
    }
    delete mesh;

    clock_t t2=clock();
    std::cout << "Done building elements.  Elapsed time: "
    	      << (double)(t2-t1)/CLOCKS_PER_SEC
    	      << std::endl;
    std::cout<<"\nNumber of elements "<<_elements.size()<<std::endl;
    // compute mechanics

  }
  
} // namespace voom
