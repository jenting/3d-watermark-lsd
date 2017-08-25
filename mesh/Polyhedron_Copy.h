/***************************************************************************
copy_poly.h  -
----------------------------------------------------------------------------
begin                : march 2006
copyright            : (C) 2006 by Celine Roudet - Liris
email                : croudet@liris.cnrs.fr
***************************************************************************/

#ifndef COPY_POLY
#define COPY_POLY

#include <iostream>

#include <CGAL/circulator.h>
#include <CGAL/basic.h>

#include "Mesh.h"
#include "Parser.h"
#include "Polyhedron_Builder.h"

template <class HDS,class Polyhedron,class kernel>
class CModifierCopyPoly : public CGAL::Modifier_base<HDS>
{
private:
	typedef typename HDS::Vertex Vertex;
	typedef typename Vertex::Point Point;
	typedef typename HDS::Face_handle Face_handle;
	typedef typename HDS::Halfedge_handle Halfedge_handle;
	typedef typename CGAL::Enriched_polyhedron_incremental_builder_3<HDS> builder;

	typedef typename kernel::FT FT;
	typedef typename kernel::Vector_3 Vector;  //ajout Céline
	typedef typename kernel::Plane_3 Plane;
	typedef typename CGAL::Triangle_3<kernel> Triangle_3;
	typedef typename Polyhedron::Vertex_handle Vertex_handle;
	typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
	typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
	typedef typename Polyhedron::Edge_iterator Edge_iterator;
	typedef typename Polyhedron::Facet_iterator Facet_iterator;

	typedef typename Polyhedron::Halfedge_around_vertex_circulator
											Halfedge_around_vertex_circulator;  // fin ajout Céline
	typedef typename Polyhedron::Halfedge_around_facet_circulator
											Halfedge_around_facet_circulator;  // fin ajout Céline

	//private fileds :
	Polyhedron *m_pMesh;

public:

	// life cycle
	CModifierCopyPoly(Polyhedron *pMesh)
	{
		CGAL_assertion(pMesh != NULL);
		m_pMesh = pMesh;
	}

	~CModifierCopyPoly() {}

 //////////////////////////////////////////////// BUILDER /////////////////////////////////////////////

	void operator()( HDS& hds)
	{
		builder B(hds,true);        //builder init
		B.begin_surface(3,1,6);     //initialisation of the new polyhedron
			add_vertices(B);            //create all vertices for the new polyhedron
			add_facets(B);
		B.end_surface();
	}

	// add vertices
	void add_vertices(builder &B)
	{
		int index = 0;
		Vertex_iterator pVertex;
		for(pVertex = m_pMesh->vertices_begin() ; pVertex != m_pMesh->vertices_end() ; pVertex++)
		{
			pVertex->tag(index);       //tag each original vertex
            B.add_vertex(pVertex->point());  //add original vertices to the new poly
			index++;
		}
	}


	void add_facets(builder &B)
	{
		Facet_iterator pFacet;
		for(pFacet = m_pMesh->facets_begin() ; pFacet != m_pMesh->facets_end() ; pFacet++)
		{
				unsigned int degree = Polyhedron::degree(pFacet);  //facet degree
				CGAL_assertion(degree >= 3);

                Halfedge_handle pHalfedge = pFacet->halfedge();
//                int tag = 0;//MT

                B.begin_facet();
                do
                {
                    B.add_vertex_to_facet(pHalfedge->vertex()->tag());
                    pHalfedge = pHalfedge->next();
                } while (pHalfedge != pFacet->halfedge());

                B.end_facet();
		}
		CGAL_assertion(!B.check_unconnected_vertices());
	}
};

template <class Polyhedron, class kernel>
class CCopyPoly
{
   typedef typename kernel::Vector_3 Vector;
   typedef typename Polyhedron::HalfedgeDS HalfedgeDS;

public:

	CCopyPoly() {}
	~CCopyPoly() {}

public:
   int copy(Polyhedron *OriginalMesh, Polyhedron *NewMesh)
   {
		CModifierCopyPoly<HalfedgeDS, Polyhedron, kernel> builder(OriginalMesh);
        NewMesh->delegate(builder);                     //calls the `operator()' of the `modifier'
        return(0);
   }
};


/*class CopyMesh
{
protected:
	MeshBuilder* m_build;

public:
	bool copy(Mesh *originalMesh, Mesh *newMesh)
	{
		RepeatingMeshBuilder b(newMesh);
		
		m_build->startSurface(1, 1, 3);
		m_build->setVtxStart(0);

		Mesh::Vertex_const_iterator vit = originalMesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = originalMesh->vertices_end();
		for (; vit != vit_end; vit++)
			m_build->addVtx(vit->point().x(), vit->point().y(), vit->point().z());

		FaceList faces;

		Mesh::Facet_const_iterator fit = originalMesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = originalMesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			Face f;

			Mesh::Halfedge_around_facet_const_circulator	fvait = fit->halfedge()->facet_begin();
			Mesh::Halfedge_around_facet_const_circulator	begin = fvait;
			CGAL_For_all(fvait, begin)
				f.add(fvait->vertex()->index());
			
			faces.add(f);
		}

		m_build->processFaces(faces);

		return true;
	}
};*/

#endif // COPY_POLY