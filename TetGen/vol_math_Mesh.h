#ifndef MESH_H
#define MESH_H
#include "vol_math_Base.h"
#include <vector>
#include <algorithm>
class Mesh
{
public:
	std::vector<Point3d> Vertices;
	std::vector<Triangle> Faces;
	std::vector<Vector> VertexNormals;
public:
	Mesh()
	{
	}
	~Mesh()
	{
	}
	MESH_INDEXTYPE AddVertex(Point3d& toAdd)
	{
		MESH_INDEXTYPE index = (MESH_INDEXTYPE)Vertices.size();
		Vertices.push_back(toAdd);
		return index;
	}
	MESH_INDEXTYPE  AddFace(Triangle& triangle)
	{
		MESH_INDEXTYPE  index = (MESH_INDEXTYPE)Faces.size();
		Faces.push_back(triangle);
		return index;
	}
};
class TetraMesh : public Mesh
{
public:
	std::vector<Triangle> InnerTriangles;
	std::vector<Tetrahedra> Tedtrahedras;
	std::vector<Edge> Edges;
	std::vector<std::vector<int> > NeighborTetraPerTriangle;
	std::vector<std::vector<int> > NeighborTetraPerVertex;
	std::vector<std::vector<int> > NeighborTrianglePerEdge;
	std::vector<std::vector<int> > NeighborVertexPerVertex;
public:
	MESH_INDEXTYPE AddEdge(Edge& e)
	{
		MESH_INDEXTYPE index = (MESH_INDEXTYPE)Edges.size();
		Edges.push_back(e);
		return index;
	}
	MESH_INDEXTYPE AddTetra(Tetrahedra& tet)
	{
		MESH_INDEXTYPE index = (MESH_INDEXTYPE)Tedtrahedras.size();
		Tedtrahedras.push_back(tet);
		return index;
	}
	MESH_INDEXTYPE AddInnerTriangle(Triangle& t)
	{
		MESH_INDEXTYPE index = (MESH_INDEXTYPE)InnerTriangles.size();
		InnerTriangles.push_back(t);
		return index;
	}
	void InitNeighborTetraPerVertex()
	{
		NeighborTetraPerVertex.resize(Vertices.size());
		for(size_t i=0;i<Tedtrahedras.size();i++)
		{
			NeighborTetraPerVertex[Tedtrahedras[i].PIndex[0]].push_back(i);
		}
	}
	void InitInnerTriangles()
	{
		if(Tedtrahedras.size()==0)
			return;
		std::vector<IntIndexTriple> templist;
		for(size_t i=0;i<Tedtrahedras.size();i++)
		{
			Tetrahedra& te=Tedtrahedras[i];
			IntIndexTriple t1(te.PIndex[0],te.PIndex[2],te.PIndex[1]);
			IntIndexTriple t2(te.PIndex[0],te.PIndex[1],te.PIndex[3]);
			IntIndexTriple t3(te.PIndex[0],te.PIndex[3],te.PIndex[2]);
			IntIndexTriple t4(te.PIndex[1],te.PIndex[2],te.PIndex[3]);
			t1.Sort();
			t2.Sort();
			t3.Sort();
			t4.Sort();
			t1.Index=i;
			t2.Index=i;
			t3.Index=i;
			t4.Index=i;
			templist.push_back(t1);
			templist.push_back(t2);
			templist.push_back(t3);
			templist.push_back(t4);
		}
		std::sort(templist.begin(),templist.end());
		int uniquesize=1;
		Triangle t0(templist[0].X,templist[0].Y,templist[0].Z);
		InnerTriangles.push_back(t0);
		std::vector<int> vec0;vec0.reserve(2);
		vec0.push_back(templist[0].Index);
		NeighborTetraPerTriangle.push_back(vec0);
		for(size_t i=1;i<templist.size();i++)
		{
			if(!(templist[i]==templist[i-1]))
			{
				uniquesize++;
				Triangle t(templist[i].X,templist[i].Y,templist[i].Z);
				InnerTriangles.push_back(t);
				std::vector<int> vec;vec.reserve(2);
				vec.push_back(templist[i].Index);
				NeighborTetraPerTriangle.push_back(vec);
			}
			else
			{
				NeighborTetraPerTriangle[uniquesize-1].push_back(templist[i].Index);
			}
		}
		size_t sum2=0;
		size_t sum1=0;
		for(size_t i=0;i<NeighborTetraPerTriangle.size();i++)
		{
			if(NeighborTetraPerTriangle[i].size()==2)
				sum2++;
			if(NeighborTetraPerTriangle[i].size()==1)
				sum1++;
		}
		if(sum1+sum2!=NeighborTetraPerTriangle.size())
			throw std::exception();
	}
	void InitSurfaceEdges()
	{
		std::vector<Triangle>& faces=this->Faces;
		if(faces.size()==0)
			return;
		std::vector<IntIndexDouble> templist;
		for(size_t i=0;i<faces.size();i++)
		{
			Triangle& te=faces[i];
			IntIndexDouble t1(te.P0Index,te.P1Index);
			IntIndexDouble t2(te.P0Index,te.P2Index);
			IntIndexDouble t3(te.P1Index,te.P2Index);
			t1.Sort();
			t2.Sort();
			t3.Sort();
			t1.Index=i;
			t2.Index=i;
			t3.Index=i;
			templist.push_back(t1);
			templist.push_back(t2);
			templist.push_back(t3);
		}
		std::sort(templist.begin(),templist.end());
		int uniquesize=1;
		Edge t0(templist[0].X,templist[0].Y);
		Edges.push_back(t0);
		std::vector<int> vec0;
		vec0.push_back(templist[0].Index);
		NeighborTrianglePerEdge.push_back(vec0);
		for(size_t i=1;i<templist.size();i++)
		{
			if(!(templist[i]==templist[i-1]))
			{
				uniquesize++;
				Edge t(templist[i].X,templist[i].Y);
				Edges.push_back(t);
				std::vector<int> vec;
				vec.push_back(templist[i].Index);
				NeighborTrianglePerEdge.push_back(vec);
			}
			else
			{
				NeighborTrianglePerEdge[uniquesize-1].push_back(templist[i].Index);
			}
		}
	}
	void InitNeighborVertexPerVertex()
	{

	}
};
#endif
