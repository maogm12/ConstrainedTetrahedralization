#include <iostream>
#include "tetgen.h"
#include "vol_math_Base.h"
#include "vol_math_IOFile.h"
#include "vol_math_Mesh.h"
#include "vol_math_ConstraintTetralization.h"
#include "vol_math_PointSurfaceGenerator.h"
#include "VertexSpliter.h"

#if 0
//void show_progress_bak(int total,int step)
//{
//	printf("\r:%d%%",(int)(((double)step*100/total)));
//}
//void Test1()
//{
//	Mesh m;
//	Point3d p0(0,0,0);
//	Point3d p1(1,0,0);
//	Point3d p2(0,1,0);
//	Point3d p3(1,1,0);
//	Point3d p4(0,0,1);
//	Point3d p5(1,0,1);
//	Point3d p6(0,1,1);
//	Point3d p7(1,1,1);
//	m.AddVertex(p0);
//	m.AddVertex(p1);
//	m.AddVertex(p2);
//	m.AddVertex(p3);
//	m.AddVertex(p4);
//	m.AddVertex(p5);
//	m.AddVertex(p6);
//	m.AddVertex(p7);
//
//	Point3d pa(0.5f,0.7f,1.0f);
//	Point3d pb(0.3f,0.5f,1.0f);
//	Point3d pc(0.9f,0.7f,1.0f);
//
//	m.AddVertex(pa);
//	m.AddVertex(pb);
//	m.AddVertex(pc);
//	TetraMesh* me=ConstraintTetrahedralizer::Tetrahedralize(m.Vertices);
//	PlyManager::Output4(*me,"Test_Surface.ply",1);
//	delete me;
//}
//void InitSample2(Mesh& m,Box3Float& outbox)
//{
//	PlyManager::ReadFile(m,"F:\\workspacecsharp\\C++Proj\\TetGen\\Debug\\onebody2.ply");
//	Box3Float box=m.GetBox3();
//	double xlen=box.GetXLength();
//	double x0=box.Min3[0]-0.2f*xlen,x1=box.Max3[0]+0.2f*xlen;
//	double ylen=box.GetYLength();
//	double y0=box.Min3[1]-0.2f*ylen,y1=box.Max3[1]+0.2f*ylen;
//	double zlen=box.GetZLength();
//	double z0=box.Min3[2]-zlen,z1=box.Max3[2]+zlen;
//	outbox=Box3Float(x0,y0,z0,x1,y1,z1);
//}
//void Test2()
//{
//	TetraMesh m;
//	Mesh* mesh=Sample::GetSampleMesh1(3,0);
//	PlyManager::Output(*mesh,"testedge.ply");
//	m.Vertices.swap(mesh->Vertices);
//	m.Faces.swap(mesh->Faces);
//	delete mesh;
//	m.InitSurfaceEdges();
//	for(size_t i=0;i<m.Edges.size();i++)
//	{
//		Edge& e=m.Edges[i];
//		printf("edge:[%d,%d],neigbtri:",e.PIndex[0],e.PIndex[1]);
//		std::vector<int>& neib=m.NeighborTrianglePerEdge[i];
//		for(size_t j=0;j<neib.size();j++)
//		{
//			Triangle &t=m.Faces[neib[j]];
//			printf("[%d,%d,%d],",t.P0Index,t.P1Index,t.P2Index);
//		}
//		printf("\n");
//	}
//	system("pause");
//}
//void Test3()
//{
//	Mesh m;
//	PlyManager::ReadFileEx(m,"tetragridtest2\\mergeall.ply");
//	Sample::CleanMeshTVertices(m);
//	system("pause");
//}
//void Test4()
//{
//	Mesh m;
//	DatManager dat;
//	dat.ReadSBDatFile("tetragridtest3\\Yhs_test_1.dat");
//	dat.GetMesh();
//	for(int i=0;i<dat.GetCount();i++)
//	{
//		char name[100];
//		sprintf(name,"tetragridtest3\\%d.ply",i);
//		PlyManager::Output(dat.meshes[i],name);
//	}
//	dat.idtomerge.push_back(-1);
//	dat.GetMergedMesh();
//	PlyManager::Output(dat.mergeMesh,"tetragridtest3\\mergeall.ply");
//	dat.GetBoxBound();
//	for(int i=0;i<6;i++)
//	{
//		char name[100];
//		sprintf(name,"tetragridtest3\\face%d.ply",i);
//		PlyManager::Output(dat.boxBound[i],name);
//	}
//}
//void Test5(char* filename)
//{
//	Mesh m;
//	PlyManager::ReadFileEx(m,filename);
//	Box3Float box(0.000000,0.000000,-5000.000000,5000.000000,5000.0000000,0.000000);
//	//Box3Float box(0,0,0,49,49,49);
//	//Box3Float box(0,0,-25,49,49,25);
//	tetgenio in,out;
//	tetgenbehavior b;
//	ConstraintTetrahedralizer::InitTetgenIO(in,out,b);
//	ConstraintTetrahedralizer::InitPoints(in,m.Vertices,box);
//	ConstraintTetrahedralizer::InitFacets(in,m.Vertices,m.Faces,box);
//	//ConstraintTetrahedralizer::InitRegion(in,m.Vertices,m.Faces,box);
//	tetrahedralize(&b,&in,&out);
//	TetraMesh* outM= ConstraintTetrahedralizer::Out2TetraMesh(out);
//	printf("%d %d %d %d\n", outM->Vertices.size(),outM->Faces.size(),outM->InnerTriangles.size(),outM->Tedtrahedras.size());
//	PlyManager::Output4(*outM,"Test_Const_2.ply",1);
//}
//void Test6()
//{
//	Mesh m;
//	Box3Float box(0,0,0,1,1,1);
//	/*Point3d p9(0.1,0.5,0.5);
//	Point3d p10(0.5,0.45,0.5);
//	Point3d p11(0.9,0.5,0.5);
//	Point3d p12(0.5,0.65,0.5);
//	m.AddVertex(p9);
//	m.AddVertex(p10);
//	m.AddVertex(p11);
//	m.AddVertex(p12);
//	m.AddFace(Triangle(0,1,2));
//	m.AddFace(Triangle(0,2,3));*/
//
//	Point3d p0(0.5f,0.5f,0.5f);
//	Point3d p1(0.5f,0,0.75f);
//	Point3d p2(0.5f,0,0.5f);
//	Point3d p3(0.5f,0,0.25f);
//	Point3d p4(0.75f,0,0.5f);
//	m.AddVertex(p0);
//	m.AddVertex(p1);
//	m.AddVertex(p2);
//	m.AddVertex(p3);
//	m.AddVertex(p4);
//	m.AddFace(Triangle(0,1,2));
//	m.AddFace(Triangle(0,2,4));
//	m.AddFace(Triangle(0,2,3));
//
//	std::vector<Point3d> insertedPoints;
//	std::vector<Int16Triple> op;
//	//insertedPoints.push_back(Point3d(0.6,0.6,1.0));
//	ConstraintTetrahedralizer::BuildInsertedPoints(box,insertedPoints,0.2f,op);
//
//	tetgenio in,out,addin;
//	tetgenbehavior b;
//	ConstraintTetrahedralizer::InitTetgenIO(in,out,b);
//	ConstraintTetrahedralizer::InitPoints(in,m.Vertices,box);
//	ConstraintTetrahedralizer::InitFacets(in,m.Vertices,m.Faces,box);
//	ConstraintTetrahedralizer::InitAddin(addin,insertedPoints);
//	tetrahedralize(&b,&in,&out,&addin);
//	TetraMesh* outM= ConstraintTetrahedralizer::Out2TetraMesh(out);
//	printf("%d %d %d %d\n", outM->Vertices.size(),outM->Faces.size(),outM->InnerTriangles.size(),outM->Tedtrahedras.size());
//	PlyManager::Output4(*outM,"Test_Const_2.ply",1);
//	system("pause");
//}
//void TestInsertPoints()
//{
//	/*Box3Float box(0,0,0,141,222,134);
//	double ulen=5.7f;
//	std::vector<Point3d> points;
//	ConstraintTetrahedralizer::BuildInsertedPoints(box,points,ulen,);
//	PlyManager::WriteXYZFile(points,"ins.xyz");*/
//}
//void Test7()
//{
//	Mesh* m=Sample::GetSampleMesh1(100,0);
//	Box3Float box=m->GetBox3();
//	ConstrainedMeshSmoother cms(m);
//	cms.InitBoundaryMaker();
//	int id=m->Vertices.size()/3;
//	int id2=id+1;
//	cms.SetMarker(id,2);
//	cms.SetMarker(id2,2);
//	m->Vertices[id].X=box.Min3[0]+0.33f*box.GetXLength();
//	m->Vertices[id].Y=box.Min3[1]+0.33f*box.GetYLength();
//	m->Vertices[id].Z=box.Min3[2]+0.33f*box.GetZLength();
//	m->Vertices[id2].X=box.Min3[0]+0.66f*box.GetXLength();
//	m->Vertices[id2].Y=box.Min3[1]+0.66f*box.GetYLength();
//	m->Vertices[id2].Z	=box.Min3[2]+0.66f*box.GetZLength();
//	PlyManager::Output(*m,"beforesm.ply");
//	cms.ExecuteConstrainedSmoothing(2000);
//	PlyManager::Output(*m,"afttersm.ply");
//	delete m;
//}
//void Test8()
//{
//	Mesh m;
//	m.AddVertex(Point3d(0.3f,0.5f,0.5f));
//	m.AddVertex(Point3d(0.4f,0.4f,0.5f));
//	m.AddVertex(Point3d(0.5f,0.5f,0.5f));
//	m.AddVertex(Point3d(0.6f,0.4f,0.5f));
//	m.AddVertex(Point3d(0.7f,0.5f,0.5f));
//	m.AddVertex(Point3d(0.5f,0.7f,0.5f));
//
//	m.AddFace(Triangle(0,1,2));
//	m.AddFace(Triangle(1,3,2));
//	m.AddFace(Triangle(2,3,4));
//	m.AddFace(Triangle(0,4,5));
//	//m.AddFace(Triangle(0,4,2));
//	Box3Float box(0,0,0,1,1,1);
//	PlyManager::Output(m,"non-v-sample.ply");
//	TetraMesh* outM= ConstraintTetrahedralizer::Execute(m,box);
//	printf("%d %d %d %d\n", outM->Vertices.size(),outM->Faces.size(),outM->InnerTriangles.size(),outM->Tedtrahedras.size());
//	PlyManager::Output4(*outM,"Test_Const.ply",1);
//	system("pause");
//}
//void Test9()
//{
//	std::vector<Point3d> vec;
//	vec.push_back(Point3d(0,0,0));
//	vec.push_back(Point3d(0,1,0));
//	vec.push_back(Point3d(1,0,0));
//	vec.push_back(Point3d(0,0,1));
//	vec.push_back(Point3d(0.5,0.5,-1));
//	//vec.push_back(Point3d(0.25,0.25,0));
//	TetraMesh* outM=ConstraintTetrahedralizer::Tetrahedralize(vec);
//	for(size_t i=0;i<outM->InnerTriangles.size();i++)
//	{
//		Triangle& t=outM->InnerTriangles[i];
//		printf("Triange : [%d,%d,%d]",t.P0Index,t.P1Index,t.P2Index);
//		std::vector<int> &adj=outM->NeighborTetraPerTriangle[i];
//		for(size_t j=0;j<adj.size();j++)
//		{
//			Tetrahedra& tet=outM->Tedtrahedras[adj[j]];
//			printf("Adj[%d] : [%d,%d,%d,%d] ",j,tet.PIndex[0],tet.PIndex[1],tet.PIndex[2],tet.PIndex[3]);
//		}
//		printf("\n");
//	}
//	printf("%d %d %d %d\n", outM->Vertices.size(),outM->Faces.size(),outM->InnerTriangles.size(),outM->Tedtrahedras.size());
//	system("pause");
//	PlyManager::Output4(*outM,"Test_Const.ply",1);
//}
//int Test10()
//{
//	tetgenio in, out;
//	tetgenio::facet *f;
//	tetgenio::polygon *p;
//	int i;
//
//	// All indices start from 1.
//	in.firstnumber = 1;
//
//	in.numberofpoints = 8;
//	in.pointlist = new REAL[in.numberofpoints * 3];
//	in.pointlist[0]  = 0;  // node 1.
//	in.pointlist[1]  = 0;
//	in.pointlist[2]  = 0;
//	in.pointlist[3]  = 2;  // node 2.
//	in.pointlist[4]  = 0;
//	in.pointlist[5]  = 0;
//	in.pointlist[6]  = 2;  // node 3.
//	in.pointlist[7]  = 2;
//	in.pointlist[8]  = 0;
//	in.pointlist[9]  = 0;  // node 4.
//	in.pointlist[10] = 2;
//	in.pointlist[11] = 0;
//	// Set node 5, 6, 7, 8.
//	for (i = 4; i < 8; i++) {
//		in.pointlist[i * 3]     = in.pointlist[(i - 4) * 3];
//		in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
//		in.pointlist[i * 3 + 2] = 4;
//	}
//
//	in.numberoffacets = 6;
//	in.facetlist = new tetgenio::facet[in.numberoffacets];
//	in.facetmarkerlist = new int[in.numberoffacets];
//
//	// Facet 1. The leftmost facet.
//	f = &in.facetlist[0];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 1;
//	p->vertexlist[1] = 2;
//	p->vertexlist[2] = 3;
//	p->vertexlist[3] = 4;
//
//	// Facet 2. The rightmost facet.
//	f = &in.facetlist[1];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 5;
//	p->vertexlist[1] = 6;
//	p->vertexlist[2] = 7;
//	p->vertexlist[3] = 8;
//
//	// Facet 3. The bottom facet.
//	f = &in.facetlist[2];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 1;
//	p->vertexlist[1] = 5;
//	p->vertexlist[2] = 6;
//	p->vertexlist[3] = 2;
//
//	// Facet 4. The back facet.
//	f = &in.facetlist[3];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 2;
//	p->vertexlist[1] = 6;
//	p->vertexlist[2] = 7;
//	p->vertexlist[3] = 3;
//
//	// Facet 5. The top facet.
//	f = &in.facetlist[4];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 3;
//	p->vertexlist[1] = 7;
//	p->vertexlist[2] = 8;
//	p->vertexlist[3] = 4;
//
//	// Facet 6. The front facet.
//	f = &in.facetlist[5];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 4;
//	p->vertexlist[1] = 8;
//	p->vertexlist[2] = 5;
//	p->vertexlist[3] = 1;
//
//	// Set 'in.facetmarkerlist'
//
//	in.facetmarkerlist[0] = -1;
//	in.facetmarkerlist[1] = -2;
//	in.facetmarkerlist[2] = 0;
//	in.facetmarkerlist[3] = 0;
//	in.facetmarkerlist[4] = 0;
//	in.facetmarkerlist[5] = 0;
//
//	tetgenbehavior b;
//	//b.quality=1;
//	b.plc=1;
//	tetrahedralize(&b, &in, &out);
//	TetraMesh *outM=ConstraintTetrahedralizer::Out2TetraMesh(out);
//	printf("%d %d %d %d\n", outM->Vertices.size(),outM->Faces.size(),outM->InnerTriangles.size(),outM->Tedtrahedras.size());
//	PlyManager::Output4(*outM,"Test_Const.ply",1);
//	delete outM;
//	return 0;
//}
//void Test11()
//{
//	tetgenio in, out;
//	tetgenio::facet *f;
//	tetgenio::polygon *p;
//	// All indices start from 1.
//	in.firstnumber = 1;
//
//	in.numberofpoints = 16;
//	in.pointlist = new REAL[in.numberofpoints * 3];
//	in.pointlist[0]  = 0;  // node 1.
//	in.pointlist[1]  = 0;
//	in.pointlist[2]  = 0;
//	in.pointlist[3]  = 1;  // node 2.
//	in.pointlist[4]  = 0;
//	in.pointlist[5]  = 0;
//	in.pointlist[6]  = 1;  // node 3.
//	in.pointlist[7]  = 1;
//	in.pointlist[8]  = 0;
//	in.pointlist[9]  = 0;  // node 4.
//	in.pointlist[10] = 1;
//	in.pointlist[11] = 0;
//	// Set node 5, 6, 7, 8.
//	for (int i = 4; i < 8; i++) {
//		in.pointlist[i * 3]     = in.pointlist[(i - 4) * 3];
//		in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
//		in.pointlist[i * 3 + 2] = 1;
//	}
//	in.pointlist[24]=0.1;
//	in.pointlist[25]=0.5;
//	in.pointlist[26]=0.5;
//
//	in.pointlist[27]=0.5;
//	in.pointlist[28]=0.45;
//	in.pointlist[29]=0.5;
//
//	in.pointlist[30]=0.9;
//	in.pointlist[31]=0.5;
//	in.pointlist[32]=0.5;
//
//	in.pointlist[33]=0.5;
//	in.pointlist[34]=0.55;
//	in.pointlist[35]=0.5;
//
//
//	in.pointlist[36]=0.25;
//	in.pointlist[37]=0.00;
//	in.pointlist[38]=0.00;
//
//	in.pointlist[39]=0.00;
//	in.pointlist[40]=0.25;
//	in.pointlist[41]=0.00;
//
//	in.pointlist[42]=0.00;
//	in.pointlist[43]=0.00;
//	in.pointlist[44]=0.8;
//
//	in.pointlist[45]=0.5;
//	in.pointlist[46]=0.5;
//	in.pointlist[47]=0.9;
//
//
//	in.numberoffacets = 9;
//	in.facetlist = new tetgenio::facet[in.numberoffacets];
//	in.facetmarkerlist = new int[in.numberoffacets];
//	// Facet 1. The leftmost facet.
//	f = &in.facetlist[0];
//	f->numberofpolygons = 2;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 6;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 1;
//	p->vertexlist[1] = 13;
//	p->vertexlist[2] = 2;
//	p->vertexlist[3] = 3;
//	p->vertexlist[4] = 4;
//	p->vertexlist[5] = 14;
//	p = &f->polygonlist[1];
//	p->numberofvertices = 2;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 14;
//
//	// Facet 2. The rightmost facet.
//	f = &in.facetlist[1];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 5;
//	p->vertexlist[1] = 6;
//	p->vertexlist[2] = 7;
//	p->vertexlist[3] = 8;
//
//	// Facet 3. The bottom facet.
//	f = &in.facetlist[2];
//	f->numberofpolygons =2;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 6;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 1;
//	p->vertexlist[1] = 15;
//	p->vertexlist[2] = 5;
//	p->vertexlist[3] = 6;
//	p->vertexlist[4] = 2;
//	p->vertexlist[5] =13;
//	p = &f->polygonlist[1];
//	p->numberofvertices = 2;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 15;
//
//
//	// Facet 4. The back facet.
//	f = &in.facetlist[3];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 2;
//	p->vertexlist[1] = 6;
//	p->vertexlist[2] = 7;
//	p->vertexlist[3] = 3;
//
//	// Facet 5. The top facet.
//	f = &in.facetlist[4];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 4;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 3;
//	p->vertexlist[1] = 7;
//	p->vertexlist[2] = 8;
//	p->vertexlist[3] = 4;
//
//	// Facet 6. The front facet.
//	f = &in.facetlist[5];
//	f->numberofpolygons = 2;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 6;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 4;
//	p->vertexlist[1] = 8;
//	p->vertexlist[2] = 5;
//	p->vertexlist[3] = 15;
//	p->vertexlist[4] = 1;
//	p->vertexlist[5] = 14;
//	p = &f->polygonlist[1];
//	p->numberofvertices = 2;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 14;
//	p->vertexlist[1] = 15;
//
//	//
//	f = &in.facetlist[6];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 9;
//	p->vertexlist[1] = 10;
//	p->vertexlist[2] = 11;
//
//	f = &in.facetlist[7];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 9;
//	p->vertexlist[1] = 11;
//	p->vertexlist[2] = 12;
//	/*p = &f->polygonlist[2];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 9;
//	p->vertexlist[1] = 11;
//	p->vertexlist[2] = 16;*/
//	/*p = &f->polygonlist[3];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 14;
//	p->vertexlist[2] = 15;
//	p = &f->polygonlist[4];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 9;
//	p->vertexlist[2] = 15;
//	p = &f->polygonlist[5];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 9;
//	p->vertexlist[2] = 10;*/
//
//	f = &in.facetlist[8];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 14;
//	p->vertexlist[2] = 15;
//
//	/*f = &in.facetlist[8];
//	f->numberofpolygons = 1;
//	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//	f->numberofholes = 0;
//	f->holelist = NULL;
//	p = &f->polygonlist[0];
//	p->numberofvertices = 3;
//	p->vertexlist = new int[p->numberofvertices];
//	p->vertexlist[0] = 13;
//	p->vertexlist[1] = 9;
//	p->vertexlist[2] = 15;*/
//
//	// Set 'in.facetmarkerlist'
//
//	in.facetmarkerlist[0] = -1;
//	in.facetmarkerlist[1] = -2;
//	in.facetmarkerlist[2] = 1;
//	in.facetmarkerlist[3] = 2;
//	in.facetmarkerlist[4] = 3;
//	in.facetmarkerlist[5] = 4;
//
//	/*in.numberofregions=2;
//	in.regionlist=new double[5*in.numberofregions];
//	in.regionlist[0]=0.9;
//	in.regionlist[1]=0.9;
//	in.regionlist[2]=0.9;
//	in.regionlist[3]=3;
//	in.regionlist[4]=1;
//
//	in.regionlist[5]=0.1;
//	in.regionlist[6]=0.1;
//	in.regionlist[7]=0.1;
//	in.regionlist[8]=2;
//	in.regionlist[9]=0.0001;*/
//
//
//	tetgenbehavior b;
//	//b.quality=1;
//	//b.regionattrib=1;
//	//b.varvolume=1;
//	b.plc=1;
//	b.nobisect=1;
//	tetrahedralize(&b, &in, &out);
//	//TetraMesh *outM=ConstraintTetrahedralizer::Out2TetraMesh(out,true);
//	//printf("%d %d %d %d\n", outM->Vertices.size(),outM->Faces.size(),outM->InnerTriangles.size(),outM->Tedtrahedras.size());
//	//PlyManager::Output4(*outM,"Test_Const.ply",1);
//	//delete outM;
//};
//void Test12()
//{
//	DatManager dat;
//	dat.ReadDatFile("tetragridtest2\\test_automodel-7df.dat");
//	printf("%d",dat.GetCount());
//	dat.GetMesh();
//	for(int i=0;i<dat.GetCount();i++)
//	{
//		char name[100];
//		sprintf(name,"tetragridtest2\\id-%d_type-%d.ply",dat.ids[i],dat.types[i]);
//		//PlyManager::Output(dat.meshes[i],name);
//	}
//	dat.idtomerge.push_back(130);
//	dat.idtomerge.push_back(137);
//	dat.idtomerge.push_back(144);
//	dat.idtomerge.push_back(6635);
//	dat.idtomerge.push_back(15877);
//	dat.idtomerge.push_back(22847);
//	dat.idtomerge.push_back(52007);
//	dat.idtomerge.push_back(52019);
//	dat.idtomerge.push_back(52031);
//	dat.GetMergedMesh();
//	PlyManager::Output(dat.mergeMesh,"mergeall.ply");
//	dat.GetBoxBound();
//	for(int i=0;i<6;i++)
//	{
//		char name[100];
//		sprintf(name,"tetragridtest2\\face%d.ply",i);
//		//PlyManager::Output(dat.boxBound[i],name);
//	}
//}
//void TestSp()
//{
//	Mesh m;
//	PlyManager::ReadFileEx(m,"tetragridtest2\\sample4.ply");
//	/*m.AddVertex(Point3d(0,0,0));
//	m.AddVertex(Point3d(1,0,0));
//	m.AddVertex(Point3d(1,1,0));
//	m.AddVertex(Point3d(0,1,0));
//	m.AddVertex(Point3d(0.5,0,0));
//
//	m.AddVertex(Point3d(0.5,1,0));
//	m.AddVertex(Point3d(1,0.5,0));
//	m.AddVertex(Point3d(0,0.5,0));
//	m.AddVertex(Point3d(0.5,0.5,-0.5));
//	m.AddVertex(Point3d(0.5,0,-0.5));
//
//	m.AddVertex(Point3d(0.5,1,-0.5));
//	m.AddVertex(Point3d(0.5,0.5,0));
//	m.AddVertex(Point3d(0.5,1,0.5));
//	m.AddVertex(Point3d(0.5,0.5,0.5));
//	m.AddVertex(Point3d(0.5,0,0.5));
//
//	m.AddFace(Triangle(0,4,11));
//	m.AddFace(Triangle(0,11,7));
//
//	m.AddFace(Triangle(7,11,5));
//	m.AddFace(Triangle(7,5,3));
//
//	m.AddFace(Triangle(11,6,2));
//	m.AddFace(Triangle(11,2,5));
//
//	m.AddFace(Triangle(11,4,1));
//	m.AddFace(Triangle(11,1,6));
//
//	m.AddFace(Triangle(11,13,14));
//	m.AddFace(Triangle(11,14,4));
//
//	m.AddFace(Triangle(11,5,12));
//	m.AddFace(Triangle(11,12,13));
//
//	m.AddFace(Triangle(11,4,9));
//	m.AddFace(Triangle(11,9,8));
//
//	m.AddFace(Triangle(11,8,10));
//	m.AddFace(Triangle(11,10,5));
//	PlyManager::Output(m,"mytest.ply");*/
//	Box3Float box;
//	for(size_t i=0;i<m.Vertices.size();i++)
//	{
//		box.UpdateRange(m.Vertices[i].X,m.Vertices[i].Y,m.Vertices[i].Z);
//	}
//	TetraMesh* tm=ConstraintTetrahedralizer::Execute(m,box,false,false,5);
//	PlyManager::Output4(*tm,"tetragridtest5_0.ply",1);
//	/*std::vector<Color> colors;colors.resize(tm->Vertices.size(),Color(255,255,255));
//	for(int i=0;i<m.Vertices.size();i++)
//	{
//	colors[i].R=255;
//	colors[i].G=0;
//	colors[i].B=0;
//	}*/
//	VertexSpliter3d vs3(*tm,m);
//	vs3.ProcessSpliting();
//
//	PlyManager::Output4(*tm,"tetragridtest5_1.ply",1);
//	//delete tm;
//	system("pause");
//}
//int plc;                                                         // '-p', 0.
//int psc;                                                         // '-s', 0.
//int refine;                                                      // '-r', 0.
//int quality;                                                     // '-q', 0.
//int nobisect;                                                    // '-Y', 0.
//int coarsen;                                                     // '-R', 0.
//int weighted;                                                    // '-w', 0.
//int brio_hilbert;                                                // '-b', 1.
//int incrflip;                                                    // '-l', 0.
//int flipinsert;                                                  // '-L', 0.
//int metric;                                                      // '-m', 0.
//int varvolume;                                                   // '-a', 0.
//int fixedvolume;                                                 // '-a', 0.
//int regionattrib;                                                // '-A', 0.
//int conforming;                                                  // '-D', 0.
//int insertaddpoints;                                             // '-i', 0.
//int diagnose;                                                    // '-d', 0.
//int convex;                                                      // '-c', 0.
//int nomergefacet;                                                // '-M', 0.
//int nomergevertex;                                               // '-M', 0.
//int noexact;                                                     // '-X', 0.
//int nostaticfilter;                                              // '-X', 0.
//int zeroindex;                                                   // '-z', 0.
//int facesout;                                                    // '-f', 0.
//int edgesout;                                                    // '-e', 0.
//int neighout;                                                    // '-n', 0.
//int voroout;                                                     // '-v', 0.
//int meditview;                                                   // '-g', 0.
//int vtkview;                                                     // '-k', 0.
//int nobound;                                                     // '-B', 0.
//int nonodewritten;                                               // '-N', 0.
//int noelewritten;                                                // '-E', 0.
//int nofacewritten;                                               // '-F', 0.
//int noiterationnum;                                              // '-I', 0.
//int nojettison;                                                  // '-J', 0.
//int reversetetori;                                               // '-R', 0.
//int docheck;                                                     // '-C', 0.
//int quiet;                                                       // '-Q', 0.
//int verbose;                                                     // '-V', 0.

//-p	Tetrahedralizes a piecewise linear complex (PLC).
//	-Y	Preserves the input surface mesh (does not modify it).
//	-r	Reconstructs a previously generated mesh.
//	-q	Refines mesh (to improve mesh quality).
//	-R	Mesh coarsening (to reduce the mesh elements).
//	-A	Assigns attributes to tetrahedra in different regions.
//	-a	Applies a maximum tetrahedron volume constraint.
//	-m	Applies a mesh sizing function.
//	-i	Inserts a list of additional points.
//	-O	Specifies the level of mesh optimization.
//	-S	Specifies maximum number of added points.
//	-T	Sets a tolerance for coplanar test (default 10−8).
//	-X	Suppresses use of exact arithmetic.
//	-M	No merge of coplanar facets or very close vertices.
//	-w	Generates weighted Delaunay (regular) triangulation.
//	-c	Retains the convex hull of the PLC.
//	-d	Detects self-intersections of facets of the PLC.
//	-z	Numbers all output items starting from zero.
//	-f	Outputs all faces to .face file.
//	-e	Outputs all edges to .edge file.
//	-n	Outputs tetrahedra neighbors to .neigh file.
//	-v	Outputs Voronoi diagram to files.
//	-g	Outputs mesh to .mesh file for viewing by Medit.
//	-k	Outputs mesh to .vtk file for viewing by Paraview.
//	-J	No jettison of unused vertices from output .node file.
//	-B	Suppresses output of boundary information.
//	-N	Suppresses output of .node file.
//	-E	Suppresses output of .ele file.
//	-F	Suppresses output of .face and .edge file.
//	-I	Suppresses mesh iteration numbers.
//	-C	Checks the consistency of the final mesh.
//	-Q	Quiet: No terminal output except errors.
//	-V	Verbose: Detailed information, more terminal output.
//	-h	Help: A brief instruction for using TetGen.
//3.2.6
//	-i Inserts additional points
//	The -i switch indicates to insert a list of additional points into a CDT (when
//	-p switch is used) or a previously generated mesh (when -r switch is used).
//	The list of additional file is read from file xxx-a.node, which xxx stands
//	for the input file name (i.e., xxx.poly or xxx.smesh, or xxx.ele, ...). This
//		switch is useful for refining a finite element or finite volume mesh using a list
//		of user-defined points. Following are some pointers that you may need be
//careful:
//• Points lie out of the mesh domain are ignored by TetGen.
//	• The mesh may not be constrained Delaunay or conforming Delaunay
//	any more after the insertion of additional points. However, in combi-
//	nation with the -q switch, TetGen will automatically add additional
//	points to ensure the conforming Delaunay property of the mesh.
#endif
void STest()
{
	DatManager dat;
	dat.ReadSBDatFile("tetragridtest4\\Yhs_test_2.dat");
	//dat.GetMesh();
	dat.GetMergedMesh();
	PlyManager::Output(dat.mergeMesh,"tetragridtest4\\mergeall.ply");

	return;
}
void Test13()
{
	Mesh m;
	PlyManager::ReadFileEx(m,"tetragridtest4\\mergeall.ply");
	Box3LFloat box(0.000000,0.000000,0.000000,7161.497070,10037.639648,6000.000000);
	//Box3Float box;
	/*for(size_t i=0;i<m.Vertices.size();i++)
	{
	box.UpdateRange(m.Vertices[i].X,m.Vertices[i].Y,m.Vertices[i].Z);
	}*/
	TetraMesh* tm=TetGen::Execute(m,box,true,true,10);
	VertexSpliter3d vs3(*tm,m);
	vs3.ProcessSpliting();
	printf("check point7 \n");
	std::vector<Color> colors;
	TetraRegionMarker trm(*tm);
	trm.MarkRegions(colors);
	for(size_t i=0;i<m.Vertices.size();i++)
	{
		colors[i].R=255;
		colors[i].G=0;
		colors[i].B=0;
	}
	printf("check point8 \n");
	PlyManager::Output_C2(*tm,colors,"tetragridtest5_5.ply");
	system("pause");
}
 void TestRC()
 {
	 std::vector<Point3d> vec;
	 /*vec.resize(5);
	 vec[0]=Point3d(0,0,0);
	 vec[1]=Point3d(1,0,0);
	 vec[2]=Point3d(0,1,0);
	 vec[3]=Point3d(0,0,0.5);
	 vec[4]=Point3d(0.45,0.45,0);*/
	 PlyManager::ReadXYZFile(vec,"test_data_2_result_ori.xyz");
	 PointSurfaceGenerator psg(vec); 
	 Mesh* m=psg.ExecuteSurfaceReconstruction(10);
	 PlyManager::Output(*m,"hehe.ply");
 }

void testFilterInsertedPoints()
{

}

int main()
{
	//Test13();
	TestRC();
	//PointSurfaceGenerator psg;

	return 0;
}
