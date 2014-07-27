#ifndef POINTSURFACEGENERATOR_H
#define POINTSURFACEGENERATOR_H
#define TETLIBRARY 1
#include <queue>
#include "tetgen.h"
#include "vol_math_Base.h"
#include "vol_math_Mesh.h"

class PointSurfaceGenerator
{
	struct SingleTetraFace
	{
		int initializedCount;
		int FIndex[4];
		SingleTetraFace()
		{
			FIndex[0]=-1;
			FIndex[1]=-1;
			FIndex[2]=-1;
			FIndex[3]=-1;
			initializedCount=0;
		}
	};
private:
	static double GetLength(Point3d& p1,Point3d& p2)
	{
		return sqrt((p1.X-p2.X)*(p1.X-p2.X)+(p1.Y-p2.Y)*(p1.Y-p2.Y)+(p1.Z-p2.Z)*(p1.Z-p2.Z));
	}
private:
    std::vector<Point3d>& points;
    TetraMesh mesh;
	std::vector<Triangle> InnerTriangles;
    std::vector<std::vector<int> > FTAdjs;
    std::vector<SingleTetraFace> TetraToFace;
    std::vector<bool> BoundaryFlags;
	std::vector<bool> TetraFlags;
	std::vector<bool> TriangleFlags;
	std::queue<int> queue;
	std::vector<int> vertexMap;
	Mesh retMesh;
public:
    PointSurfaceGenerator(std::vector<Point3d>& pointList):points(pointList){}
    ~PointSurfaceGenerator(){}
private:
    void Tetrahedralize()
	{
		tetgenbehavior b;
		tetgenio in,out;
		char pstr[]="nn";
		b.parse_commandline(pstr);
		in.numberofpoints=points.size();
		in.pointlist=new REAL[points.size()*3];
		for(size_t i=0;i<points.size();i++)
		{
			in.pointlist[3*i]=points[i].X;
			in.pointlist[3*i+1]=points[i].Y;
			in.pointlist[3*i+2]=points[i].Z;
		}
		tetrahedralize(&b, &in, &out, 0, 0);
		for(int i=0;i<out.numberofpoints;i++)
		{
			Point3d p((double)(out.pointlist[3*i]),(double)(out.pointlist[3*i+1]),(double)(out.pointlist[3*i+2]));
			mesh.AddVertex(p);
		}
		for(int i=0;i<out.numberoftrifaces;i++)
		{
			Triangle t(out.trifacelist[3*i],out.trifacelist[3*i+1],out.trifacelist[3*i+2]);
			mesh.AddFace(t);
		}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			Tetrahedra te(out.tetrahedronlist[4*i],out.tetrahedronlist[4*i+1],out.tetrahedronlist[4*i+2],out.tetrahedronlist[4*i+3]);
			mesh.AddTetra(te);
		}
	}
	void InitBoundaryFlags()
	{
		if(InnerTriangles.size()==0||FTAdjs.size()==0)
			throw std::exception("inner trianlge zero!");
		this->BoundaryFlags.resize(InnerTriangles.size(),false);
		for(size_t i=0;i<InnerTriangles.size();i++)
		{
			size_t neigborCount=FTAdjs[i].size();
			if(neigborCount==1)
				BoundaryFlags[i]=true;
			else if(neigborCount==2)
				BoundaryFlags[i]=false;
			else
				throw std::exception("odd FT adj info!");
		}
	}
	void InitTetraToFaceRelationship()
	{
		if(FTAdjs.size()==0)
			throw std::exception("ft adj=0!");
		TetraToFace.resize(mesh.Tedtrahedras.size());
		std::vector<bool> visited(mesh.Tedtrahedras.size(),false);
		for(size_t i=0;i<InnerTriangles.size();i++)
		{
			std::vector<int>& adjts=FTAdjs[i];
			for(size_t j=0;j<adjts.size();j++)
			{
				SingleTetraFace& stf=TetraToFace[adjts[j]];
				stf.FIndex[stf.initializedCount]=i;
				stf.initializedCount++;
			}

		}
	}
	bool HasValidTetra(bool findex,std::vector<bool>& tetraFlags)
	{
		std::vector<int>& adjts=FTAdjs[findex];
		for(size_t i=0;i<adjts.size();i++)
		{
			if(tetraFlags[adjts[i]])
				return true;
		}
		return false;
	}
	void DeleteOldAndSaveNew(int findex,std::vector<int>& newtempfindex)
	{
		newtempfindex.clear();
		TriangleFlags[findex]=false;
		std::vector<int>& adjts=FTAdjs[findex];
		for(size_t j=0;j<adjts.size();j++)
		{
			TetraFlags[adjts[j]]=false;
		}
		for(size_t j=0;j<adjts.size();j++)
		{
			for(int k=0;k<4;k++)
			{
				int f2index=TetraToFace[adjts[j]].FIndex[k];
				if(HasValidTetra(f2index,TetraFlags)&&TriangleFlags[f2index])
				{
					newtempfindex.push_back(f2index);
				}
				else
				{
					TriangleFlags[f2index]=false;
				}
			}
		}
	}
	Mesh* GenerateNewMesh()
	{
		retMesh.Vertices.clear();
		retMesh.Faces.clear();
		Mesh* m=new Mesh();
		vertexMap.resize(mesh.Vertices.size(),-1);
		for(size_t i=0;i<InnerTriangles.size();i++)
		{
			if(BoundaryFlags[i])
			{
				Triangle t=InnerTriangles[i];
				if(vertexMap[t.P0Index]==-1)
				{
					Point3d p=mesh.Vertices[t.P0Index];
					vertexMap[t.P0Index]=retMesh.AddVertex(p);
				}
				if(vertexMap[t.P1Index]==-1)
				{
					Point3d p=mesh.Vertices[t.P1Index];
					vertexMap[t.P1Index]=retMesh.AddVertex(p);
				}
				if(vertexMap[t.P2Index]==-1)
				{
					Point3d p=mesh.Vertices[t.P2Index];
					vertexMap[t.P2Index]=retMesh.AddVertex(p);
				}
				t.P0Index=vertexMap[t.P0Index];
				t.P1Index=vertexMap[t.P1Index];
				t.P2Index=vertexMap[t.P2Index];
				retMesh.AddFace(t);
			}
		}
		m->Vertices.swap(retMesh.Vertices);
		m->Faces.swap(retMesh.Faces);
		return m;
	}
public:
    Mesh* ExecuteSurfaceReconstruction(double maxlength)
    {
		if(InnerTriangles.size()==0)
		{
			Tetrahedralize();
			mesh.InitInnerTriangles();
			InnerTriangles.assign(mesh.InnerTriangles.begin(),mesh.InnerTriangles.end());
			FTAdjs.assign(mesh.NeighborTetraPerTriangle.begin(),mesh.NeighborTetraPerTriangle.end());
		}
		InitBoundaryFlags();
		InitTetraToFaceRelationship();
		TriangleFlags.resize(InnerTriangles.size(),true);
		TetraFlags.resize(mesh.Tedtrahedras.size(),true);
		for(size_t i=0;i<BoundaryFlags.size();i++)
		{
			if(BoundaryFlags[i])
			{
				queue.push(i);
				TriangleFlags[i]=false;
			}
		}
		std::vector<int> newtempfindex;
		newtempfindex.reserve(4);
		while(!queue.empty())
		{
			int findex=queue.front();
			queue.pop();
			if(!CheckTriangleValid(InnerTriangles[findex],maxlength))
			{
				DeleteOldAndSaveNew(findex,newtempfindex);
				BoundaryFlags[findex]=false;
				for(size_t i=0;i<newtempfindex.size();i++)
				{
					queue.push(newtempfindex[i]);
					BoundaryFlags[newtempfindex[i]]=true;
				}
			}
		}
		return GenerateNewMesh();
    }
protected:
	bool CheckTriangleValid(Triangle& tri,double maxlength)
	{
		Point3d& p0=mesh.Vertices[tri.P0Index];
		Point3d& p1=mesh.Vertices[tri.P1Index];
		Point3d& p2=mesh.Vertices[tri.P2Index];
		double p01=GetLength(p0,p1);
		if(p01>maxlength)
			return false;
		double p12=GetLength(p1,p2);
		if(p12>maxlength)
			return false;
		double p02=GetLength(p0,p2);
		if(p02>maxlength)
			return false;
		return true;
	}
};

#endif

