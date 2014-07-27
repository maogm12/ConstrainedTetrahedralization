#ifndef IOFILE_H
#define IOFILE_H
#include <vector>
#include <stdio.h>
#include <algorithm>
#include "vol_math_Base.h"
#include "vol_math_Mesh.h"
#include "vol_math_ConstraintTetralization.h"

class PlyManager
{
public:
	static void ReadXYZFile(std::vector<Point3d>& list, const char* fileName)
	{
		FILE * nfile = fopen(fileName,"r");
		Point3d p;
		while(fscanf(nfile,"%lf %lf %lf \n",&p.X,&p.Y,&p.Z)!=-1)
		{
			list.push_back(p);
		}
		fclose(nfile);
	}
	static void WriteXYZFile(std::vector<Int16Triple>& list, const char* fileName);
	static void WriteXYZFile(std::vector<Point3d>& list, const char* fileName)
	{
		FILE * nfile = fopen(fileName,"wb");
		fprintf(nfile,"%d\n",list.size());
		for(size_t i=0;i<list.size();i++)
		{
			fprintf(nfile,"%f %f %f\n",list[i].X,list[i].Y,list[i].Z);
		}
		fclose(nfile);
	}
	static void Output(Mesh& mesh,const char* filename);
	static void Output_C(Mesh& mesh,std::vector<Color> & colors,const char* filename);
	static void Output_C2(TetraMesh& mesh,std::vector<Color> & colors,const char* filename);
	static void Output2(std::vector<Point3d> &points,const char* fileName);
	static void ReadFile(Mesh& mesh,const char* fileName);
	static void ReadFileEx(Mesh& mesh,const char* fileName);
	static void Output3(std::vector<Point3d> &point,std::vector<Vector> &normals,const char* fileName);
	static void Output4(TetraMesh& mesh,const char* fileName,int type);
private:
	static void AWriteV(FILE* sw, double v1, double v2, double v3,unsigned char r,unsigned char g,unsigned char b);
	static void AWriteF(FILE* sw, MESH_INDEXTYPE i1, MESH_INDEXTYPE i2, MESH_INDEXTYPE i3);
	static void AWriteE(FILE* sw,MESH_INDEXTYPE i1,MESH_INDEXTYPE i2);
};
class DatManager
{
	struct LFloatTriple
	{
		double X;
		double Y;
		double Z;
		int Index;
		bool operator<(const LFloatTriple &right) const
		{
			if(X < right.X)
				return true;
			else if(X>right.X)
				return false;
			else if(Y < right.Y)
				return true;
			else if(Y >right.Y)
				return false;
			else if(Z<right.Z)
				return true;
			else if(Z>right.Z)
				return false;
			else
				return false;
		}
	};
	struct OringinalTriangle
	{
		Point3d P[3];
		OringinalTriangle(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2)
		{
			P[0].X=x0;
			P[0].Y=y0;
			P[0].Z=z0;
			P[1].X=x1;
			P[1].Y=y1;
			P[1].Z=z1;
			P[2].X=x2;
			P[2].Y=y2;
			P[2].Z=z2;
		}
	};
public:
	std::vector<std::vector<OringinalTriangle> > alltriangles;
	std::vector<OringinalTriangle> merges;
	Mesh mergeMesh;
	std::vector<int> idtomerge;
	std::vector<int> ids;
	std::vector<int> types;
	Box3LFloat box;
public:
	std::vector<Mesh> meshes;
	Mesh boxBound[6];
public:
    void ReadDatFile(const char* filename)
	{
		char line[4096];
		FILE* file=fopen(filename,"r");
		int facecount;
		fscanf(file,"%d\n",&facecount);
		alltriangles.resize(facecount);
		ids.resize(facecount,-1);
		types.resize(facecount,-1);
		fscanf(file,"%lf,%lf,%lf,%lf,%lf,%lf\n",&box.Min3[0],&box.Min3[1],&box.Min3[2],&box.Max3[0],&box.Max3[1],&box.Max3[2]);
		for(int k=0;k<facecount;k++)
		{
			std::vector<OringinalTriangle> otriangles;
			fgets(line,sizeof(line),file);
			int id=-1,type=-1;
			fscanf(file,"id %d, type %d\n",&id,&type);
			int tcount=0;
			fscanf(file,"%d\n",&tcount);
			ids[k]=id;types[k]=type;
			printf("type:%d,count:%d\n",type,tcount);
			otriangles.reserve(tcount);
			for(int i=0;i<tcount;i++)
			{
				double x0,y0,z0,x1,y1,z1,x2,y2,z2;
				fscanf(file,"TV %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &x0,&y0,&z0,&x1,&y1,&z1,&x2,&y2,&z2);
				OringinalTriangle ot(x0,y0,z0,x1,y1,z1,x2,y2,z2);
				otriangles.push_back(ot);
			}
			alltriangles[k].swap(otriangles);
		}
	}
	void ReadSBDatFile(const char* filename)
	{
		FILE* file=fopen(filename,"r");
		fscanf(file,"range:%lf,%lf,%lf,%lf,%lf,%lf\n",&box.Min3[0],&box.Min3[1],&box.Min3[2],&box.Max3[0],&box.Max3[1],&box.Max3[2]);
		int facecount;
		fscanf(file,"nface = %d\n",&facecount);
		alltriangles.resize(facecount);
		ids.resize(facecount,-1);
		types.resize(facecount,-1);
		for(int k=0;k<facecount;k++)
		{
			std::vector<OringinalTriangle> otriangles;
			int tcount=0;
			char name[30];
			fscanf(file,"name = %s = %d",name,&tcount);
			otriangles.reserve(tcount);
			for(int i=0;i<tcount;i++)
			{
				double x0,y0,z0,x1,y1,z1,x2,y2,z2;
				fscanf(file,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &x0,&y0,&z0,&x1,&y1,&z1,&x2,&y2,&z2);
				OringinalTriangle ot(x0,y0,z0,x1,y1,z1,x2,y2,z2);
				otriangles.push_back(ot);
			}
			alltriangles[k].swap(otriangles);
		}
	}
	void GetMesh()
	{
		meshes.resize(GetCount());
		for(int k=0;k<GetCount();k++)
		{
			std::vector<OringinalTriangle>& otriangles=alltriangles[k];
			Mesh *retMesh=&meshes[k];
			std::vector<LFloatTriple> vec;
			std::vector<Triangle> triangles;
			triangles.reserve(otriangles.size());
			vec.reserve(otriangles.size()*3);
			for(size_t i=0;i<otriangles.size();i++)
			{
				LFloatTriple ft0,ft1,ft2;
				ft0.X=otriangles[i].P[0].X;
				ft0.Y=otriangles[i].P[0].Y;
				ft0.Z=otriangles[i].P[0].Z;
				ft1.X=otriangles[i].P[1].X;
				ft1.Y=otriangles[i].P[1].Y;
				ft1.Z=otriangles[i].P[1].Z;
				ft2.X=otriangles[i].P[2].X;
				ft2.Y=otriangles[i].P[2].Y;
				ft2.Z=otriangles[i].P[2].Z;
				int count=vec.size();
				ft0.Index=count;
				ft1.Index=count+1;
				ft2.Index=count+2;
				vec.push_back(ft0);
				vec.push_back(ft1);
				vec.push_back(ft2);
				Triangle t(count,count+1,count+2);
				triangles.push_back(t);
			}
			std::vector<int> map1,map2;
			map1.resize(otriangles.size()*3);
			map2.resize(otriangles.size()*3);
			std::sort(vec.begin(),vec.end());
			for(size_t i=0;i<vec.size();i++)
				map1[vec[i].Index]=i;
			int index=0;
			std::vector<Point3d> points;
			Point3d p0(vec[0].X,vec[0].Y,vec[0].Z);
			points.push_back(p0);
			for(size_t i=1;i<vec.size();i++)
			{
				if(!(vec[i]<vec[i-1])&&!(vec[i-1]<vec[i]))
				{
					map2[i]=index;
				}
				else
				{
					index++;
					map2[i]=index;
					Point3d pi(vec[i].X,vec[i].Y,vec[i].Z);
					points.push_back(pi);
				}
			}
			for(size_t i=0;i<triangles.size();i++)
			{
				triangles[i].P0Index=map2[map1[triangles[i].P0Index]];
				triangles[i].P1Index=map2[map1[triangles[i].P1Index]];
				triangles[i].P2Index=map2[map1[triangles[i].P2Index]];
			}
			retMesh->Vertices.swap(points);
			retMesh->Faces.swap(triangles);
		}
	}
	void GetMergedMesh()
	{
		for(size_t i=0;i<alltriangles.size();i++)
		{
			//if(Contains(idtomerge,ids[i]))
			{
				merges.insert(merges.end(),alltriangles[i].begin(),alltriangles[i].end());
			}
		}
		std::vector<OringinalTriangle>& otriangles=merges;
		Mesh *retMesh=&mergeMesh;
		std::vector<LFloatTriple> vec;
		std::vector<Triangle> triangles;
		triangles.reserve(otriangles.size());
		vec.reserve(otriangles.size()*3);
		for(size_t i=0;i<otriangles.size();i++)
		{
			LFloatTriple ft0,ft1,ft2;
			ft0.X=otriangles[i].P[0].X;
			ft0.Y=otriangles[i].P[0].Y;
			ft0.Z=otriangles[i].P[0].Z;
			ft1.X=otriangles[i].P[1].X;
			ft1.Y=otriangles[i].P[1].Y;
			ft1.Z=otriangles[i].P[1].Z;
			ft2.X=otriangles[i].P[2].X;
			ft2.Y=otriangles[i].P[2].Y;
			ft2.Z=otriangles[i].P[2].Z;
			int count=vec.size();
			ft0.Index=count;
			ft1.Index=count+1;
			ft2.Index=count+2;
			vec.push_back(ft0);
			vec.push_back(ft1);
			vec.push_back(ft2);
			Triangle t(count,count+1,count+2);
			triangles.push_back(t);
		}
		std::vector<int> map1,map2;
		map1.resize(otriangles.size()*3);
		map2.resize(otriangles.size()*3);
		std::sort(vec.begin(),vec.end());
		for(size_t i=0;i<vec.size();i++)
			map1[vec[i].Index]=i;
		int index=0;
		std::vector<Point3d> points;
		Point3d p0(vec[0].X,vec[0].Y,vec[0].Z);
		points.push_back(p0);
		for(size_t i=1;i<vec.size();i++)
		{
			if(!(vec[i]<vec[i-1])&&!(vec[i-1]<vec[i]))
			{
				map2[i]=index;
			}
			else
			{
				index++;
				map2[i]=index;
				Point3d pi(vec[i].X,vec[i].Y,vec[i].Z);
				points.push_back(pi);
			}
		}
		for(size_t i=0;i<triangles.size();i++)
		{
			triangles[i].P0Index=map2[map1[triangles[i].P0Index]];
			triangles[i].P1Index=map2[map1[triangles[i].P1Index]];
			triangles[i].P2Index=map2[map1[triangles[i].P2Index]];
		}
		retMesh->Vertices.swap(points);
		retMesh->Faces.swap(triangles);
	}
	void GetBoxBound()
	{
//		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
//		Point3d p1(x0,y0,z0),p2(x1,y0,z0),p3(x1,y1,z0),p4(x0,y1,z0);
//		Point3d p5(x0,y0,z1),p6(x1,y0,z1),p7(x1,y1,z1),p8(x0,y1,z1);
//		{
//			Mesh& m=boxBound[0];
//			m.AddVertex(p1);m.AddVertex(p2);m.AddVertex(p3);m.AddVertex(p4);
//			m.AddFace(Triangle(1-1,4-1,3-1));
//			m.AddFace(Triangle(1-1,3-1,2-1));
//		}
//		{
//			Mesh& m=boxBound[1];
//			m.AddVertex(p5);m.AddVertex(p6);m.AddVertex(p7);m.AddVertex(p8);
//			m.AddFace(Triangle(1-1,4-1,3-1));
//			m.AddFace(Triangle(1-1,3-1,2-1));
//		}
//		{
//			Mesh& m=boxBound[2];
//			m.AddVertex(p1);m.AddVertex(p4);m.AddVertex(p8);m.AddVertex(p5);
//			m.AddFace(Triangle(1-1,4-1,3-1));
//			m.AddFace(Triangle(1-1,3-1,2-1));
//		}
//		{
//			Mesh& m=boxBound[3];
//			m.AddVertex(p2);m.AddVertex(p3);m.AddVertex(p7);m.AddVertex(p6);
//			m.AddFace(Triangle(1-1,4-1,3-1));
//			m.AddFace(Triangle(1-1,3-1,2-1));
//		}
//		{
//			Mesh& m=boxBound[4];
//			m.AddVertex(p1);m.AddVertex(p2);m.AddVertex(p6);m.AddVertex(p5);
//			m.AddFace(Triangle(1-1,4-1,3-1));
//			m.AddFace(Triangle(1-1,3-1,2-1));
//		}
//		{
//			Mesh& m=boxBound[5];
//			m.AddVertex(p3);m.AddVertex(p4);m.AddVertex(p8);m.AddVertex(p7);
//			m.AddFace(Triangle(1-1,4-1,3-1));
//			m.AddFace(Triangle(1-1,3-1,2-1));
//		}
	}
	Box3LFloat GetBox()
	{
		return box;
	}
	int GetCount()
	{
		return alltriangles.size();
	}
	bool Contains(std::vector<int> ids, int id)
	{
		return std::find(ids.begin(),ids.end(),id)!=ids.end();
	}
};

#endif
