#ifndef BOXSEGMENTTRIANGULATOR_H
#define BOXSEGMENTTRIANGULATOR_H
#include "vol_math_Mesh.h"
#include <queue>
//class Sample
//{
//	struct IntDouble
//	{
//		int X;
//		int Y;
//		int Index;
//		bool operator<(const IntDouble &right)
//		{
//			if(X < right.X)
//				return true;
//			else if(X>right.X)
//				return false;
//			else if(Y < right.Y)
//				return true;
//			else if(Y >right.Y)
//				return false;
//		}
//	};
//	struct FloatTriple
//	{
//		double X;
//		double Y;
//		double Z;
//		int Index;
//		bool operator<(const FloatTriple &right)
//		{
//			if(X < right.X)
//				return true;
//			else if(X>right.X)
//				return false;
//			else if(Y < right.Y)
//				return true;
//			else if(Y >right.Y)
//				return false;
//			else if(Z<right.Z)
//				return true;
//			else if(Z>right.Z)
//				return false;
//			else
//				return false;
//		}
//	};
//	struct OringinalTriangle
//	{
//		Point3d P[3];
//		OringinalTriangle(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2)
//		{
//			P[0].X=x0;
//			P[0].Y=y0;
//			P[0].Z=z0;
//			P[1].X=x1;
//			P[1].Y=y1;
//			P[1].Z=z1;
//			P[2].X=x2;
//			P[2].Y=y2;
//			P[2].Z=z2;
//		}
//	};
//	struct TRegion
//	{
//		std::vector<Point3d>* points;
//		std::vector<int> pindices;
//		int Id;
//		Box3LFloat regionBox;
//		bool isTHole;
//		int EndIndex[2];
//		void InitEnd()
//		{
//
//		}
//		void Sort()
//		{
//
//		}
//		Point3d GetCenter()
//		{
//			Point3d av(0,0,0);
//			for(size_t i=0;i<pindices.size();i++)
//			{
//				av.X+=(*points)[pindices[i]].X;
//				av.Y+=(*points)[pindices[i]].Y;
//				av.Z+=(*points)[pindices[i]].Z;
//			}
//			av.X/=pindices.size();
//			av.Y/=pindices.size();
//			av.Z/=pindices.size();
//			return av;
//		}
//	};
//
//public:
//	static Mesh* GetSampleMesh1(int width,double z)
//	{
//		Mesh* m = new Mesh();
//		int** hash=new int*[width];
//		for(int i=0;i<width;i++)
//			hash[i]=new int[width];
//		for(int i=0;i<width;i++)
//		{
//			for (int j=0;j<width;j++)
//			{
//				Point3d p((double)i,(double)j,(double)z);
//				hash[i][j]=m->AddVertex(p);
//			}
//		}
//		for(int i=0;i<width-1;i++)
//		{
//			for (int j = 0; j < width-1; j++)
//			{
//				Triangle t1(hash[i][j],hash[i+1][j],hash[i][j+1]);
//				Triangle t2(hash[i][j+1], hash[i + 1][ j], hash[i+1][ j + 1]);
//				m->AddFace(t1);
//				m->AddFace(t2);
//			}
//		}
//		for(int i=0;i<width;i++)
//			delete[] hash[i];
//		delete[] hash;
//	    return m;
//	}
//	static Mesh* GetSampleMesh1_v(int width,double y)
//	{
//		Mesh* m = new Mesh();
//		int** hash=new int*[width];
//		for(int i=0;i<width;i++)
//			hash[i]=new int[width];
//		for(int i=0;i<width;i++)
//		{
//			for (int j=0;j<width;j++)
//			{
//				Point3d p((double)i,(double)y,(double)j);
//				hash[i][j]=m->AddVertex(p);
//			}
//		}
//		for(int i=0;i<width-1;i++)
//		{
//			for (int j = 0; j < width-1; j++)
//			{
//				Triangle t1(hash[i][j],hash[i+1][j],hash[i][j+1]);
//				Triangle t2(hash[i][j+1], hash[i + 1][ j], hash[i+1][ j + 1]);
//				m->AddFace(t1);
//				m->AddFace(t2);
//			}
//		}
//		for(int i=0;i<width;i++)
//			delete[] hash[i];
//		delete[] hash;
//		return m;
//	}
//	static Mesh* GetSampleMesh1_v2(int width,double x)
//	{
//		Mesh* m = new Mesh();
//		int** hash=new int*[width];
//		for(int i=0;i<width;i++)
//			hash[i]=new int[width];
//		for(int i=0;i<width;i++)
//		{
//			for (int j=0;j<width;j++)
//			{
//				Point3d p((double)x,(double)i,(double)j);
//				hash[i][j]=m->AddVertex(p);
//			}
//		}
//		for(int i=0;i<width-1;i++)
//		{
//			for (int j = 0; j < width-1; j++)
//			{
//				Triangle t1(hash[i][j],hash[i+1][j],hash[i][j+1]);
//				Triangle t2(hash[i][j+1], hash[i + 1][ j], hash[i+1][ j + 1]);
//				m->AddFace(t1);
//				m->AddFace(t2);
//			}
//		}
//		for(int i=0;i<width;i++)
//			delete[] hash[i];
//		delete[] hash;
//		return m;
//	}
//	static Mesh* GetSampleMesh2()
//	{
//		Mesh *m=new Mesh();
//		Point3d p0(0,0,1);
//
//		Point3d p1(0,0,0);
//		Point3d p2(1,0,0);
//		Point3d p3(3,3,0);
//		Point3d p4(0,1,0);
//
//		m->AddVertex(p0);
//		m->AddVertex(p1);
//		m->AddVertex(p2);
//		m->AddVertex(p3);
//		m->AddVertex(p4);
//
//
//		/**
//		 \brief             
//		 p4............p3
//		 
//		 p1............p2
//		 \return	.
//		 */
//		Triangle t0(1,2,3);
//		Triangle t1(1,3,4);
//		m->AddFace(t0);
//		m->AddFace(t1);
//
//
//		return m;
//	}
//	static Mesh* MergeMesh(std::vector<Mesh*>& meshes)
//	{
//		std::vector<OringinalTriangle> otriangles;
//		int fsum=0;
//		for(size_t i=0;i<meshes.size();i++)
//			fsum+=meshes[i]->Faces.size();
//		otriangles.reserve(fsum);
//		for(size_t i=0;i<meshes.size();i++)
//		{
//			for(size_t j=0;j<meshes[i]->Faces.size();j++)
//			{
//				Triangle& t=meshes[i]->Faces[j];
//				Point3d& p0=meshes[i]->Vertices[t.P0Index];
//				Point3d& p1=meshes[i]->Vertices[t.P1Index];
//				Point3d& p2=meshes[i]->Vertices[t.P2Index];
//				OringinalTriangle ot(p0.X,p0.Y,p0.Z,p1.X,p1.Y,p1.Z,p2.X,p2.Y,p2.Z);
//				otriangles.push_back(ot);
//			}
//		}
//		return ConvertToMesh(otriangles);
//	}
//	static Mesh* ConvertToMesh(std::vector<OringinalTriangle>& otriangles)
//	{
//		Mesh *retMesh=new Mesh();
//		std::vector<FloatTriple> vec;
//		std::vector<Triangle> triangles;
//		triangles.reserve(otriangles.size());
//		vec.reserve(otriangles.size()*3);
//		for(size_t i=0;i<otriangles.size();i++)
//		{
//			FloatTriple ft0,ft1,ft2;
//			ft0.X=otriangles[i].P[0].X;
//			ft0.Y=otriangles[i].P[0].Y;
//			ft0.Z=otriangles[i].P[0].Z;
//			ft1.X=otriangles[i].P[1].X;
//			ft1.Y=otriangles[i].P[1].Y;
//			ft1.Z=otriangles[i].P[1].Z;
//			ft2.X=otriangles[i].P[2].X;
//			ft2.Y=otriangles[i].P[2].Y;
//			ft2.Z=otriangles[i].P[2].Z;
//			int count=vec.size();
//			ft0.Index=count;
//			ft1.Index=count+1;
//			ft2.Index=count+2;
//			vec.push_back(ft0);
//			vec.push_back(ft1);
//			vec.push_back(ft2);
//			Triangle t(count,count+1,count+2);
//			triangles.push_back(t);
//		}
//		std::vector<int> map1,map2;
//		map1.resize(otriangles.size()*3);
//		map2.resize(otriangles.size()*3);
//		std::sort(vec.begin(),vec.end());
//		for(size_t i=0;i<vec.size();i++)
//			map1[vec[i].Index]=i;
//		int index=0;
//		std::vector<Point3d> points;
//		Point3d p0(vec[0].X,vec[0].Y,vec[0].Z);
//		points.push_back(p0);
//		for(size_t i=1;i<vec.size();i++)
//		{
//			if(!(vec[i]<vec[i-1])&&!(vec[i-1]<vec[i]))
//			{
//				map2[i]=index;
//			}
//			else
//			{
//				index++;
//				map2[i]=index;
//				Point3d pi(vec[i].X,vec[i].Y,vec[i].Z);
//				points.push_back(pi);
//			}
//		}
//		for(size_t i=0;i<triangles.size();i++)
//		{
//			triangles[i].P0Index=map2[map1[triangles[i].P0Index]];
//			triangles[i].P1Index=map2[map1[triangles[i].P1Index]];
//			triangles[i].P2Index=map2[map1[triangles[i].P2Index]];
//		}
//		retMesh->Vertices.swap(points);
//		retMesh->Faces.swap(triangles);
//		return retMesh;
//	}
//	static void CleanMeshTVertices(Mesh& m)
//	{
//		TetraMesh tm;
//		tm.Vertices.swap(m.Vertices);
//		tm.Faces.swap(m.Faces);
//		tm.InitSurfaceEdges();
//		std::vector<int> emark;
//		std::vector<int> vmark;
//		vmark.resize(tm.Vertices.size(),0);
//		emark.resize(tm.Edges.size(),0);
//		std::vector<std::vector<int>> nv;
//		nv.resize(tm.Vertices.size());
//		for(size_t i=0;i<tm.Edges.size();i++)
//		{
//			int nc=tm.NeighborTrianglePerEdge[i].size();
//			if(nc==1)
//			{
//				emark[i]=-1;
//				vmark[tm.Edges[i].PIndex[0]]=-1;
//				vmark[tm.Edges[i].PIndex[1]]=-1;
//				nv[tm.Edges[i].PIndex[0]].push_back(tm.Edges[i].PIndex[1]);
//				nv[tm.Edges[i].PIndex[1]].push_back(tm.Edges[i].PIndex[0]);
//			}
//			if(nc==0)
//				throw std::exception();
//		}
//		std::queue<int> queue; 
//		int rec=1;
//		for(size_t i=0;i<tm.Vertices.size();i++)
//		{
//			if(vmark[i]==-1)
//			{
//				vmark[i]=rec;
//				queue.push(i);
//				while(!queue.empty())
//				{
//					int index=queue.front();
//					queue.pop();
//					std::vector<int> &neigbors=nv[index];
//					for(size_t k=0;k<neigbors.size();k++)
//					{
//						if(vmark[neigbors[k]]==-1)
//						{		
//							vmark[neigbors[k]]=rec;
//							queue.push(neigbors[k]);				
//						}
//					}
//				}
//				rec++;
//			}
//		}
//		std::vector<TRegion> regions;
//		regions.resize(rec);
//		for(size_t i=0;i<tm.Vertices.size();i++)
//		{
//			regions[vmark[i]].Id=vmark[i];
//			regions[vmark[i]].pindices.push_back(i);
//			regions[vmark[i]].regionBox.UpdateRange(tm.Vertices[i].X,tm.Vertices[i].Y,tm.Vertices[i].Z);
//		}
//		for(size_t i=0;i<regions.size();i++)
//		{
//			regions[i].points=&tm.Vertices;
//			Box3LFloat& box=regions[i].regionBox;
//			double boxv=box.GetXLength()*box.GetYLength()*box.GetZLength();
//			regions[i].isTHole=true;
//			regions[i].InitEnd();
//		    Point3d p=regions[i].GetCenter();
//			printf("region:[%d,center[%f,%f,%f],count:%d]\n",i,p.X,p.Y,p.Z,regions[i].pindices.size());
//		}
//		std::vector<Color> colors;
//		colors.resize(tm.Vertices.size());
//		for(size_t i=0;i<tm.Vertices.size();i++)
//		{
//			if(vmark[i]!=0)
//				colors[i]=IdToColor(vmark[i]);
//		}
//		m.Vertices.swap(tm.Vertices);
//		m.Faces.swap(tm.Faces);
//		PlyManager::Output_C(m,colors,"hole_dec.ply");
//	}
//	static Color IdToColor(int index)
//	{
//		index=(index*23)*(index*index)%255;
//		return GrayToC256((unsigned char)index);
//	}
//	static Color GrayToC256(unsigned char f)
//	{
//		unsigned char r = 0, g = 0, b = 0;
//		if (f >= 0 && f < 63)
//		{
//			r = 0; g = (unsigned char)(254 - 4 * f); b = 255;
//		}
//		if (f >= 64 && f < 127)
//		{
//			r = 0; g = (unsigned char)(4 * f - 254); b = (unsigned char)(510 - 4 * f);
//		}
//		if (f >= 128 && f <= 191)
//		{
//			r = (unsigned char)(4 * f - 510);
//			g = (unsigned char)(255);
//			b = (unsigned char)0;
//		}
//		if (f >= 192 && f <= 255)
//		{
//			r = (unsigned char)255;
//			g = (unsigned char)(1022 - 4 * f);
//			b = (unsigned char)0;
//		}
//		return Color(r, g, b);
//	}
//};
//class ConstrainedMeshSmoother
//{
//private:
//	Mesh* mesh;
//	int* marker;
//public:
//	ConstrainedMeshSmoother(Mesh* m)
//	{
//		this->mesh=m;
//		this->marker=new int[mesh->Vertices.size()];
//		memset(marker,0,sizeof(int)*mesh->Vertices.size());
//	}
//	~ConstrainedMeshSmoother()
//	{
//		delete[] marker;
//	}
//	void ExecuteConstrainedSmoothing(int iteration)
//	{
//		size_t vcount=mesh->Vertices.size();
//		size_t fcount=mesh->Faces.size();
//		double* vd1array = new (nothrow)double[vcount];
//		if(vd1array==NULL){return;}
//		double* vd2array = new (nothrow)double[vcount];
//		if(vd2array==NULL){delete[] vd1array;return;}
//		double* vd3array = new (nothrow)double[vcount];
//		if(vd3array==NULL){delete[] vd1array;delete[] vd2array;return;}
//		short* numbers=new (nothrow)short[vcount];
//		if(numbers==NULL){delete[] vd1array;delete[] vd2array;delete[] vd3array;return;}
//		if(vd1array==NULL||vd2array==NULL||vd3array==NULL||numbers==NULL){return;}
//		for(int c=0;c<iteration;c++)
//		{
//			memset(vd1array,0,sizeof(double)*vcount);
//			memset(vd2array,0,sizeof(double)*vcount);
//			memset(vd3array,0,sizeof(double)*vcount);
//			memset(numbers,0,sizeof(short)*vcount);
//			for(size_t i=0;i<fcount;i++)
//			{
//				Triangle& t=mesh->Faces[i];
//				Point3d p0=mesh->Vertices[t.P0Index];
//				Point3d p1=mesh->Vertices[t.P1Index];
//				Point3d p2=mesh->Vertices[t.P2Index];
//				vd1array[t.P0Index]+=p1.X;
//				vd2array[t.P0Index]+=p1.Y;
//				vd3array[t.P0Index]+=p1.Z;
//				vd1array[t.P0Index]+=p2.X;
//				vd2array[t.P0Index]+=p2.Y;
//				vd3array[t.P0Index]+=p2.Z;
//
//				vd1array[t.P1Index]+=p0.X;
//				vd2array[t.P1Index]+=p0.Y;
//				vd3array[t.P1Index]+=p0.Z;
//				vd1array[t.P1Index]+=p2.X;
//				vd2array[t.P1Index]+=p2.Y;
//				vd3array[t.P1Index]+=p2.Z;
//
//				vd1array[t.P2Index]+=p1.X;
//				vd2array[t.P2Index]+=p1.Y;
//				vd3array[t.P2Index]+=p1.Z;
//				vd1array[t.P2Index]+=p0.X;
//				vd2array[t.P2Index]+=p0.Y;
//				vd3array[t.P2Index]+=p0.Z;
//
//				numbers[t.P0Index]+=2;
//				numbers[t.P1Index]+=2;
//				numbers[t.P2Index]+=2;
//			}
//			for(size_t i=0;i<vcount;i++)
//			{
//				if(numbers[i]!=0&&marker[i]==0)
//				{
//					mesh->Vertices[i].X=vd1array[i]/numbers[i];
//					mesh->Vertices[i].Y=vd2array[i]/numbers[i];
//					mesh->Vertices[i].Z=vd3array[i]/numbers[i];
//				}
//			}
//		}
//		delete[] vd1array;
//		delete[] vd2array;
//		delete[] vd3array;
//		delete[] numbers;
//	}
//	void InitBoundaryMaker()
//	{
//		int* vn = new int[mesh->Vertices.size()];
//		int* fn = new int[mesh->Vertices.size()];
//		memset(marker,0,sizeof(int)*mesh->Vertices.size());
//		memset(vn,0,sizeof(int)*mesh->Vertices.size());
//		memset(fn,0,sizeof(int)*mesh->Vertices.size());
//		for(size_t i=0;i<mesh->Faces.size();i++)
//		{
//			Triangle& t=mesh->Faces[i];
//			fn[t.P0Index]++;
//			fn[t.P1Index]++;
//			fn[t.P2Index]++;
//			vn[t.P0Index]+=2;
//			vn[t.P1Index]+=2;
//			vn[t.P2Index]+=2;
//		}
//		int count=0;
//		for(size_t i=0;i<mesh->Vertices.size();i++)
//		{
//			if(2*fn[i]!=vn[i])
//			{
//				count++;
//				marker[i]=1;
//			}
//		}
//		printf("%d",count);
//		delete[] vn;
//		delete[] fn;
//	}
//	void SetMarker(int index,int value)
//	{
//		marker[index]=value;
//	}
//private:
//
//};
#endif