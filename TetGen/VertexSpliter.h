#ifndef VERTEXSPLITER_H
#define VERTEXSPLITER_H
#include "vol_math_Mesh.h"
#include "vol_math_ConstraintTetralization.h"
#define WIDTH 10
struct Slice
{
	std::vector<Point3d> Points;
	std::vector<Int16Double> Segments;
	int AddPoints(double x,double y)
	{
		Points.push_back(Point3d(x,y,0));
		return Points.size()-1;
	}
	int AddPoints(int x,int y)
	{
		Points.push_back(Point3d((double)x,(double)y,0));
		return Points.size()-1;
	}
	void MakeSegments(int hash[WIDTH][WIDTH])
	{
		for(size_t i=0;i<Points.size()-1;i++)
		{
			int index1=hash[(int)Points[i].X][(int)Points[i].Y];
			int index2=hash[(int)Points[i+1].X][(int)Points[i+1].Y];
			Int16Double seg(index1,index2);
			Segments.push_back(seg);
		}
	}
};
struct SampleData2d
{
	int hash[WIDTH][WIDTH];
	Mesh mesh;
	std::vector<Slice> Slices;
	std::vector<Int16Double> Segments;
	std::vector<Color> Colors;
	void GetSample()
	{
		for(int i=0;i<WIDTH;i++)
		{
			for (int j=0;j<WIDTH;j++)
			{
				Point3d p((double)i,(double)j,0);
				hash[i][j]=mesh.AddVertex(p);
			}
		}
		for(int i=0;i<WIDTH-1;i++)
		{
			for (int j = 0; j < WIDTH-1; j++)
			{
				Triangle t1(hash[i][j],hash[i+1][j],hash[i][j+1]);
				Triangle t2(hash[i][j+1], hash[i + 1][ j], hash[i+1][ j + 1]);
				mesh.AddFace(t1);
				mesh.AddFace(t2);
			}
		}

	}
	void GetSlices()
	{
		Colors.push_back(Color(255,0,0));
		Colors.push_back(Color(255,165,0));
		Colors.push_back(Color(255,255,0));
		Colors.push_back(Color(0,255,0));
		Colors.push_back(Color(0,127,255));
		Colors.push_back(Color(0,0,255));
		Colors.push_back(Color(139,0,255));
		Colors.push_back(Color(139,239,139));
		Slice slice1;
		slice1.AddPoints(7,1);
		slice1.AddPoints(7,2);
		slice1.AddPoints(7,3);
		slice1.AddPoints(7,4);
		slice1.AddPoints(7,5);
		slice1.AddPoints(7,6);
		slice1.AddPoints(7,7);
		slice1.MakeSegments(hash);
		Slice slice2;
		slice2.AddPoints(1,8);
		slice2.AddPoints(2,7);
		slice2.AddPoints(3,6);
		slice2.AddPoints(4,5);
		slice2.AddPoints(5,4);
		slice2.AddPoints(6,3);
		slice2.AddPoints(7,2);
		slice2.MakeSegments(hash);
		Slice slice3;
		for(int i=0;i<=9;i++)
			slice3.AddPoints(4,i);
		slice3.MakeSegments(hash);
		Slice slice4;
		for(int i=0;i<=9;i++)
			slice4.AddPoints(5,i);
		slice4.MakeSegments(hash);
		Slice slice5;
		for(int i=0;i<=5;i++)
			slice5.AddPoints(2,i);
		slice5.MakeSegments(hash);
		Slice slice6;
		for(int i=1;i<=7;i++)
			slice6.AddPoints(i,1);
		slice6.MakeSegments(hash);
		Slice slice7;
		slice7.AddPoints(7,5);
		slice7.AddPoints(8,4);
		slice7.AddPoints(9,3);
		slice7.MakeSegments(hash);
		Slice slice8;
		for(int i=1;i<=9;i++)
			slice8.AddPoints(i,5);
		slice8.MakeSegments(hash);
		Slices.push_back(slice1);
		Slices.push_back(slice2);
		Slices.push_back(slice3);
		Slices.push_back(slice4);
		Slices.push_back(slice5);
		Slices.push_back(slice6);
		Slices.push_back(slice7);
		Slices.push_back(slice8);
	}
	void GetSegments()
	{
		for(size_t i=0;i<Slices.size();i++)
		{
			for(size_t j=0;j<Slices[i].Segments.size();j++)
			{
				Segments.push_back(Slices[i].Segments[j]);
			}
		}
	}
	void Output()
	{
		std::vector<Color> cm;cm.resize(mesh.Vertices.size(),Color(255,255,255));
		for(size_t i=0;i<Slices.size();i++)
		{
			Color c=Colors[i];
			for(size_t j=0;j<Slices[i].Segments.size();j++)
			{
				int pi1=Slices[i].Segments[j].X;
				int pi2=Slices[i].Segments[j].Y;
				cm[pi1]=c;
				cm[pi2]=c;
			}
		}
		PlyManager::Output_C(mesh,cm,"VSsample.ply");
	}
};
class VertexSpliter2d
{
	struct ConnectedVF
	{
		int pcenterIndex;
		std::vector<int> faces;
		int targetpIndex;
		ConnectedVF()
		{
			pcenterIndex=-1;
			targetpIndex=-1;
		}
	};
private:
	static void AppendPoints(std::vector<Point3d>&points,int count,Point3d v)
	{
		if(count==0)
			return;
		for(int i=0;i<count;i++)
		{
			points.push_back(v);
		}
	}
	static void ChangeTrianglePoints(Triangle* t,int orindex,int newindex)
	{
		if(t->P0Index==orindex)
		{
			t->P0Index=newindex;
			return;
		}
		if(t->P1Index==orindex)
		{
			t->P1Index=newindex;
			return;
		}
		if(t->P2Index==orindex)
		{
			t->P2Index=newindex;
			return;
		}
	}
private:
	Mesh& mesh;
	std::vector<Int16Double>& Segments;
	std::vector<int> SegMapToEdgeIndex;
	std::vector<int> EdgeMapToSegIndex;
	std::vector<bool> VBoundFlag;
	std::vector<int> VTypeFlag;
	std::vector<std::vector<int> > VEAdjs;
	std::vector<std::vector<int> > FFAdjs;
	std::vector<Edge> Edges;
	std::vector<std::vector<int> > EFAdjs;
	std::vector<std::vector<int> > VFAdjs;
private:
	std::vector<bool> tempmark;
	std::vector<ConnectedVF> regions;
	std::queue<int> queue;
	MassiveIndexedHash hashtable;
public:
	VertexSpliter2d(SampleData2d& da)
		:mesh(da.mesh),Segments(da.Segments),hashtable(da.mesh.Faces.size(),-1)
	{
		VBoundFlag.resize(da.mesh.Vertices.size(),false);
		VTypeFlag.resize(da.mesh.Vertices.size(),0);
	}
	~VertexSpliter2d()
	{

	}
private:
	void InitVF()
	{
		VFAdjs.resize(mesh.Vertices.size());
		for(size_t i=0;i<mesh.Faces.size();i++)
		{
			Triangle& t=mesh.Faces[i];
			if(VTypeFlag[t.P0Index]>0)
				VFAdjs[t.P0Index].push_back(i);
			if(VTypeFlag[t.P1Index]>0)
				VFAdjs[t.P1Index].push_back(i);
			if(VTypeFlag[t.P2Index]>0)
				VFAdjs[t.P2Index].push_back(i);
		}
	}
	void InitTypeFlags()
	{
		VTypeFlag.resize(mesh.Vertices.size(),0);
		for(size_t i=0;i<Segments.size();i++)
		{
			VTypeFlag[Segments[i].X]++;
			VTypeFlag[Segments[i].Y]++;
		}
	}
	void InitBoundFlags()
	{
		VBoundFlag.resize(mesh.Vertices.size(),false);
		for(size_t i=0;i<EFAdjs.size();i++)
		{
			if(EFAdjs[i].size()==1)
			{
				Edge &e=Edges[i];
				VBoundFlag[e.PIndex[0]]=true;
				VBoundFlag[e.PIndex[1]]=true;
			}
		}
	}
	void InitEdge_EF()
	{
		std::vector<IntIndexDouble> templist;
		for(size_t i=0;i<mesh.Faces.size();i++)
		{
			Triangle& te=mesh.Faces[i];
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
		FFAdjs.resize(mesh.Faces.size());
		int uniquesize=1;
		Edge t0(templist[0].X,templist[0].Y);
		Edges.push_back(t0);
		std::vector<int> vec0;vec0.reserve(2);
		vec0.push_back(templist[0].Index);
		EFAdjs.push_back(vec0);
		for(size_t i=1;i<templist.size();i++)
		{
			if(!(templist[i]==templist[i-1]))
			{
				uniquesize++;
				Edge t(templist[i].X,templist[i].Y);
				Edges.push_back(t);
				std::vector<int> vec;vec.reserve(2);
				vec.push_back(templist[i].Index);
				EFAdjs.push_back(vec);
			}
			else
			{
				EFAdjs[uniquesize-1].push_back(templist[i].Index);
			}
		}
	}
	void InitVE()
	{
		VEAdjs.resize(mesh.Vertices.size());
		for(size_t i=0;i<Edges.size();i++)
		{
			VEAdjs[Edges[i].PIndex[0]].push_back(i);
			VEAdjs[Edges[i].PIndex[1]].push_back(i);
		}
	}
	void InitSegMap()
	{
		EdgeMapToSegIndex.resize(Edges.size(),-1);
		for(size_t i=0;i<Segments.size();i++)
		{
			Int16Double& seg=Segments[i];
			std::vector<int>& adjs=VEAdjs[seg.X];
			for(size_t j=0;j<adjs.size();j++)
			{
				int eindex=adjs[j];
				Edge &e=Edges[eindex];
				if(e.PIndex[0]==seg.Y||e.PIndex[1]==seg.Y)
				{
					EdgeMapToSegIndex[eindex]=i;
				}
			}
		}
		SegMapToEdgeIndex.resize(Segments.size(),-1);
		for(size_t i=0;i<EdgeMapToSegIndex.size();i++)
		{
			int indexInSeg=EdgeMapToSegIndex[i];
			if(indexInSeg!=-1)
			{
				SegMapToEdgeIndex[indexInSeg]=i;
			}
		}
	}
	void InitFF()
	{
		for(size_t i=0;i<Edges.size();i++)
		{
			std::vector<int>& ef=EFAdjs[i];
			if(ef.size()==2)
			{
				int fi1=ef[0];
				int fi2=ef[1];
				FFAdjs[fi1].push_back(fi2);
				FFAdjs[fi2].push_back(fi1);
			}
			if(ef.size()==0||ef.size()>2)
				throw std::exception();
		}
	}
	void InitFF_Type()
	{
		for(size_t i=0;i<Edges.size();i++)
		{
			if(EdgeMapToSegIndex[i]!=-1)
				continue;
			else
			{
				std::vector<int>& ef=EFAdjs[i];
				if(ef.size()==2)
				{
					int fi1=ef[0];
					int fi2=ef[1];
					FFAdjs[fi1].push_back(fi2);
					FFAdjs[fi2].push_back(fi1);
				}
				if(ef.size()==0||ef.size()>2)
					throw std::exception();
			}
		}
	}
	int GetInsertedNum()
	{
		int sum=0;
		for(size_t i=0;i<VTypeFlag.size();i++)
		{
			if(!VBoundFlag[i])
				sum+=(VTypeFlag[i]-1);
			else
				sum+=VTypeFlag[i];
		}
		return sum;
	}
	void InitRegions(int i)
	{
		regions.clear();regions.reserve(VTypeFlag[i]);
		std::vector<int>& allFace=VFAdjs[i];
		tempmark.clear();tempmark.resize(allFace.size(),false);
		hashtable.ClearTable();
		for(size_t j=0;j<allFace.size();j++)
			hashtable.SetValue(allFace[j],j);
		for(size_t j=0;j<allFace.size();j++)
		{
			if(!tempmark[j])
			{
				ConnectedVF region;
				region.pcenterIndex=i;
				region.faces.push_back(allFace[j]);
				queue.push(j);
				tempmark[j]=true;
				while(!queue.empty())
				{
					int k=queue.front();
					queue.pop();
					int findex=allFace[k];
					std::vector<int>& FFNeighb=FFAdjs[findex];
					for(size_t n=0;n<FFNeighb.size();n++)
					{
						int nj=hashtable.GetValue(FFNeighb[n]);
						if(nj!=-1)
						{
							if(!tempmark[nj])
							{
								tempmark[nj]=true;
								region.faces.push_back(allFace[nj]);
								queue.push(nj);
							}
						}
					}
				}
				regions.push_back(region);
			}
		}
	}
	void MarshalNewPoints(int index)
	{
		int orsize=mesh.Vertices.size()-regions.size()+1;
		regions[0].targetpIndex=index;
		for(size_t i=1;i<regions.size();i++)
		{
			regions[i].targetpIndex=orsize+i-1;
			std::vector<int>& facesInthisregion=regions[i].faces;
			for(size_t j=0;j<facesInthisregion.size();j++)
			{
				ChangeTrianglePoints(&(mesh.Faces[facesInthisregion[j]]),regions[i].pcenterIndex,regions[i].targetpIndex);
			}
		}
	}
	void ProcessSingleVertex(int i)
	{
		InitRegions(i);
		AppendPoints(mesh.Vertices,regions.size()-1,mesh.Vertices[i]);
		MarshalNewPoints(i);
		printf("%d done\n",i);
	}
public:
	void ProcessSpliting()
	{
		InitTypeFlags();
		InitVF();
		InitEdge_EF();
		InitBoundFlags();
		InitVE();
		InitSegMap();
		InitFF_Type();
		mesh.Vertices.reserve(mesh.Vertices.size()+GetInsertedNum());
		size_t orcount=mesh.Vertices.size();
		for(size_t i=0;i<orcount;i++)
		{
			if(VTypeFlag[i]>0)
			{
				ProcessSingleVertex(i);
			}
		}
	}
public:
	static void Test()
	{
		SampleData2d data;
		data.GetSample();
		data.GetSlices();
		data.GetSegments();
		VertexSpliter2d vs(data);
		vs.ProcessSpliting();
		data.Output();
	}
};


#endif
