#ifndef CONSTRAINTTETRALIZATION_H
#define CONSTRAINTTETRALIZATION_H
#define TETLIBRARY 1
#include <queue>
#include "tetgen.h"
#include "vol_math_Mesh.h"
class RandomVertexGenerator
{
private:
	Box3LFloat box;
	double unitLen;
public:
	RandomVertexGenerator(double unitLength,Box3LFloat box)
	{
		this->unitLen=unitLength;
		this->box=box;
	}
	~RandomVertexGenerator(){}
public:
	int GenerateInnerPoints(std::vector<Point3d>& plist,std::vector<Int16Triple>& typelist,Int16Triple v)
	{
		int count=0;
		int xlenCount=(int)(box.GetXLength()/unitLen);
		int ylenCount=(int)(box.GetYLength()/unitLen);
		int zlenCount=(int)(box.GetZLength()/unitLen);
		double xlen=box.GetXLength()/xlenCount;
		double ylen=box.GetYLength()/ylenCount;
		double zlen=box.GetZLength()/zlenCount;

		for(int i=0;i<xlenCount-1;i++)
		{
			double xmin=(i+0.5f)*xlen+box.Min3[0];
			double xmax=((xmin+xlen)>box.Max3[0]?(box.Max3[0]):(xmin+xlen));
			for(int j=0;j<ylenCount-1;j++)
			{
				double ymin=(j+0.5f)*ylen+box.Min3[1];
				double ymax=((ymin+ylen)>box.Max3[1]?(box.Max3[1]):(ymin+ylen));
				for (int k=0;k<zlenCount-1;k++)
				{
					double zmin=(k+0.5f)*zlen+box.Min3[2];
					double zmax=((zmin+zlen)>box.Max3[2]?(box.Max3[2]):(zmin+zlen));
					Point3d p=GetClusterRandomPoints(xmin,xmax,ymin,ymax,zmin,zmax);
					plist.push_back(p);
					typelist.push_back(v);
					count++;
				}
			}
		}
		return count;
	}
	int GeneratePlanePoints(std::vector<Point3d>& plist,int dem,double value,std::vector<Int16Triple>& typelist,Int16Triple v)
	{
		int count=0;
		int xlenCount=0,ylenCount=0;
		double xlen=0,ylen=0,xst=0,xed=0,yst=0,yed=0;
		if(dem==2)
		{
			xlenCount=(int)(box.GetXLength()/unitLen);
			ylenCount=(int)(box.GetYLength()/unitLen);
			xlen=box.GetXLength()/xlenCount;
			ylen=box.GetYLength()/ylenCount;
			xst=box.Min3[0];
			xed=box.Max3[0];
			yst=box.Min3[1];
			yed=box.Max3[1];
		}
		if(dem==1)
		{
			xlenCount=(int)(box.GetXLength()/unitLen);
			ylenCount=(int)(box.GetZLength()/unitLen);
			xlen=box.GetXLength()/xlenCount;
			ylen=box.GetZLength()/ylenCount;
			xst=box.Min3[0];
			xed=box.Max3[0];
			yst=box.Min3[2];
			yed=box.Max3[2];
		}
		if(dem==0)
		{
			xlenCount=(int)(box.GetYLength()/unitLen);
			ylenCount=(int)(box.GetZLength()/unitLen);
			xlen=box.GetYLength()/xlenCount;
			ylen=box.GetZLength()/ylenCount;
			xst=box.Min3[1];
			xed=box.Max3[1];
			yst=box.Min3[2];
			yed=box.Max3[2];
		}

		for(int i=0;i<xlenCount-1;i++)
		{
			double xmin=(i+0.5f)*xlen+xst;
			double xmax=((xmin+xlen)>xed?(xed):(xmin+xlen));
			for(int j=0;j<ylenCount-1;j++)
			{
				double ymin=(j+0.5f)*ylen+yst;
				double ymax=((ymin+ylen)>yed?(yed):(ymin+ylen));
				double x=GetClusterRandom(xmin,xmax);
				double y=GetClusterRandom(ymin,ymax);
				if(dem==2)
				{
					Point3d p(x,y,value);plist.push_back(p);
					typelist.push_back(v);
					count++;
				}
				if(dem==1)
				{
					Point3d p(x,value,y);plist.push_back(p);
					typelist.push_back(v);
					count++;
				}
				if(dem==0)
				{
					Point3d p(value,x,y);plist.push_back(p);
					typelist.push_back(v);
					count++;
				}
			}
		}
		return count;
	}
	int GenerateBoxPoints(std::vector<Point3d>& plist,std::vector<Int16Triple>& typelist,Int16Triple v)
	{
		//double xlen=box.GetXLength();
		double x0=box.Min3[0],x1=box.Max3[0];
		//double ylen=box.GetYLength();
		double y0=box.Min3[1],y1=box.Max3[1];
		//double zlen=box.GetZLength();
		double z0=box.Min3[2],z1=box.Max3[2];
		Point3d p0(x0,y0,z0);
		Point3d p1(x1,y0,z0);
		Point3d p2(x0,y1,z0);
		Point3d p3(x1,y1,z0);
		Point3d p4(x0,y0,z1);
		Point3d p5(x1,y0,z1);
		Point3d p6(x0,y1,z1);
		Point3d p7(x1,y1,z1);
		plist.push_back(p0);
		plist.push_back(p1);
		plist.push_back(p2);
		plist.push_back(p3);
		plist.push_back(p4);
		plist.push_back(p5);
		plist.push_back(p6);
		plist.push_back(p7);
		for(int i=0;i<8;i++)
			typelist.push_back(v);
		return 8;
	}
	int GenerateLinePoints(std::vector<Point3d>& plist,int demParalTo,double d1,double d2,std::vector<Int16Triple>& typelist,Int16Triple v)
	{
		int count=0;
		int xlenCount=0;
		double xlen=0,xst=0,xed=0;
		if(demParalTo==0)
		{
			xlenCount=(int)(box.GetXLength()/unitLen);
			xlen=box.GetXLength()/xlenCount;
			xst=box.Min3[0];
			xed=box.Max3[0];
		}
		if(demParalTo==1)
		{
			xlenCount=(int)(box.GetYLength()/unitLen);
			xlen=box.GetYLength()/xlenCount;
			xst=box.Min3[1];
			xed=box.Max3[1];
		}
		if(demParalTo==2)
		{
			xlenCount=(int)(box.GetZLength()/unitLen);
			xlen=box.GetZLength()/xlenCount;
			xst=box.Min3[2];
			xed=box.Max3[2];
		}
		for(int i=0;i<xlenCount-1;i++)
		{
			double xmin=(i+0.5f)*xlen+xst;
			double xmax=((xmin+xlen)>xed?(xed):(xmin+xlen));
			if(demParalTo==0)
			{
				Point3d p(GetClusterRandom(xmin,xmax),d1,d2);
				plist.push_back(p);
				typelist.push_back(v);
				count++;
			}
			if(demParalTo==1)
			{
				Point3d p(d1,GetClusterRandom(xmin,xmax),d2);
				plist.push_back(p);
				typelist.push_back(v);
				count++;
			}
			if(demParalTo==2)
			{
				Point3d p(d1,d2,GetClusterRandom(xmin,xmax));
				plist.push_back(p);
				typelist.push_back(v);
				count++;
			}

		}
		return count;
	}
private:
	static double GetRandom(double st,double ed)
	{
		/*int r=rand()%1000;
		double x=r*(ed-st)/1000+st;
		return x;*/
		return (st+ed)/2.0f;
	}
	static double GetClusterRandom(double st,double ed,double clusterRate=0.35f)
	{
		if(clusterRate<0||clusterRate>=0.49f)
			return (st+ed)/2;
		else
		{
			double len=ed-st;
			st+=clusterRate*len;
			ed-=clusterRate*len;
			return GetRandom(st,ed);
		}
	}
	static Point3d GetClusterRandomPoints(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
	{
		return Point3d(GetClusterRandom(xmin,xmax), GetClusterRandom(ymin,ymax), GetClusterRandom(zmin,zmax) );
	}
};
class MassiveIndexedHash
{
private:
	std::vector<int> records;
	std::vector<int> table;
	int defaultValue;
public:
	MassiveIndexedHash(int size,int v)
	{
		table.resize(size,v);
		this->defaultValue=v;
	}
	~MassiveIndexedHash(){}
	void SetValue(int key,int v)
	{
		if(key<0||(size_t)key>table.size())
			throw std::exception();
		if(defaultValue!=v)
		{
			records.push_back(key);
			table[key]=v;
		}
	}
	int GetValue(int key)
	{
		if(key<0||(size_t)key>table.size())
			throw std::exception();
		return table[key];
	}
	void ClearTable()
	{
		for(size_t i=0;i<records.size();i++)
			table[records[i]]=defaultValue;
		records.clear();
	}
	int GetCount()
	{
		return (int)records.size();
	}
};
class VertexSpliter3d
{
	struct ConnectedVF
	{
		int pcenterIndex;
		std::vector<int> teds;
		int targetpIndex;
		ConnectedVF()
		{
			pcenterIndex=-1;
			targetpIndex=-1;
		}
	};
private:
	static inline int Max(int a,int b,int c)
	{
		int ab=(a>b?a:b);
		return ab>c?ab:c;
	}
	static inline int Min(int a,int b,int c)
	{
		int ab=(a<b?a:b);
		return ab<c?ab:c;
	}
	static void AppendPoints(std::vector<Point3d>&points,int count,Point3d v)
	{
		if(count==0)
			return;
		for(int i=0;i<count;i++)
		{
			points.push_back(v);
		}
	}
	static void ChangeTetraPoints(Tetrahedra* t,int orindex,int newindex)
	{
		for(int i=0;i<4;i++)
		{
			if(t->PIndex[i]==orindex)
			{
				t->PIndex[i]=newindex;
				return;
			}
		}
		printf("error!\n");
		throw std::exception();
	}
	static bool TriangleEqual(Triangle& t1,Triangle&t2)
	{
		int sum1=t1.P0Index+t1.P1Index+t1.P2Index;
		int sum2=t2.P0Index+t2.P1Index+t2.P2Index;
		if(sum1!=sum2)
			return false;
		int max1=Max(t1.P0Index,t1.P1Index,t1.P2Index);
		int max2=Max(t2.P0Index,t2.P1Index,t2.P2Index);
		if(max1!=max2)
			return false;
		int min1=Min(t1.P0Index,t1.P1Index,t1.P2Index);
		int min2=Min(t2.P0Index,t2.P1Index,t2.P2Index);
		if(min1!=min2)
			return false;
		return true;
	}
private:
	TetraMesh& mesh;
	std::vector<Triangle>& Constrains;
	std::vector<int> ConstrainsMapToTriangleIndex;
	std::vector<int> TriangleMapToConstrainsIndex;
	std::vector<bool> VBoundFlag;
	std::vector<int> VTypeFlag;
	std::vector<int> VRegionRes;
	std::vector<std::vector<int> > VFAdjs;
	std::vector<std::vector<int> > TTAdjs;
	std::vector<Triangle> InnerTriangles;
	std::vector<std::vector<int> > FTAdjs;
	std::vector<std::vector<int> > VTAdjs;
private:
	std::vector<bool> tempmark;
	std::vector<ConnectedVF> regions;
	std::queue<int> queue;
	MassiveIndexedHash hashtable;
public:
	VertexSpliter3d(TetraMesh& m,Mesh& constrains):mesh(m),Constrains(constrains.Faces),hashtable(m.Tedtrahedras.size(),-1)
	{
	}
	~VertexSpliter3d()
	{

	}
private:
	void InitVT()
	{
		VTAdjs.resize(mesh.Vertices.size());
		for(size_t i=0;i<mesh.Tedtrahedras.size();i++)
		{
			Tetrahedra& t=mesh.Tedtrahedras[i];
			if(VTypeFlag[t.PIndex[0]]>0)
				VTAdjs[t.PIndex[0]].push_back(i);
			if(VTypeFlag[t.PIndex[1]]>0)
				VTAdjs[t.PIndex[1]].push_back(i);
			if(VTypeFlag[t.PIndex[2]]>0)
				VTAdjs[t.PIndex[2]].push_back(i);
			if(VTypeFlag[t.PIndex[3]]>0)
				VTAdjs[t.PIndex[3]].push_back(i);
		}
	}
	void InitTypeFlags()
	{
		VTypeFlag.resize(mesh.Vertices.size(),0);
		for(size_t i=0;i<Constrains.size();i++)
		{
			VTypeFlag[Constrains[i].P0Index]++;
			VTypeFlag[Constrains[i].P1Index]++;
			VTypeFlag[Constrains[i].P2Index]++;
		}
		VRegionRes.resize(mesh.Vertices.size(),0);
	}
	void InitBoundFlags()
	{
		VBoundFlag.resize(mesh.Vertices.size(),false);
		for(size_t i=0;i<FTAdjs.size();i++)
		{
			if(FTAdjs[i].size()==1)
			{
				Triangle &t=InnerTriangles[i];
				VBoundFlag[t.P0Index]=true;
				VBoundFlag[t.P1Index]=true;
				VBoundFlag[t.P2Index]=true;
			}
		}
	}
	void InitInnerTriangle_FT()
	{
		if(mesh.InnerTriangles.size()!=0&&mesh.NeighborTetraPerTriangle.size()!=0)
		{
			this->FTAdjs.swap(mesh.NeighborTetraPerTriangle);
			this->InnerTriangles.swap(mesh.InnerTriangles);
		}
		else
		{
			mesh.InnerTriangles.clear();
			mesh.NeighborTetraPerTriangle.clear();
			mesh.InitInnerTriangles();
			this->FTAdjs.swap(mesh.NeighborTetraPerTriangle);
			this->InnerTriangles.swap(mesh.InnerTriangles);
		}
	}
	void InitVF()
	{
		VFAdjs.resize(mesh.Vertices.size());
		for(size_t i=0;i<InnerTriangles.size();i++)
		{
			VFAdjs[InnerTriangles[i].P0Index].push_back(i);
			VFAdjs[InnerTriangles[i].P1Index].push_back(i);
			VFAdjs[InnerTriangles[i].P2Index].push_back(i);
		}
	}
	void InitConstrainsMap()
	{
		TriangleMapToConstrainsIndex.resize(InnerTriangles.size(),-1);
		for(size_t i=0;i<Constrains.size();i++)
		{
			Triangle& con=Constrains[i];
			std::vector<int>& adjs=VFAdjs[con.P0Index];
			for(size_t j=0;j<adjs.size();j++)
			{
				int findex=adjs[j];
				Triangle &f=InnerTriangles[findex];
				if(TriangleEqual(con,f))
				{
					TriangleMapToConstrainsIndex[findex]=i;
				}
			}
		}
		ConstrainsMapToTriangleIndex.resize(Constrains.size(),-1);
		for(size_t i=0;i<TriangleMapToConstrainsIndex.size();i++)
		{
			int indexInCon=TriangleMapToConstrainsIndex[i];
			if(indexInCon!=-1)
			{
				ConstrainsMapToTriangleIndex[indexInCon]=i;
			}
		}
	}
	void InitTT()
	{
		TTAdjs.resize(mesh.Tedtrahedras.size());
		for(size_t i=0;i<InnerTriangles.size();i++)
		{
			std::vector<int>& fts=FTAdjs[i];
			if(fts.size()==2)
			{
				int ti1=fts[0];
				int ti2=fts[1];
				TTAdjs[ti1].push_back(ti2);
				TTAdjs[ti2].push_back(ti1);
			}
			if(fts.size()==0||fts.size()>2)
				throw std::exception();
		}
	}
	void InitTT_Type()
	{
		TTAdjs.resize(mesh.Tedtrahedras.size());
		for(size_t i=0;i<InnerTriangles.size();i++)
		{
			if(TriangleMapToConstrainsIndex[i]!=-1)
				continue;
			else
			{
				std::vector<int>& fts=FTAdjs[i];
				if(fts.size()==2)
				{
					int ti1=fts[0];
					int ti2=fts[1];
					TTAdjs[ti1].push_back(ti2);
					TTAdjs[ti2].push_back(ti1);
				}
				if(fts.size()==0||fts.size()>2)
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
	void InitRegions(int i,bool check=false)
	{
		regions.clear();regions.reserve(VTypeFlag[i]);
		std::vector<int>& allTetra=VTAdjs[i];
		tempmark.clear();tempmark.resize(allTetra.size(),false);
		hashtable.ClearTable();
		for(size_t j=0;j<allTetra.size();j++)
			hashtable.SetValue(allTetra[j],j);
		for(size_t j=0;j<allTetra.size();j++)
		{
			if(!tempmark[j])
			{
				ConnectedVF region;
				region.pcenterIndex=i;
				region.teds.push_back(allTetra[j]);
				queue.push(j);
				tempmark[j]=true;
				while(!queue.empty())
				{
					int k=queue.front();
					queue.pop();
					int tindex=allTetra[k];
					std::vector<int>& TTNeighb=TTAdjs[tindex];
					for(size_t n=0;n<TTNeighb.size();n++)
					{
						int nj=hashtable.GetValue(TTNeighb[n]);
						if(nj!=-1)
						{
							if(!tempmark[nj])
							{
								tempmark[nj]=true;
								region.teds.push_back(allTetra[nj]);
								queue.push(nj);
							}
						}
					}
				}
				regions.push_back(region);
			}
		}
		if(check)
		{
			size_t sum=0;
			for(size_t j=0;j<regions.size();j++)
			{
				sum+=regions[j].teds.size();
			}
			if(sum!=VTAdjs[i].size())
			{
				printf("region fail!\n");
				throw std::exception();
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
			std::vector<int>& facesInthisregion=regions[i].teds;
			for(size_t j=0;j<facesInthisregion.size();j++)
			{
				ChangeTetraPoints(&(mesh.Tedtrahedras[facesInthisregion[j]]),regions[i].pcenterIndex,regions[i].targetpIndex);
			}
		}
	}
	void ProcessSingleVertex(int i)
	{
		InitRegions(i,true);
		VRegionRes[i]=regions.size();
		AppendPoints(mesh.Vertices,regions.size()-1,mesh.Vertices[i]);
		MarshalNewPoints(i);
	}
public:
	void ProcessSpliting()
	{
		InitTypeFlags();
		InitVT();
		InitInnerTriangle_FT();
		InitBoundFlags();
		InitVF();
		InitConstrainsMap();
		InitTT_Type();
		mesh.Vertices.reserve(mesh.Vertices.size()+GetInsertedNum());
		size_t orcount=mesh.Vertices.size();
		for(size_t i=0;i<orcount;i++)
		{
			if(VTypeFlag[i]>0)
			{
				ProcessSingleVertex(i);
			}
		}
		mesh.InnerTriangles.swap(this->InnerTriangles);
		mesh.NeighborTetraPerTriangle.swap(this->FTAdjs);
		mesh.InnerTriangles.clear();
		mesh.NeighborTetraPerTriangle.clear();
		mesh.InitInnerTriangles();
	}
};
class TetGen
{
private:
	enum
	{
		X0=1,
		X1,
		Y0,
		Y1,
		Z0,
		Z1
	};
	static bool Equal(double a,double b)
	{
		double del=a-b>0?(a-b):(b-a);
		return del<0.00000001f;
	}
	static Point3d GetBoxCenter(Box3LFloat box)
	{
		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
		return Point3d(x0/2+x1/2,y0/2+y1/2,z0/2+z1/2);
	}
	static double GetLength(Point3d& p1,Point3d& p2)
	{
		return sqrt((p1.X-p2.X)*(p1.X-p2.X)+(p1.Y-p2.Y)*(p1.Y-p2.Y)+(p1.Z-p2.Z)*(p1.Z-p2.Z));
	}
	static double GetUlen(std::vector<Point3d>& points,std::vector<Triangle>& faces,double rate=2.0f)
	{
		double sum=0;
		for(size_t i=0;i<faces.size();i++)
		{
			double len1=GetLength(points[faces[i].P0Index],points[faces[i].P1Index]);
			double len2=GetLength(points[faces[i].P0Index],points[faces[i].P2Index]);
			double len3=GetLength(points[faces[i].P1Index],points[faces[i].P2Index]);
			sum+=len1;
			sum+=len2;
			sum+=len3;
			//printf("%f\n",len1);
		}
		return sum/(faces.size()*3);
	}
public:
	static void BuildInsertedPoints(Box3LFloat box,std::vector<Point3d>& insertedPoints,double ulen,std::vector<Int16Triple>& outFlags)
	{
		insertedPoints.clear();
		outFlags.clear();
		bool insertPoints_Inner=true;
		bool insertPoints_Plane=true;
		bool insertPoints_Line=true;
		bool insertPoints_Box=false;
		RandomVertexGenerator rvg(ulen,box);
		if(insertPoints_Inner)
		{
			rvg.GenerateInnerPoints(insertedPoints,outFlags,Int16Triple(1,1,1));
		}
		if(insertPoints_Plane)
		{
			rvg.GeneratePlanePoints(insertedPoints,0,box.Min3[0],outFlags,Int16Triple(0,1,1));
			rvg.GeneratePlanePoints(insertedPoints,0,box.Max3[0],outFlags,Int16Triple(0,1,1));
			rvg.GeneratePlanePoints(insertedPoints,1,box.Min3[1],outFlags,Int16Triple(1,0,1));
			rvg.GeneratePlanePoints(insertedPoints,1,box.Max3[1],outFlags,Int16Triple(1,0,1));
			rvg.GeneratePlanePoints(insertedPoints,2,box.Min3[2],outFlags,Int16Triple(1,1,0));
			rvg.GeneratePlanePoints(insertedPoints,2,box.Max3[2],outFlags,Int16Triple(1,1,0));
		}

		if(insertPoints_Line)
		{
			rvg.GenerateLinePoints(insertedPoints,0,box.Min3[1],box.Min3[2],outFlags,Int16Triple(1,0,0));
			rvg.GenerateLinePoints(insertedPoints,0,box.Max3[1],box.Min3[2],outFlags,Int16Triple(1,0,0));
			rvg.GenerateLinePoints(insertedPoints,0,box.Min3[1],box.Max3[2],outFlags,Int16Triple(1,0,0));
			rvg.GenerateLinePoints(insertedPoints,0,box.Max3[1],box.Max3[2],outFlags,Int16Triple(1,0,0));

			rvg.GenerateLinePoints(insertedPoints,1,box.Min3[0],box.Min3[2],outFlags,Int16Triple(0,1,0));
			rvg.GenerateLinePoints(insertedPoints,1,box.Max3[0],box.Min3[2],outFlags,Int16Triple(0,1,0));
			rvg.GenerateLinePoints(insertedPoints,1,box.Min3[0],box.Max3[2],outFlags,Int16Triple(0,1,0));
			rvg.GenerateLinePoints(insertedPoints,1,box.Max3[0],box.Max3[2],outFlags,Int16Triple(0,1,0));

			rvg.GenerateLinePoints(insertedPoints,2,box.Min3[0],box.Min3[1],outFlags,Int16Triple(0,0,1));
			rvg.GenerateLinePoints(insertedPoints,2,box.Max3[0],box.Min3[1],outFlags,Int16Triple(0,0,1));
			rvg.GenerateLinePoints(insertedPoints,2,box.Min3[0],box.Max3[1],outFlags,Int16Triple(0,0,1));
			rvg.GenerateLinePoints(insertedPoints,2,box.Max3[0],box.Max3[1],outFlags,Int16Triple(0,0,1));
		}
		if(insertPoints_Box)
			rvg.GenerateBoxPoints(insertedPoints,outFlags,Int16Triple(1,1,1));
		if(insertedPoints.size()!=outFlags.size())
			throw std::exception();
	}
	static TetraMesh* Execute(Mesh& m,Box3LFloat box,bool insertion=false,bool innersmoothing=true,int iter=1)
	{
		tetgenio in,out;
		tetgenbehavior b;
		tetgenio addin;
		TetGen::InitTetgenIO(in,out,b);
		TetGen::InitPoints(in,m.Vertices,box);
		TetGen::InitFacets(in,m.Vertices,m.Faces,box);
		std::vector<Point3d> insertedPoints;
		std::vector<Int16Triple> insertedoptypes;
		if(insertion)
		{
			double ulen=GetUlen(m.Vertices,m.Faces);
			TetGen::BuildInsertedPoints(box,insertedPoints,ulen,insertedoptypes);
			TetGen::InitAddin(addin,insertedPoints);
		}
		printf("check point1 \n");
		tetrahedralize(&b,&in,&out,&addin);
		printf("check point2 \n");
		TetraMesh* outM= TetGen::Out2TetraMesh(out);
		printf("check point3 \n");
		outM->InitNeighborTetraPerVertex();
		printf("check point4 \n");
		outM->InitInnerTriangles();
		printf("check point5 \n");
		if(insertion&&innersmoothing)
		{
			ExecuteInnerSmoothing(*outM,m,iter,insertedoptypes);
			printf("check point5.5 \n");
		}
		printf("check point6 \n");
		return outM;
	}
	static TetraMesh* Tetrahedralize(std::vector<Point3d>& points)
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
		return Out2TetraMesh(out);
	}
	static TetraMesh* Out2TetraMesh(tetgenio& out)
	{
		TetraMesh* m=new TetraMesh();
		for(int i=0;i<out.numberofpoints;i++)
		{
			Point3d p((double)(out.pointlist[3*i]),(double)(out.pointlist[3*i+1]),(double)(out.pointlist[3*i+2]));
			m->AddVertex(p);
		}
		for(int i=0;i<out.numberoftrifaces;i++)
		{
			Triangle t(out.trifacelist[3*i],out.trifacelist[3*i+1],out.trifacelist[3*i+2]);
			m->AddFace(t);
		}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			Tetrahedra te(out.tetrahedronlist[4*i],out.tetrahedronlist[4*i+1],out.tetrahedronlist[4*i+2],out.tetrahedronlist[4*i+3]);
			m->AddTetra(te);
		}

		return m;
	}
	static void InitTetgenIO(tetgenio& in,tetgenio& out, tetgenbehavior& b)
	{
		//b.quality=1;
		//b.regionattrib=1;
		//b.varvolume=1;

		//b.parse_commandline("-nn");
		b.plc=1;
		b.insertaddpoints=1;
		b.nobisect=1;
		b.verbose=1;
		//b.
	}
	static void InitAddin(tetgenio& addin,std::vector<Point3d>& insertedPoints)
	{
		addin.numberofpoints=insertedPoints.size();
		addin.pointlist=new REAL[insertedPoints.size()*3];
		for(size_t i=0;i<insertedPoints.size();i++)
		{
			addin.pointlist[3*i]=insertedPoints[i].X;
			addin.pointlist[3*i+1]=insertedPoints[i].Y;
			addin.pointlist[3*i+2]=insertedPoints[i].Z;
		}
	}
	static void InitPoints(tetgenio& in,std::vector<Point3d>& points,Box3LFloat box)
	{
		in.numberofpoints=points.size()+8;
		in.pointlist=new REAL[in.numberofpoints*3];
		for(size_t i=0;i<points.size();i++)
		{
			in.pointlist[3*i]=points[i].X;
			in.pointlist[3*i+1]=points[i].Y;
			in.pointlist[3*i+2]=points[i].Z;
		}
		InitBoxPoints(in,points,box);
	}
	static void InitFacets(tetgenio& in,std::vector<Point3d>& points,std::vector<Triangle>& faces,Box3LFloat box)
	{
			in.firstnumber = 0;
			in.numberoffacets = faces.size()+6;
			in.facetlist = new tetgenio::facet[in.numberoffacets];
			in.facetmarkerlist = new int[in.numberoffacets];
			memset(in.facetmarkerlist,0,sizeof(int)*in.numberoffacets);
			tetgenio::facet *f;tetgenio::polygon *p;
			for(size_t i=0;i<faces.size();i++)
			{
				f = &in.facetlist[i];
				f->numberofpolygons = 1;
				f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
				f->numberofholes = 0;
				f->holelist = NULL;
				p = &f->polygonlist[0];
				p->numberofvertices = 3;
				p->vertexlist = new int[p->numberofvertices];
				p->vertexlist[0] = faces[i].P0Index;
				p->vertexlist[1] = faces[i].P1Index;
				p->vertexlist[2] = faces[i].P2Index;
			}
			std::vector<std::vector<Int16Double> > segmentsInbound;
			segmentsInbound.resize(6);
			std::vector<std::vector<int> > pointsInbound;
			pointsInbound.resize(12);
			InitBoxEdgePoints(in,points,faces,box,pointsInbound);
			InitBoxSegments(in,points,faces,box,segmentsInbound);
			InitBoxFacets(in,points,faces,box,segmentsInbound,pointsInbound);
	}
	static void InitBoxPoints(tetgenio& in,std::vector<Point3d>& points,Box3LFloat box)
	{
		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
		int count=points.size()*3;
		in.pointlist[count]=x0;in.pointlist[count+1]=y0;in.pointlist[count+2]=z0;count+=3;
		in.pointlist[count]=x1;in.pointlist[count+1]=y0;in.pointlist[count+2]=z0;count+=3;
		in.pointlist[count]=x1;in.pointlist[count+1]=y1;in.pointlist[count+2]=z0;count+=3;
		in.pointlist[count]=x0;in.pointlist[count+1]=y1;in.pointlist[count+2]=z0;count+=3;
		in.pointlist[count]=x0;in.pointlist[count+1]=y0;in.pointlist[count+2]=z1;count+=3;
		in.pointlist[count]=x1;in.pointlist[count+1]=y0;in.pointlist[count+2]=z1;count+=3;
		in.pointlist[count]=x1;in.pointlist[count+1]=y1;in.pointlist[count+2]=z1;count+=3;
		in.pointlist[count]=x0;in.pointlist[count+1]=y1;in.pointlist[count+2]=z1;count+=3;
	}
	static void InitBoxFacets(tetgenio& in,std::vector<Point3d>& points,std::vector<Triangle>& faces,Box3LFloat box,std::vector<std::vector<Int16Double> >& segInBox,std::vector<std::vector<int> >& pb)
	{
		tetgenio::facet *f;
		tetgenio::polygon *p;
		std::vector<Int16Double>* segs;
		int count=faces.size();
		in.facetmarkerlist[count]=-1;
		in.facetmarkerlist[count+1]=-2;
		in.facetmarkerlist[count+2]=-3;
		in.facetmarkerlist[count+3]=-4;
		in.facetmarkerlist[count+4]=-5;
		in.facetmarkerlist[count+5]=-6;

		int p1=in.numberofpoints-8;
		int p2=in.numberofpoints-8+1;
		int p3=in.numberofpoints-8+2;
		int p4=in.numberofpoints-8+3;
		int p5=in.numberofpoints-8+4;
		int p6=in.numberofpoints-8+5;
		int p7=in.numberofpoints-8+6;
		int p8=in.numberofpoints-8+7;

		//if(p.X==x0&&p.Y==y0)  1 5    pointsInbound[0]
		//if(p.X==x1&&p.Y==y0)  2 6    pointsInbound[1]
		//if(p.X==x1&&p.Y==y1)  3 7     pointsInbound[2]
		//if(p.X==x0&&p.Y==y1)  4 8    pointsInbound[3]
		//
		//if(p.Y==y0&&p.Z==z0)  1 2    pointsInbound[4]
		//if(p.Y==y1&&p.Z==z0)  3 4     pointsInbound[5]
		//if(p.Y==y1&&p.Z==z1)  7 8    pointsInbound[6]
		//if(p.Y==y0&&p.Z==z1)  5 6    pointsInbound[7]
		//
		//if(p.X==x0&&p.Z==z0)  1 4    pointsInbound[8]
		//if(p.X==x1&&p.Z==z0)  2 3    pointsInbound[9]
		//if(p.X==x1&&p.Z==z1)  6 7    pointsInbound[10]
		//if(p.X==x0&&p.Z==z1)  5 8    pointsInbound[11]
		// Facet  x=x0
		f = &in.facetlist[count];
		segs=&(segInBox[0]);
		f->numberofpolygons =1+segs->size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		/*p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[1] =p5;
		p->vertexlist[2] =p8;
		p->vertexlist[3] =p4;
		p->vertexlist[4] =p1;*/
		ProcessBoxIntersectedPoints(p1,pb[0],true,p5,pb[11],true,p8,pb[3],false,p4,pb[8],false,p,points);
		ProcessBoxSegments(f,p,segs);


		// Facet x=x2
		f = &in.facetlist[count+1];
		segs=&(segInBox[1]);
		f->numberofpolygons =1+segs->size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		/*p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = p2;
		p->vertexlist[1] = p3;
		p->vertexlist[2] = p7;
		p->vertexlist[3] = p6;*/
		ProcessBoxIntersectedPoints(p2,pb[9],true,p3,pb[2],true,p7,pb[10],false,p6,pb[1],false,p,points);
		ProcessBoxSegments(f,p,segs);

		// Facet  y=y0
		f = &in.facetlist[count+2];
		segs=&(segInBox[2]);
		f->numberofpolygons =1+segs->size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		//p->numberofvertices = 4;
		//p->vertexlist = new int[p->numberofvertices];
		//p->vertexlist[0] = p1;
		//p->vertexlist[1] = p2;
		//p->vertexlist[2] = p6;
		//p->vertexlist[3] = p5;
		ProcessBoxIntersectedPoints(p1,pb[4],true,p2,pb[1],true,p6,pb[7],false,p5,pb[0],false,p,points);
		ProcessBoxSegments(f,p,segs);

		// Facet  y=y2
		f = &in.facetlist[count+3];
		segs=&(segInBox[3]);
		f->numberofpolygons =1+segs->size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
	/*	p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = p3;
		p->vertexlist[1] = p4;
		p->vertexlist[2] = p8;
		p->vertexlist[3] = p7;*/
		ProcessBoxIntersectedPoints(p3,pb[5],false,p4,pb[3],true,p8,pb[6],true,p7,pb[2],false,p,points);
		ProcessBoxSegments(f,p,segs);

		// Facet z=z0
		f = &in.facetlist[count+4];
		segs=&(segInBox[4]);
		f->numberofpolygons =1+segs->size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		//p->numberofvertices = 4;
		//p->vertexlist = new int[p->numberofvertices];
		//p->vertexlist[0] = p1;
		//p->vertexlist[1] = p2;
		//p->vertexlist[2] = p3;
		//p->vertexlist[3] = p4;
		ProcessBoxIntersectedPoints(p1,pb[4],true,p2,pb[9],true,p3,pb[5],false,p4,pb[8],false,p,points);
		ProcessBoxSegments(f,p,segs);

		// Facet  z=z2
		f = &in.facetlist[count+5];
		segs=&(segInBox[5]);
		f->numberofpolygons =1+segs->size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
	/*	p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = p5;
		p->vertexlist[1] = p6;
		p->vertexlist[2] = p7;
		p->vertexlist[3] = p8;*/
		ProcessBoxIntersectedPoints(p5,pb[7],true,p6,pb[10],true,p7,pb[6],false,p8,pb[11],false,p,points);
		ProcessBoxSegments(f,p,segs);

	}
	static void InitBoxEdgePoints(tetgenio& in,std::vector<Point3d>& points,std::vector<Triangle>& faces,Box3LFloat box,std::vector<std::vector<int> >& boundPoints)
	{
		int ret[12];
		for(size_t i=0;i<points.size();i++)
		{
			if(!InBox(points[i],box))
			{
				printf("%f %f %f",points[i].X,points[i].Y,points[i].Z);
				//throw new exception();
			}
			PointAlignBoxEdge(points[i],box,ret);
			for(int j=0;j<12;j++)
			{
				if(ret[j]==1)
				{
					boundPoints[j].push_back(i);
				}
			}
		}
	}
	static void InitBoxSegments(tetgenio& in,std::vector<Point3d>& points,std::vector<Triangle>& faces,Box3LFloat box,std::vector<std::vector<Int16Double> >& boundSegs)
	{
		int ret[6];
		for(size_t i=0;i<faces.size();i++)
		{
			Triangle& t=faces[i];
			Point3d& p0=points[t.P0Index];
			Point3d& p1=points[t.P1Index];
			Point3d& p2=points[t.P2Index];
			if(!InBox(p0,box)||!InBox(p1,box)||!InBox(p2,box))
				throw std::exception();
			SegAlignBoxFace(p0,p1,box,ret);
			for(int j=0;j<6;j++)
			{
				if(ret[j]==1)
				{
					Int16Double seg(t.P0Index,t.P1Index);
					boundSegs[j].push_back(seg);
				}
			}
			SegAlignBoxFace(p1,p2,box,ret);
			for(int j=0;j<6;j++)
			{
				if(ret[j]==1)
				{
					Int16Double seg(t.P1Index,t.P2Index);
					boundSegs[j].push_back(seg);
				}
			}
			SegAlignBoxFace(p2,p0,box,ret);
			for(int j=0;j<6;j++)
			{
				if(ret[j]==1)
				{
					Int16Double seg(t.P2Index,t.P0Index);
					boundSegs[j].push_back(seg);
				}
			}
		}
	}
	static bool InBox(Point3d& p,Box3LFloat& box)
	{
		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
		return p.X>=x0&&p.X<=x1&&p.Y>=y0&&p.Y<=y1&&p.Z>=z0&&p.Z<=z1;
	}
	static bool IsBoxVertex(Point3d& p,Box3LFloat& box)
	{
		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
		if(p.X!=x0&&p.X!=x1)
			return false;
		if(p.Y!=y0&&p.Y!=y1)
			return false;
		if(p.Z!=z0&&p.Z!=z1)
			return false;
		return true;
	}
	static void SegAlignBoxFace(Point3d& pa,Point3d& pb,Box3LFloat& box,int ret[])
	{
		memset(ret,0,sizeof(int)*6);
		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
		if(pa.X==x0&&pb.X==x0)
			ret[0]=1;
		if(pa.X==x1&&pb.X==x1)
			ret[1]=1;
		if(pa.Y==y0&&pb.Y==y0)
			ret[2]=1;
		if(pa.Y==y1&&pb.Y==y1)
			ret[3]=1;
		if(pa.Z==z0&&pb.Z==z0)
			ret[4]=1;
		if(pa.Z==z1&&pb.Z==z1)
			ret[5]=1;
	}
	static void PointAlignBoxEdge(Point3d& p,Box3LFloat& box,int ret[])
	{
		memset(ret,0,sizeof(int)*12);
		double x0=box.Min3[0],y0=box.Min3[1],z0=box.Min3[2],x1=box.Max3[0],y1=box.Max3[1],z1=box.Max3[2];
	    if(IsBoxVertex(p,box))
			throw new exception();
		if(p.X==x0&&p.Y==y0)
			ret[0]=1;
		if(p.X==x1&&p.Y==y0)
			ret[1]=1;
		if(p.X==x1&&p.Y==y1)
			ret[2]=1;
		if(p.X==x0&&p.Y==y1)
			ret[3]=1;

		if(p.Y==y0&&p.Z==z0)
			ret[4]=1;
		if(p.Y==y1&&p.Z==z0)
			ret[5]=1;
		if(p.Y==y1&&p.Z==z1)
			ret[6]=1;
		if(p.Y==y0&&p.Z==z1)
			ret[7]=1;

		if(p.X==x0&&p.Z==z0)
			ret[8]=1;
		if(p.X==x1&&p.Z==z0)
			ret[9]=1;
		if(p.X==x1&&p.Z==z1)
			ret[10]=1;
		if(p.X==x0&&p.Z==z1)
			ret[11]=1;
	}
	static void ProcessBoxSegments(tetgenio::facet *f,tetgenio::polygon *p,std::vector<Int16Double>* segs)
	{
		for(int i=1;i<f->numberofpolygons;i++)
		{
			p = &f->polygonlist[i];
			p->numberofvertices = 2;
			p->vertexlist = new int[p->numberofvertices];
			p->vertexlist[0] =(*segs)[i-1].X;
			p->vertexlist[1] =(*segs)[i-1].Y;
		}
	}
	static void ProcessBoxIntersectedPoints(int fvid1,std::vector<int>& ins1,bool sorttype1,int fvid2,std::vector<int>& ins2,bool sorttype2,int fvid3,std::vector<int>& ins3,bool sorttype3,int fvid4,std::vector<int>& ins4,bool sorttype4,tetgenio::polygon *p,std::vector<Point3d>& points)
	{
		std::vector<int> ret;
		int index=0;
		p->numberofvertices=4+ins1.size()+ins2.size()+ins3.size()+ins4.size();
		p->vertexlist=new int[p->numberofvertices];

		p->vertexlist[index]=fvid1;
		index++;
		ret.clear();
		SortPoints(points,ins1,ret,sorttype1);
		for(size_t i=0;i<ret.size();i++)
			p->vertexlist[index+i]=ret[i];
		index+=ret.size();

		p->vertexlist[index]=fvid2;
		index++;
		ret.clear();
		SortPoints(points,ins2,ret,sorttype2);
		for(size_t i=0;i<ret.size();i++)
			p->vertexlist[index+i]=ret[i];
		index+=ret.size();

		p->vertexlist[index]=fvid3;
		index++;
		ret.clear();
		SortPoints(points,ins3,ret,sorttype3);
		for(size_t i=0;i<ret.size();i++)
			p->vertexlist[index+i]=ret[i];
		index+=ret.size();

		p->vertexlist[index]=fvid4;
		index++;
		ret.clear();
		SortPoints(points,ins4,ret,sorttype4);
		for(size_t i=0;i<ret.size();i++)
			p->vertexlist[index+i]=ret[i];
		index+=ret.size();

	}
	static void SortPoints(std::vector<Point3d> &points,std::vector<int>& ins,std::vector<int>& ret,bool type)
	{
		if(ret.size()!=0)
			throw new exception();
		ret.resize(ins.size());
		std::copy(ins.begin(),ins.end(),ret.begin());
		int n=ins.size();
		for(int j=0;j<n-1;j++)
		{
			for(int i=0;i<n-1-j;i++)
			{
				double sum1 = points[ret[i]].X+points[ret[i]].Y+points[ret[i]].Z;
				double sum2 = points[ret[i+1]].X+points[ret[i+1]].Y+points[ret[i+1]].Z;
				if(type)
				{
					if(sum1>sum2)
					{
						int temp=ret[i];
						ret[i]=ret[i+1];
						ret[i+1]=temp;
					}
				}
				else
				{
					if(sum1<sum2)
					{
						int temp=ret[i];
						ret[i]=ret[i+1];
						ret[i+1]=temp;
					}
				}
			}
		}
	}
	static void DoInsertionFilting(Mesh& m,std::vector<Point3d>& points,double thres,std::vector<Point3d>& result)
	{

	}
	static void CheckConstrains(std::vector<Triangle>& innertriangles,std::vector<Triangle>& surfaces,std::vector<int>& errors,void (*show_progress)(int total,int step)=0)
	{
		for(size_t i=0;i<surfaces.size();i++)
		{
			if(show_progress!=0)
			{
				show_progress(surfaces.size(),i);
			}
			Triangle t=surfaces[i];
			bool contains=false;
			for(size_t j=0;j<innertriangles.size();j++)
			{
				Triangle tj=innertriangles[j];
				if(TriangleEqual(t,tj))
				{
					contains=true;
					break;
				}
			}
			if(contains)
				continue;
			else
				errors.push_back(i);
		}
	}
	static bool TriangleEqual(Triangle t1,Triangle t2)
	{
		if(t1.P0Index==t2.P0Index&&t1.P1Index==t2.P1Index&&t1.P2Index==t2.P2Index)
			return true;
		IntIndexTriple it1(t1.P0Index,t1.P1Index,t1.P2Index);
		IntIndexTriple it2(t2.P0Index,t2.P1Index,t2.P2Index);
		it1.Sort();
		it2.Sort();
		if(it1.X==it2.X&&it1.Y==it2.Y&&it1.Z==it2.Z)
			return true;
		else
			return false;
	}
	static void ExecuteInnerSmoothing(TetraMesh& mesh,Mesh& surface,int iteration,std::vector<Int16Triple>& optypes)
	{
		std::vector<Int16Triple> alltypes;
		size_t fixcount=surface.Vertices.size()+8;
		if(fixcount+optypes.size()==mesh.Vertices.size())
		{
			alltypes.resize(mesh.Vertices.size());
			for(size_t i=0;i<mesh.Vertices.size();i++)
			{
				if(i<fixcount)
					alltypes[i]=Int16Triple(0,0,0);
				else
					alltypes[i]=optypes[i-fixcount];
			}
		}
		else
		{
			printf("woca!\n");
			throw std::exception();
		}
		std::vector<double> vd1array;
		std::vector<double> vd2array;
		std::vector<double> vd3array;
		std::vector<short> numbers;
		size_t vcount=mesh.Vertices.size();
		size_t tcount=mesh.Tedtrahedras.size();
		vd1array.resize(mesh.Vertices.size(),0.0f);
		vd2array.resize(mesh.Vertices.size(),0.0f);
		vd3array.resize(mesh.Vertices.size(),0.0f);
		numbers.resize(mesh.Vertices.size(),0);
		for(int c=0;c<iteration;c++)
		{
			for(size_t i=0;i<tcount;i++)
			{
				Tetrahedra& t=mesh.Tedtrahedras[i];
				Point3d p0=mesh.Vertices[t.PIndex[0]];
				Point3d p1=mesh.Vertices[t.PIndex[1]];
				Point3d p2=mesh.Vertices[t.PIndex[2]];
				Point3d p3=mesh.Vertices[t.PIndex[3]];
				int P0Index=t.PIndex[0];
				int P1Index=t.PIndex[1];
				int P2Index=t.PIndex[2];
				int P3Index=t.PIndex[3];

				vd1array[P0Index]+=p1.X;vd2array[P0Index]+=p1.Y;vd3array[P0Index]+=p1.Z;
				vd1array[P0Index]+=p2.X;vd2array[P0Index]+=p2.Y;vd3array[P0Index]+=p2.Z;
				vd1array[P0Index]+=p3.X;vd2array[P0Index]+=p3.Y;vd3array[P0Index]+=p3.Z;

				vd1array[P1Index]+=p0.X;vd2array[P1Index]+=p0.Y;vd3array[P1Index]+=p0.Z;
				vd1array[P1Index]+=p2.X;vd2array[P1Index]+=p2.Y;vd3array[P1Index]+=p2.Z;
				vd1array[P1Index]+=p3.X;vd2array[P1Index]+=p3.Y;vd3array[P1Index]+=p3.Z;

				vd1array[P2Index]+=p1.X;vd2array[P2Index]+=p1.Y;vd3array[P2Index]+=p1.Z;
				vd1array[P2Index]+=p0.X;vd2array[P2Index]+=p0.Y;vd3array[P2Index]+=p0.Z;
				vd1array[P2Index]+=p3.X;vd2array[P2Index]+=p3.Y;vd3array[P2Index]+=p3.Z;

				vd1array[P3Index]+=p1.X;vd2array[P3Index]+=p1.Y;vd3array[P3Index]+=p1.Z;
				vd1array[P3Index]+=p0.X;vd2array[P3Index]+=p0.Y;vd3array[P3Index]+=p0.Z;
				vd1array[P3Index]+=p2.X;vd2array[P3Index]+=p2.Y;vd3array[P3Index]+=p2.Z;

				numbers[P0Index]+=3;
				numbers[P1Index]+=3;
				numbers[P2Index]+=3;
				numbers[P3Index]+=3;
			}
			for(size_t i=0;i<vcount;i++)
			{
				if(numbers[i]!=0)
				{
					if(alltypes[i].X==1)
						mesh.Vertices[i].X=vd1array[i]/numbers[i];
					if(alltypes[i].Y==1)
						mesh.Vertices[i].Y=vd2array[i]/numbers[i];
					if(alltypes[i].Z==1)
						mesh.Vertices[i].Z=vd3array[i]/numbers[i];
				}
				vd1array[i]=0.0f;
				vd2array[i]=0.0f;
				vd3array[i]=0.0f;
				numbers[i]=0;
			}
		}
	}
};
class TetraRegionMarker
{
public:
	static Color IdToColor(int index)
	{
		index=(index*23)*(index*index)%255;
		return GrayToC256((unsigned char)index);
	}
	static Color GrayToC256(unsigned char f)
	{
		unsigned char r = 0, g = 0, b = 0;
		if (f >= 0 && f < 63)
		{
			r = 0; g = (unsigned char)(254 - 4 * f); b = 255;
		}
		if (f >= 64 && f < 127)
		{
			r = 0; g = (unsigned char)(4 * f - 254); b = (unsigned char)(510 - 4 * f);
		}
		if (f >= 128 && f <= 191)
		{
			r = (unsigned char)(4 * f - 510);
			g = (unsigned char)(255);
			b = (unsigned char)0;
		}
		if (f >= 192 && f <= 255)
		{
			r = (unsigned char)255;
			g = (unsigned char)(1022 - 4 * f);
			b = (unsigned char)0;
		}
		return Color(r, g, b);
	}
private:
	TetraMesh& mesh;
	std::vector<std::vector<int> > VTAdjs;
	std::queue<int> queue;
public:
	TetraRegionMarker(TetraMesh& m):mesh(m){}
	~TetraRegionMarker(){}
public:
	void MarkRegions(std::vector<int>& vertexMark)
	{
		vertexMark.clear();
		vertexMark.resize(mesh.Vertices.size(),-1);
		VTAdjs.resize(mesh.Vertices.size());
		for(size_t i=0;i<mesh.Tedtrahedras.size();i++)
		{
			Tetrahedra& te=mesh.Tedtrahedras[i];
			VTAdjs[te.PIndex[0]].push_back(i);
			VTAdjs[te.PIndex[1]].push_back(i);
			VTAdjs[te.PIndex[2]].push_back(i);
			VTAdjs[te.PIndex[3]].push_back(i);
		}
		int regionCount=0;
		for(size_t i=0;i<mesh.Vertices.size();i++)
		{
			if(vertexMark[i]==-1)
			{
				queue.push(i);
				vertexMark[i]=regionCount;
				while(!queue.empty())
				{
					int index=queue.front();
					queue.pop();
					std::vector<int>& adjs=VTAdjs[index];
					for(size_t j=0;j<adjs.size();j++)
					{
						Tetrahedra& te=mesh.Tedtrahedras[adjs[j]];
						for(int k=0;k<4;k++)
						{
							int adjindex=te.PIndex[k];
							if(vertexMark[adjindex]==-1)
							{
								queue.push(adjindex);
								vertexMark[adjindex]=regionCount;
							}
						}
					}
				}
				regionCount++;
			}
		}
	}
	void MarkRegions(std::vector<Color>& vertexColors)
	{
		std::vector<int> marks;
		Color defa(255,255,255);
		vertexColors.resize(mesh.Vertices.size(),defa);
		MarkRegions(marks);
		for(size_t i=0;i<marks.size();i++)
		{
			vertexColors[i]=IdToColor(marks[i]);
		}
	}
};
class ConstraintTetrahedralizer
{
private:
	//static timer;
public:
	struct Parameters
	{
		bool insertion;
		bool innersmoothing;
		bool iterationTime;
	};
public:
	static bool CheckDataValid(Mesh& m,Box3LFloat box)
	{
		return true;
	}
	static TetraMesh* Execute(Mesh& m,Box3LFloat box,Parameters& parms)
	{
		if(!CheckDataValid(m,box))
			return NULL;
		TetraMesh* tm=TetGen::Execute(m,box,parms.insertion,parms.innersmoothing,parms.iterationTime);
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
		return tm;
	};
};
#endif
