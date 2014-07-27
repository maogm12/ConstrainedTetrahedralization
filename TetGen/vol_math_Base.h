#ifndef BASE_H
#define BASE_H
#include <vector>
#include <math.h>
//#define NULL 0
#define MESH_INDEXTYPE int
#define COORD_INDEXTYPE int
#define M_FLOAT_RADIUS 0.0000001f
using namespace std;
struct LFloatDouble
{
	double X;
	double Y;
	LFloatDouble(double x,double y):X(x),Y(y){}
	LFloatDouble():X(0),Y(0){}
};
struct Color
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
	Color(unsigned char r,unsigned char g,unsigned char b)
	{
		R=r;
		G=g;
		B=b;
	}
	Color(){}
};
struct Int16Double
{
public:
	int X;
	int Y;
	Int16Double(int x,int y):X(x),Y(y){}
	Int16Double():X(0),Y(0){}
};
struct Int16Triple
{
	COORD_INDEXTYPE X;
	COORD_INDEXTYPE Y;
	COORD_INDEXTYPE Z;
	Int16Triple(COORD_INDEXTYPE x,  COORD_INDEXTYPE y,  COORD_INDEXTYPE z):X(x),Y(y),Z(z){}
	Int16Triple():X(0),Y(0),Z(0){}
};
struct IntTriple
{
	int X;
	int Y;
	int Z;
	IntTriple( int x,  int y,  int z):X(x),Y(y),Z(z){}
	IntTriple():X(0),Y(0),Z(0){}
};
struct Box3LFloat
{
public:
	double Min3[3];
	double Max3[3];
	Box3LFloat()
	{
		Min3[0]=9999999.0f;
		Min3[1]=9999999.0f;
		Min3[2]=9999999.0f;
		Max3[0]=-9999999.0f;
		Max3[1]=-9999999.0f;
		Max3[2]=-9999999.0f;
	}
	~Box3LFloat()
	{

	}
	bool IsValid()
	{
		return Min3[0]<=Max3[0]&&Min3[1]<=Max3[1]&&Min3[2]<=Max3[2];
	}
	Box3LFloat(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
	{
		Min3[0] = minX;
		Min3[1] = minY;
		Min3[2] = minZ;
		Max3[0] = maxX;
		Max3[1] = maxY;
		Max3[2] = maxZ;
	}
	void UpdateRange(double x, double y, double z)
	{
		if (x < Min3[0])
			Min3[0] = x;
		if (y < Min3[1])
			Min3[1] = y;
		if (z < Min3[2])
			Min3[2] = z;

		if (x > Max3[0])
			Max3[0] = x;
		if (y > Max3[1])
			Max3[1] = y;
		if (z > Max3[2])
			Max3[2] = z;
	}
	double GetXLength()
	{
		return Max3[0]-Min3[0];
	}
	double GetYLength()
	{
		return Max3[1]-Min3[1];
	}
	double GetZLength()
	{
		return Max3[2]-Min3[2];
	}
};
struct Vector
{
	double X;
	double Y;
	double Z;
	Vector():X(0.0f),Y(0.0f),Z(0.0f){}
	Vector(double x, double y, double z):X(x),Y(y),Z(z){};
};
struct Point3d
{
public:
	double X;
	double Y;
	double Z;
	Point3d():X(0.0f),Y(0.0f),Z(0.0f){}
	Point3d(double x, double y, double z):X(x),Y(y),Z(z){}
	bool operator==(const Point3d &right) const
	{
		if(X == right.X && Y == right.Y && Z == right.Z)
			return true;
		return false;
	}
	bool operator<(const Point3d &right) const
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
struct Triangle
{
public :
	MESH_INDEXTYPE P0Index;
	MESH_INDEXTYPE P1Index;
	MESH_INDEXTYPE P2Index;
	Triangle(MESH_INDEXTYPE p0index, MESH_INDEXTYPE p1index, MESH_INDEXTYPE p2index):P0Index(p0index),P1Index(p1index),P2Index(p2index){}
	Triangle():P0Index(-1),P1Index(-1),P2Index(-1){};
};
struct Tetrahedra
{
	int PIndex[4];
	Tetrahedra(int p0,int p1,int p2,int p3)
	{
		PIndex[0]=p0;
		PIndex[1]=p1;
		PIndex[2]=p2;
		PIndex[3]=p3;
	}
	Tetrahedra()
	{
		PIndex[0]=-1;
		PIndex[1]=-1;
		PIndex[2]=-1;
		PIndex[3]=-1;
	}
};
struct Edge
{
public :
	int PIndex[2];
	Edge(int p0,int p1)
	{
		PIndex[0]=p0;
		PIndex[1]=p1;
	}
	Edge()
	{

	}
};
struct IntIndexTriple
{
	int X;
	int Y;
	int Z;
	int Index;
	IntIndexTriple(int x,int y,int z):X(x),Y(y),Z(z){}
	void Sort()
	{
		int max,min;
		if(X>Y)
		{
			max=X;
			min=Y;
		}
		else
		{
			max=Y;
			min=X;
		}
		if(Z>max)
		{
			X=min;
			Y=max;
		}
		else if(Z<min)
		{
			X=Z;
			Y=min;
			Z=max;
		}
		else
		{
			X=min;
			Y=Z;
			Z=max;
		}
	}
	bool operator<(const IntIndexTriple &right) const
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
	bool operator==(const IntIndexTriple& right) const
	{
		return X==right.X&&Y==right.Y&&Z==right.Z;
	}
};
struct IntIndexDouble
{
	int X;
	int Y;
	int Index;
	IntIndexDouble(int x,int y):X(x),Y(y){}
	void Sort()
	{
		if(X>Y)
		{
			int temp=X;
			X=Y;
			Y=temp;
		}
	}
	bool operator<(const IntIndexDouble &right) const
	{
		if(X < right.X)
			return true;
		else if(X>right.X)
			return false;
		else if(Y < right.Y)
			return true;
		else if(Y >right.Y)
			return false;
		return false;
	}
	bool operator==(const IntIndexDouble& right) const
	{
		return X==right.X&&Y==right.Y;
	}
};
#endif
