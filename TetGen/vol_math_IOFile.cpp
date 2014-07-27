#include "vol_math_IOFile.h"
void PlyManager::WriteXYZFile(std::vector<Int16Triple>& list, const char* fileName)
{
	FILE * nfile = fopen(fileName,"wb");
	fprintf(nfile,"%d\n",list.size());
	for(size_t i=0;i<list.size();i++)
	{
		fprintf(nfile,"%d %d %d\n",list[i].X,list[i].Y,list[i].Z);
	}
	fclose(nfile);
}
void PlyManager::AWriteV(FILE* file, double v1, double v2, double v3,unsigned char r,unsigned char g,unsigned char b)
{
	fprintf(file,"%lf %lf %lf %d %d %d\n",v1,v2,v3,r,g,b);
}
void PlyManager::AWriteF(FILE* file, MESH_INDEXTYPE i1, MESH_INDEXTYPE i2, MESH_INDEXTYPE i3)
{
	fprintf(file,"%d %d %d %d\n",3,i1,i2,i3);
}
void PlyManager::AWriteE(FILE* file,MESH_INDEXTYPE i1,MESH_INDEXTYPE i2)
{
	fprintf(file,"%d %d\n",i1,i2);
}
void PlyManager::ReadFile(Mesh& mesh,const char* fileName)
{
	MESH_INDEXTYPE vcount=0;
	MESH_INDEXTYPE fcount=0;
	FILE * nfile = fopen(fileName,"r");
	fscanf(nfile,"ply\nformat ascii 1.0\ncomment VCGLIB generated\nelement vertex %d\n",&vcount);
	fscanf(nfile,"property double x\nproperty double y\nproperty double z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\nelement face %d\n",&fcount);
	fscanf(nfile,"property list uchar int vertex_indices\nend_header\n");
	double v1=0,v2=0,v3=0;
	int r=0,g=0,b=0,alpha;
	MESH_INDEXTYPE i1=0,i2=0,i3=0;
	for(MESH_INDEXTYPE i=0;i<vcount;i++)
	{
		fscanf(nfile,"%lf %lf %lf %d %d %d %d\n",&v1,&v2,&v3,&r,&g,&b,&alpha);
		//qym
		//mesh.AddVertex(Point3d(v1,v2,v3));
		Point3d p3d(v1,v2,v3);
		mesh.AddVertex(p3d);
	}
	for(MESH_INDEXTYPE j=0;j<fcount;j++)
	{
		fscanf(nfile,"3 %d %d %d\n",&i1,&i2,&i3);
		Triangle t(i1,i2,i3);
		mesh.AddFace(t);
	}
	fclose(nfile);
}
void PlyManager::ReadFileEx(Mesh& mesh,const char* fileName)
{
	MESH_INDEXTYPE vcount=0;
	MESH_INDEXTYPE fcount=0;
	FILE * nfile = fopen(fileName,"r");
	fscanf(nfile,"ply\nformat ascii 1.0\ncomment VCGLIB generated\nelement vertex %d\n",&vcount);
	fscanf(nfile,"property double x\nproperty double y\nproperty double z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face %d\n",&fcount);
	fscanf(nfile,"property list int int vertex_indices\nend_header\n");
	double v1=0,v2=0,v3=0;
	int r=0,g=0,b=0;
	MESH_INDEXTYPE i1=0,i2=0,i3=0;
	for(MESH_INDEXTYPE i=0;i<vcount;i++)
	{
		fscanf(nfile,"%lf %lf %lf %d %d %d\n",&v1,&v2,&v3,&r,&g,&b);
		//qym
		//mesh.AddVertex(Point3d(v1,v2,v3));
		Point3d p3d(v1,v2,v3);
		mesh.AddVertex(p3d);
	}
	for(MESH_INDEXTYPE j=0;j<fcount;j++)
	{
		fscanf(nfile,"3 %d %d %d\n",&i1,&i2,&i3);
		Triangle t(i1,i2,i3);
		mesh.AddFace(t);
	}
	fclose(nfile);
}
void PlyManager::Output(Mesh& mesh,const char* filename)
{
	FILE * nfile = fopen(filename,"wb");
	fprintf(nfile,"ply\n");
	fprintf(nfile,"format ascii 1.0\n");
	fprintf(nfile,"comment VCGLIB generated\n");
	fprintf(nfile,"element vertex %d\n",mesh.Vertices.size());
	fprintf(nfile,"property double x\n");
	fprintf(nfile,"property double y\n");
	fprintf(nfile,"property double z\n");
	fprintf(nfile,"property uchar red\n");
	fprintf(nfile,"property uchar green\n");
	fprintf(nfile,"property uchar blue\n");
	fprintf(nfile,"element face %d\n",mesh.Faces.size());
	fprintf(nfile,"property list int int vertex_indices\n");
	fprintf(nfile,"end_header\n");
	for(size_t i=0;i<mesh.Vertices.size();i++)
	{
		AWriteV(nfile,mesh.Vertices[i].X,mesh.Vertices[i].Y,mesh.Vertices[i].Z,255,255,255);
	}
	for(size_t i=0;i<mesh.Faces.size();i++)
	{
		AWriteF(nfile,mesh.Faces[i].P0Index,mesh.Faces[i].P1Index,mesh.Faces[i].P2Index);
	}
	fclose(nfile);
}
void PlyManager::Output_C(Mesh& mesh,std::vector<Color> & colors,const char* filename)
{
	FILE * nfile = fopen(filename,"wb");
	fprintf(nfile,"ply\n");
	fprintf(nfile,"format ascii 1.0\n");
	fprintf(nfile,"comment VCGLIB generated\n");
	fprintf(nfile,"element vertex %d\n",mesh.Vertices.size());
	fprintf(nfile,"property double x\n");
	fprintf(nfile,"property double y\n");
	fprintf(nfile,"property double z\n");
	fprintf(nfile,"property uchar red\n");
	fprintf(nfile,"property uchar green\n");
	fprintf(nfile,"property uchar blue\n");
	fprintf(nfile,"element face %d\n",mesh.Faces.size());
	fprintf(nfile,"property list int int vertex_indices\n");
	fprintf(nfile,"end_header\n");
	for(size_t i=0;i<mesh.Vertices.size();i++)
	{
		AWriteV(nfile,mesh.Vertices[i].X,mesh.Vertices[i].Y,mesh.Vertices[i].Z,colors[i].R,colors[i].G,colors[i].B);
	}
	for(size_t i=0;i<mesh.Faces.size();i++)
	{
		AWriteF(nfile,mesh.Faces[i].P0Index,mesh.Faces[i].P1Index,mesh.Faces[i].P2Index);
	}
	fclose(nfile);
}
void PlyManager::Output_C2(TetraMesh& mesh,std::vector<Color> & colors,const char* filename)
{
	mesh.Faces.swap(mesh.InnerTriangles);
	Mesh* m=&mesh;
	Output_C(*m,colors,filename);
}
void PlyManager::Output2(std::vector<Point3d> &points,const char* fileName)
{
	FILE * nfile = fopen(fileName,"wb");
	fprintf(nfile,"ply\n");
	fprintf(nfile,"format ascii 1.0\n");
	fprintf(nfile,"comment VCGLIB generated\n");
	fprintf(nfile,"element vertex %d\n",points.size());
	fprintf(nfile,"property double x\n");
	fprintf(nfile,"property double y\n");
	fprintf(nfile,"property double z\n");
	fprintf(nfile,"property uchar red\n");
	fprintf(nfile,"property uchar green\n");
	fprintf(nfile,"property uchar blue\n");
	fprintf(nfile,"element edge %d\n",points.size());
	fprintf(nfile,"property int vertex1\n");
	fprintf(nfile,"property int vertex2\n");
	fprintf(nfile,"end_header\n");
	for(size_t i=0;i<points.size();i++)
	{
		AWriteV(nfile,points[i].X,points[i].Y,points[i].Z,255,255,255);
	}
	for(size_t i=0;i<points.size()-1;i++)
	{
		AWriteE(nfile,i,i+1);
	}
	AWriteE(nfile,points.size()-1,0);
	fclose(nfile);
}
void PlyManager::Output3(std::vector<Point3d> &point,std::vector<Vector> &normals,const char* fileName )
{
	if(point.size()!=normals.size())
		return;
	FILE * nfile = fopen(fileName,"wb");
	fprintf(nfile,"%d\n",point.size());
	for(size_t i=0;i<point.size();i++)
	{
		fprintf(nfile,"%.2f %.2f %.2f %.2f %.2f %.2f\n",point[i].X,point[i].Y,point[i].Z,normals[i].X,normals[i].Y,normals[i].Z);
	}
	fclose(nfile);
}
void PlyManager::Output4(TetraMesh& mesh,const char* fileName ,int type)
{
	if(type==0)
	{
		Mesh* m=&mesh;
		Output(*m,fileName);
	}
	if(type==1)
	{
		mesh.Faces.swap(mesh.InnerTriangles);
		Mesh* m=&mesh;
		Output(*m,fileName);
		mesh.Faces.swap(mesh.InnerTriangles);
	}
	if(type==2)
	{

	}
}

