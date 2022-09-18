//////////////////////////////////////////////////////////
// File: materialTriangles.cpp
// Author: Court Cutting, MD
// Date: 2/26/2015
// Purpose: Triangle storage class using only uniquely linked xyz position
//    and uv texture data for vertices.  Normals are not used
//    and texture seams are not allowed thereby making all
//    vertices unique. Triangle listings include not only
//    the 3 vertex indices, but also an integer material identifier.
//    An auxilliary graphics class
//    with its included vertex shaders provides vertex normal
//    doubling and tripling for graphics and user purposes.
//    Other auxilliary classes provide
//    a tissue specific fragment shader to procedurally texture
//    the model for graphics purposes.
//////////////////////////////////////////////////////////

#include <assert.h>
#include <fstream>
#include <algorithm>
#include <exception>
#include <string.h>
#include <array>
#include <sstream>
#include "materialTriangles.h"
#include "math3d.h"
#include "boundingBox.h"
#include "Mat3x3f.h"
#include "Vec3f.h"
#include <exception>

int materialTriangles::readObjFile(const char *fileName)
{ // returned error codes: 0=no error, 1=can't open file, 2=non-triangle primitive,
	// 3=bad 3D vertex line, 4=bad 3D texture line, 5=bad uvw face line, 6=exceeds 0x3fffffff vertex limit
	// Uses "s [smoothingGroup]" separators in front of face groups to separate materials. [smoothingGroup] must start at 1 as 0 means off.
    std::ifstream fin(fileName);
    if(!fin.is_open())
        return 1;
	_xyz.clear();
	_uv.clear();
	_tris.clear();
	std::vector<float> tuv;
	std::string unparsedLine;
	std::vector<std::string> parsedLine;
	int vtx[2],vertexNumber=0;
	int triangleNumber=0;
	matTriangle triNow;
	triNow.material = 0;
	int i,j,k,l;
	std::string str;
	char s[100];
	while(parseNextInputFileLine(&fin,unparsedLine,parsedLine))
	{
		if(parsedLine.empty())
			continue;
		// ignore all other .obj lines except the following for now
		if(parsedLine[0]=="v")
		{
			if(parsedLine.size()!=4)
				return 3;
			for(i=1; i<4; ++i)
				_xyz.push_back((float)atof(parsedLine[i].c_str()));
			++vertexNumber;
		}
		else if (parsedLine[0] == "usemtl")
		{
			if (parsedLine.size() != 2)
				return 3;
			triNow.material = atoi(parsedLine[1].c_str());
		}
		else if (parsedLine[0] == "vt")
		{
			if(parsedLine.size()!=3)
				return 4;
			_uv.push_back((float)atof(parsedLine[1].c_str()));
			_uv.push_back((float)atof(parsedLine[2].c_str()));
		}
		else if(parsedLine[0]=="f")
		{	// always in vertexPosition/vertexTexture format. If vP/vT may skip normal. If vP//vN, texture is skipped.
			// triangles always input last after vertex positions and textures
			if (_uv.empty()) {  // create unique texture for each vertex
				std::cout << "Error reading .obj file: " << fileName << " . No texture coordinates specified.\n";
				return 4;
			}
			int numVerts = (int)parsedLine.size();
			if(numVerts!=4)
				return 2;
			for(i=1; i<4; ++i)
			{
				strncpy(s,parsedLine[i].c_str(),99);
				k=0; l=0;
				for(j=0; j<2; ++j)  // ignore normals
				{
					while(s[l]!='/' && s[l]!='\0')
						++l;
					str.clear();
					while(k!=l)
					{
						str.push_back(s[k]);
						++k;
					}
					vtx[j] = atoi(str.c_str()) - 1;	// remember indexes in obj files start at 1
					++k; ++l;
				}
				triNow.v[i - 1] = vtx[0];
				triNow.tex[i - 1] = vtx[1];
			}
			_tris.push_back(triNow);
			++triangleNumber;
		}
		else
			continue;
	}
	fin.close();
	if(	vertexNumber>0x3fffffff) {
		return 6;
	}
	// next line minimizes number of program switches on graphics card.
	// only done on startup as later triangle indices must remain unique for incision processing
	partitionTriangleMaterials();
	// trim excess capacity?  Maybe not.  Only going to grow requiring realloc
//	_tris.shrink_to_fit();
//	_xyz.shrink_to_fit();
//	_uv.shrink_to_fit();
	_adjacenciesComputed = false;
	return 0;
}

bool materialTriangles::parseNextInputFileLine(std::ifstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine)
{
	if(infile->eof())
	{
		infile->close();
		return false;
	}
	if(!infile->is_open())
		return false;
	char s[400];
	infile->getline(s,399);
	unparsedLine.assign(s);
	if(unparsedLine=="")
	{
		parsedLine.clear();
		return true;
	}
	parsedLine.clear();
	std::string substr;
	int start=0,next=0;
	while(next<399 && s[next]!= '\0')
	{
		// advance to whitespace
		while(next<399 && s[next]!='\0' && s[next]!=' ' && s[next]!='\t')
			++next;
		substr.assign(s+start,s+next);
		parsedLine.push_back(substr);
		while(s[next]==' ' || s[next]=='\t')
			++next;
		start = next;
	}
	return true;
}

bool materialTriangles::writeObjFile(const char *fileName, const char* materialFileName)
{
	std::string title(fileName);
	if(title.rfind(".obj")>title.size())
		title.append(".obj");
	std::ofstream fout(title.c_str());
    if(!fout.is_open())
        return false;
	char s[400];
	std::string line;
	if (materialFileName != nullptr) {
		line = "mtllib ";
		line.append(materialFileName);
		if (line.rfind(".mtl") > title.size())
			line.append(".mtl");
		line.append("\n");
		fout << line;
	}
	std::vector<int> vIdx, tIdx;
	vIdx.assign(_xyz.size() / 3, -1);
	tIdx.assign(_uv.size()>>1, -1);
	std::vector<matTriangle> tris;
	tris.reserve(_tris.size());
	for (size_t i = 0; i < _tris.size(); ++i) {
		if (_tris[i].material < 0)  // skip deleted triangles
			continue;
		tris.push_back(_tris[i]);
		for (int j = 0; j < 3; ++j) {
			vIdx[_tris[i].v[j]] = 1;
			tIdx[_tris[i].tex[j]] = 1;
		}
	}
	int vn = 1, tn = 1;  // .obj indices start at 1 , not 0
	for (size_t n = vIdx.size(), i = 0; i < n; ++i) {
		if (vIdx[i] > 0)
			vIdx[i] = vn++;
	}
	for (size_t n = tIdx.size(), i = 0; i < n; ++i) {
		if (tIdx[i] > 0)
			tIdx[i] = tn++;
	}
	for (size_t n=vIdx.size(), i = 0; i < n; ++i) {
		if (vIdx[i] < 0)
			continue;
		float* fp = &_xyz[i * 3];
		sprintf(s, "v %f %f %f\n", fp[0], fp[1], fp[2]);
		line.assign(s);
		fout.write(line.c_str(), line.size());
	}
	for (size_t n=tIdx.size(), i = 0; i < n; ++i) {
		if (tIdx[i] < 0)
			continue;
		float* fp = &_uv[i<<1];
		sprintf(s, "vt %f %f\n", fp[0], fp[1]);
		line.assign(s);
		fout.write(line.c_str(), line.size());
	}
	std::sort(tris.begin(), tris.end(),
		[](const matTriangle& a, const matTriangle& b) -> bool
		{
			return a.material < b.material;
		});
	int smoothNum=1, lastMaterial = -1;
	for(size_t n=tris.size(), i=0; i<n; ++i)	{
		if (tris[i].material > lastMaterial) {  // new material/shading group
			lastMaterial = tris[i].material;
			sprintf(s, "usemtl %d\n", lastMaterial);
			line.assign(s);
			fout.write(line.c_str(), line.size());
			sprintf(s, "s %d\n", smoothNum++);
			line.assign(s);
			fout.write(line.c_str(), line.size());
		}
		sprintf(s, "f %d/%d %d/%d %d/%d\n", vIdx[tris[i].v[0]], tIdx[tris[i].tex[0]], vIdx[tris[i].v[1]], tIdx[tris[i].tex[1]], vIdx[tris[i].v[2]], tIdx[tris[i].tex[2]]);
		line.assign(s);
		fout.write(line.c_str(),line.size());	}
	fout.close();
	return true;
}

void materialTriangles::getVertexCoordinate(unsigned int vertex, float(&xyz)[3]) const
{	// type safe version
	const float *v = &_xyz[vertex*3];
	xyz[0]=v[0]; xyz[1]=v[1]; xyz[2]=v[2];
}

bool materialTriangles::getBarycentricProjection(const int triangle, const float(&xyz)[3], float(&uv)[2])
{	// for position xyz return barycentric uv projection into triangle
	float *p,*q;
	int *t = &_tris[triangle].v[0];
	p = vertexCoordinate(t[0]);
	q = vertexCoordinate(t[1]);
	Vec3f u,v,xmp;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	q = vertexCoordinate(t[2]);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	xmp.set(xyz[0]-p[0],xyz[1]-p[1],xyz[2]-p[2]);
	float a,b,c,d;
	a=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
	b=u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	c=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	if(fabs(d=b*b-a*c)<1e-16f)	{	// degenerate triangle
		uv[0]=0.0; uv[1]=0.0f;
		return false; }
	uv[1] = ((u*b - v*a)*xmp)/d;
	uv[0] = (xmp*u - uv[1]*b)/a;
	return true;
}

void materialTriangles::getBarycentricTexture(const int triangle, const float (&uv)[2], float (&texture)[2])
{
	int *tr = &_tris[triangle].tex[0];
	float p=1.0f-uv[0]-uv[1],*t0=getTexture(tr[0]),*t1=getTexture(tr[1]),*t2=getTexture(tr[2]);
	for(int i=0; i<2; ++i)
		texture[i] = t0[i]*p + uv[0]*t1[i] + uv[1]*t2[i];
}

void materialTriangles::getBarycentricPosition(const int triangle, const float (&uv)[2], float (&xyz)[3])
{	// for barycentric uv in triangle returns position in xyz
	float *p,*q;
	int *t = &_tris[triangle].v[0];
	p = vertexCoordinate(t[0]);
	q = vertexCoordinate(t[1]);
	Vec3f u,v,r;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	q = vertexCoordinate(t[2]);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	r = u*uv[0] + v*uv[1];
	xyz[0]=r.x()+p[0]; xyz[1]=r.y()+p[1]; xyz[2]=r.z()+p[2];
}

void materialTriangles::getBarycentricNormal(const int triangle, const float(&uv)[2], float(&nrm)[3])
{
	// look for uniform surface w triangle
	int *tr = triangleVertices(triangle);
	Vec3f vNorm[3],bNorm;
	for (int i = 0; i < 3; ++i) 
		getMeanVertexNormal(tr[i], vNorm[i]._v, _tris[triangle].material);
	bNorm = vNorm[1] * uv[0];
	bNorm += vNorm[2] * uv[1];
	bNorm += vNorm[0] * (1.0f - uv[0] - uv[1]);
	bNorm.normalize();
	nrm[0] = bNorm.x(); nrm[1] = bNorm.y(); nrm[2] = bNorm.z();
}

int materialTriangles::rayIntersect(const float *rayStart, const float *rayDirection, std::vector<int> &triangles, std::vector<float> &params)
{ // lineStart and lineDirection are both 3 element vectors
	std::map<float,int> hitMap;
	std::map<float,int>::iterator hit;
	int j,bigAxis=0;
	if(fabs(rayDirection[1])>fabs(rayDirection[0]))
		bigAxis =1;
	if(fabs(rayDirection[2])>fabs(rayDirection[bigAxis]))
		bigAxis =2;
	int *tr;
	unsigned int i,tNum=(unsigned int)_tris.size();
	float *v[3];
	float minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<tNum; ++i)	{
		tr = &_tris[i].v[0];
		if(_tris[i].material<0)
			continue;
		v[0]=vertexCoordinate(tr[0]);
		v[1]=vertexCoordinate(tr[1]);
		v[2]=vertexCoordinate(tr[2]);
		minimax[0]=minimax[1]=v[0][0];
		minimax[2]=minimax[3]=v[0][1];
		minimax[4]=minimax[5]=v[0][2];
		for(j=1; j<3; ++j) {
			if(minimax[0]>v[j][0])
				minimax[0]=v[j][0];
			if(minimax[1]<v[j][0])
				minimax[1]=v[j][0];
			if(minimax[2]>v[j][1])
				minimax[2]=v[j][1];
			if(minimax[3]<v[j][1])
				minimax[3]=v[j][1];
			if(minimax[4]>v[j][2])
				minimax[4]=v[j][2];
			if(minimax[5]<v[j][2])
				minimax[5]=v[j][2];
		}
		tMin = (minimax[bigAxis<<1]-rayStart[bigAxis])/rayDirection[bigAxis];
		tMax = (minimax[(bigAxis<<1)+1]-rayStart[bigAxis])/rayDirection[bigAxis];
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=rayStart[j]+rayDirection[j]*tMin;
			rMax=rayStart[j]+rayDirection[j]*tMax;
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		// now look for triangle intersection
		float b[3],m[9],r[3];
		for(j=0; j<3; ++j) {
			b[j]=rayStart[j]-v[0][j];
			m[(j<<1)+j] = -rayDirection[j];
			m[(j<<1)+j+1] = v[1][j]-v[0][j];
			m[(j<<1)+j+2] = v[2][j]-v[0][j];
		}
		if(!m3dSolveLinearSystem3D(m,b,r))
			continue;
		if(r[1]<0.0f || r[2]<0.0f || r[1]+r[2]>1.0f)
			continue;
		hitMap.insert(std::make_pair(r[0],i));
	}
	triangles.clear();	params.clear();
	triangles.reserve(hitMap.size());
	params.reserve(hitMap.size());
	for(hit=hitMap.begin(); hit!=hitMap.end(); ++hit)	{
		params.push_back(hit->first);
		triangles.push_back(hit->second);
	}
	return (int)hitMap.size();
}

int materialTriangles::findAdjacentTriangles(bool forceCompute, bool fullManifoldTest)
{	// computes all the adjacent triangles from raw triangle input
	// returns false if non-manifold surface is input
	if (_adjacenciesComputed && !(fullManifoldTest || forceCompute))
		return true;
	typedef std::set<edge, edgeTest> edgeSet;
	typedef edgeSet::iterator edgeIt;
	std::pair <edgeIt,bool> P;
	edge E;
	edgeIt ei;
	edgeSet M;
	M.clear();
	int *tnow;
	unsigned int *adjNow;
	unsigned int i,j,tcode,numtris=(unsigned int)_tris.size();
	if (numtris < 1)
		return true;
	_adjs.clear();
	_adjs.assign(numtris*3,0x00000003);
	for(i=0; i<numtris; ++i)
	{
		if(_tris[i].material<0)	// signals a deleted triangle
			continue;
		tnow = &_tris[i].v[0];
		adjNow = &(_adjs[(i<<1)+i]); //_adjs[i*3]
		for(j=0; j<3; j++) {
			if(adjNow[j]!=0x00000003)	// adjacency already computed
				continue;
			int tmp;
			if((tmp=tnow[(j+1)%3])<tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j];
				E.reversed = 1;	}
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp;
				E.reversed = 0;	}
			E.matched = 0;
			E.adjCode = (i<<2) + j; //(4i+j)
			P = M.insert(E);
			// if P.second is true, no match so edge inserted
			if(P.second==false)	// edge match found
			{
				if (fullManifoldTest && P.first->matched) {
					throw(std::logic_error("Failed fullManifoldTest with tripled edge in findAdjacentTraiangles()."));
					return 3;  // third edge with these coordinates so not manifold
				}
				if (P.first->reversed == E.reversed && E.vtxMin != E.vtxMax) {
					throw(std::logic_error("triangle ordering error"));
					return 2;  // triangle ordering error
				}
				tcode = P.first->adjCode;
				adjNow[j] = tcode;
				_adjs[(tcode>>2)*3+(tcode&0x00000003)] = E.adjCode;
				M.erase(P.first);
				if(fullManifoldTest) {  // manifold test for tripled edge
					E.matched = 1;
					M.insert(E); }
			}
			else
				adjNow[j] = 0x00000003;
		}
	}
	makeVertexToTriangleMap();
	_adjacenciesComputed = true;
	if (M.size() > 0)
		return 1;
	else
		return 0;
}

void materialTriangles::makeVertexToTriangleMap()
{
	int i, j, numtris = (int)_tris.size();
	_vertexFace.clear();
	if (_xyz.size() < 1){  // allows processing of only topology
		int maxV = -1;
		for (i = 0; i < numtris; ++i){
			for (j = 0; j < 3; ++j){
				if (_tris[i].v[j]>maxV)
					maxV = _tris[i].v[j];
			}
		}
		assert(maxV > -1);
		_vertexFace.assign(maxV+1, 0x80000000);	// initially deleted
	}
	else
		_vertexFace.assign(_xyz.size() / 3, 0x80000000);	// initially deleted
	unsigned int vnow;
	int *tnow;
	// provide each vertex with a face it is a member of
	for(i=0; i<numtris; ++i)
	{
		tnow = &(_tris[i].v[0]);
		if (_tris[i].material<0)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
		{
			vnow = tnow[j];
			if(_vertexFace[vnow]&0x40000000)
				continue;	// vertex first on free edge, don't change
			_vertexFace[vnow] = i;
			if(_adjs[i*3+j]==0x00000003)	// vertex first on free edge, lock it for easy neighbor find
				_vertexFace[vnow] |= 0x40000000;
		}
	}
}

void materialTriangles::getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors)
{
	unsigned int trNum,adj,triStart;
	triStart = _vertexFace[vertex];
	neighbors.clear();
	if(triStart&0x80000000)	// unconnected vertex
		return;
	trNum = triStart&0x3fffffff;
	unsigned int *adjs;
	int *tnow = &(_tris[trNum].v[0]);
	assert(_tris[trNum].material>-1);	// deleted triangle
	int j;
	for(j=0; j<3; ++j)
		if(tnow[j]==vertex)
			break;
	assert(j<3);
	adjs = &(_adjs[(trNum<<1)+trNum]);
	// set triStart to the end adjacency code for counterclockwise traversal
	neighborNode n;
	if(triStart&0x40000000)	// started on a free edge, will end on one
	{
		triStart = 0x00000003;
		n.triangle = -1;	// code for open ring neighbors
		n.vertex = tnow[(j+1)%3];
		neighbors.push_back(n);
	}
	else	// create adjacency code of starting edge
		triStart = (trNum<<2)+j;
	n.triangle = trNum;
	n.vertex = tnow[(j+2)%3];
	neighbors.push_back(n);
	adj = adjs[(j+2)%3];
	while(adj!=triStart)
	{
		n.triangle = adj>>2;
		tnow = &(_tris[n.triangle].v[0]);
		adjs = &(_adjs[(n.triangle<<1)+n.triangle]);
		j = adj&0x00000003;
		n.vertex = tnow[(j+2)%3];
		neighbors.push_back(n);
		adj = adjs[(j+2)%3];
	}
}

materialTriangles::materialTriangles(const materialTriangles& x)
{
	_tris.assign(x._tris.begin(),x._tris.end());
	_xyz.assign(x._xyz.begin(), x._xyz.end());
	_uv.assign(x._uv.begin(), x._uv.end());
//	_matEnds.insert(x._matEnds.begin(), x._matEnds.end());
	_adjacenciesComputed = x._adjacenciesComputed;
	_adjs.assign(x._adjs.begin(), x._adjs.end());
	_vertexFace.assign(x._vertexFace.begin(), x._vertexFace.end());
	_name = x._name;
}

materialTriangles::materialTriangles(void) :_adjacenciesComputed(false)
{
}


materialTriangles::~materialTriangles(void)
{
}

materialTriangles::matTriangle* materialTriangles::getTriangleArray(int &numberOfTriangles)
{
	numberOfTriangles = (int)_tris.size();
	if (numberOfTriangles < 1)
		return nullptr;
	return &_tris[0];
}

float* materialTriangles::getPositionArray(int &numberOfVertices)
{
	numberOfVertices = (int)_xyz.size()/3;
	if (numberOfVertices < 1)
		return nullptr;
	return &_xyz[0];
}

float* materialTriangles::getTextureArray(int &numberOfVertices)
{
	numberOfVertices = (int)_uv.size()>>1;
	if (numberOfVertices < 1)
		return nullptr;
	return &_uv[0];
}

/*bool materialTriangles::localPick(const float *lineStart, const float *lineDirection, float(&position)[3], int &triangle, float &param, const int onlyMaterial)
{ // lineStart and lineDirection are both 3 element vectors
	std::map<float, lineHit> hits;
	rayHits(lineStart, lineDirection, hits);
	for (auto it = hits.begin(); it != hits.end(); ++it) {
		if (it->first < -1e-8f)
			continue;
		if (it->first > 1e-6f || (onlyMaterial > -1 && _tris[it->second.triangle].material != onlyMaterial))
			return false;
		triangle = it->second.triangle;
		param = it->first;
		position[0] = it->second.v._v[0];
		position[1] = it->second.v._v[1];
		position[2] = it->second.v._v[2];
	}
	return false;


	bool picked=false;
	triangle = -1;
	param = 1e30f;
	int j,bigAxis=0;
	if(fabs(lineDirection[1])>fabs(lineDirection[0]))
		bigAxis =1;
	if(fabs(lineDirection[2])>fabs(lineDirection[bigAxis]))
		bigAxis =2;
	int *tr,i;
	float *v[3];
	float t,vtx[3],minimax[6],tMin,tMax,rMin,rMax;
	bool matTest = true;
	if (onlyMaterial<0)
		matTest = false;
	for (i = 0; i<(int)_tris.size(); ++i)	{
		if (matTest && _tris[i].material != onlyMaterial)
			continue;
		tr = &_tris[i].v[0];
		if (_tris[i].material<0)
			continue;
		v[0]=vertexCoordinate(tr[0]);
		v[1]=vertexCoordinate(tr[1]);
		v[2]=vertexCoordinate(tr[2]);
		minimax[0]=minimax[1]=v[0][0];
		minimax[2]=minimax[3]=v[0][1];
		minimax[4]=minimax[5]=v[0][2];
		for(j=1; j<3; ++j) {
			if(minimax[0]>v[j][0])
				minimax[0]=v[j][0];
			if(minimax[1]<v[j][0])
				minimax[1]=v[j][0];
			if(minimax[2]>v[j][1])
				minimax[2]=v[j][1];
			if(minimax[3]<v[j][1])
				minimax[3]=v[j][1];
			if(minimax[4]>v[j][2])
				minimax[4]=v[j][2];
			if(minimax[5]<v[j][2])
				minimax[5]=v[j][2];
		}
		tMin = (minimax[bigAxis<<1]-lineStart[bigAxis])/lineDirection[bigAxis];
		tMax = (minimax[(bigAxis<<1)+1]-lineStart[bigAxis])/lineDirection[bigAxis];
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=lineStart[j]+lineDirection[j]*tMin;
			rMax=lineStart[j]+lineDirection[j]*tMax;
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		if(m3dRayTriangleIntersection(lineStart,lineDirection,v[0],v[1],v[2],t,vtx)) {
			if(fabs(t)<param){
				param = fabs(t);
				picked=true;
				triangle = i;
				position[0]=vtx[0]; position[1]=vtx[1]; position[2]=vtx[2];
			}
		}
	}
	return picked;
} */

int materialTriangles::linePick(const float *lineStart, const float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params, const int onlyMaterial)
{
	std::map<float, lineHit> hits;
	rayHits(lineStart, lineDirection, hits);
	positions.clear();	triangles.clear();	params.clear();
	positions.reserve(hits.size() * 3);
	triangles.reserve(hits.size());
	params.reserve(hits.size());
	for (auto it = hits.begin(); it != hits.end(); ++it) {
		if (onlyMaterial>-1 && _tris[it->second.triangle].material != onlyMaterial)
			continue;
		params.push_back(it->first);
		triangles.push_back(it->second.triangle);
		positions.push_back(it->second.v[0]);
		positions.push_back(it->second.v[1]);
		positions.push_back(it->second.v[2]);
	}
	return (int)params.size();
}

int materialTriangles::linePick(const float *lineStart, const float *lineDirection, std::vector<float> &rayParams, std::vector<int> &triangles, std::vector<Vec2f> &triangleParams)
{
	std::map<float, lineHit> hits;
	rayHits(lineStart, lineDirection, hits);
	rayParams.clear();   triangles.clear();	triangleParams.clear();
	rayParams.reserve(hits.size() * 3);
	triangles.reserve(hits.size());
	triangleParams.reserve(hits.size());
	for (auto it = hits.begin(); it != hits.end(); ++it) {
		rayParams.push_back(it->first);
		triangles.push_back(it->second.triangle);
		triangleParams.push_back(it->second.uv);
	}
	return (int)rayParams.size();
}

bool materialTriangles::localPick(const float *lineStart, const float *lineDirection, float(&position)[3], int &triangle, float (&triangleParam)[2], const int onlyMaterial)
{ // lineStart and lineDirection are both 3 element vectors
	std::map<float, lineHit> hits;
	rayHits(lineStart, lineDirection, hits);
	for (auto it = hits.begin(); it != hits.end(); ++it) {
		if (it->first < -1e-8f)
			continue;
		if (onlyMaterial > -1 && _tris[it->second.triangle].material != onlyMaterial)
			continue;;
		triangle = it->second.triangle;
		triangleParam[0] = it->second.uv._v[0];
		triangleParam[1] = it->second.uv._v[1];
		position[0] = it->second.v._v[0];
		position[1] = it->second.v._v[1];
		position[2] = it->second.v._v[2];
		return true;
	}
	return false;
}

int materialTriangles::rayHits(const float *rayStart, const float *rayDirection, std::map<float, lineHit> &hits)
{ // lineStart and lineDirection are both 3 element vectors
	findAdjacentTriangles();  // may have already been done so don't force it
	Vec3f lS(rayStart[0], rayStart[1], rayStart[2]),lD(rayDirection[0], rayDirection[1], rayDirection[2]);
	hits.clear();
	lineHit pT;
	std::map<float,lineHit>::iterator hit,hit2;
	int j,bigAxis=0;
	if(fabs(rayDirection[1])>fabs(rayDirection[0]))
		bigAxis =1;
	if(fabs(rayDirection[2])>fabs(rayDirection[bigAxis]))
		bigAxis =2;
	int *tr,i;
	float *v[3];
	float t,minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<(int)_tris.size(); ++i)	{
		tr = &_tris[i].v[0];
		if (_tris[i].material<0)
			continue;
		v[0]=vertexCoordinate(tr[0]);
		v[1]=vertexCoordinate(tr[1]);
		v[2]=vertexCoordinate(tr[2]);
		minimax[0]=minimax[1]=v[0][0];
		minimax[2]=minimax[3]=v[0][1];
		minimax[4]=minimax[5]=v[0][2];
		for(j=1; j<3; ++j) {
			if(minimax[0]>v[j][0])
				minimax[0]=v[j][0];
			if(minimax[1]<v[j][0])
				minimax[1]=v[j][0];
			if(minimax[2]>v[j][1])
				minimax[2]=v[j][1];
			if(minimax[3]<v[j][1])
				minimax[3]=v[j][1];
			if(minimax[4]>v[j][2])
				minimax[4]=v[j][2];
			if(minimax[5]<v[j][2])
				minimax[5]=v[j][2];
		}
		tMin = (minimax[bigAxis<<1]-rayStart[bigAxis])/rayDirection[bigAxis];
		tMax = (minimax[(bigAxis<<1)+1]-rayStart[bigAxis])/rayDirection[bigAxis];
		float pad = (tMax-tMin)*0.1f; // fix roundoff error problem
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=rayStart[j]+rayDirection[j]*(tMin - pad);
			rMax=rayStart[j]+rayDirection[j]*(tMax + pad);
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		if(rayTriangleIntersection(lS, lD, i, t, pT.uv._v, pT.v)) {
			pT.triangle = i;
			hits.insert(std::make_pair(t,pT));
		}
	}
	std::set<int> neiSet;
	auto addVertexNeighbors = [&](int vert) {
		std::vector<materialTriangles::neighborNode> nei;
		getNeighbors(vert, nei);
		auto nit = nei.begin();
		if (nit->triangle < 0)
			++nit;
		while (nit != nei.end()){
			neiSet.insert(nit->triangle);
			++nit;
		}
	};
	// remove duplicated triangle edge and vertex repeats
	for(hit=hits.begin(); hit!=hits.end(); ++hit)	{
		hit2 = hit;
		++hit2;
		while (hit2 != hits.end()) {
			if (hit2->first - hit->first < 1e-4f) {
				unsigned int atr;
				neiSet.clear();
				if (hit2->second.uv[0] < 1e-5f){
					if (hit2->second.uv[1] < 1e-5f)
						addVertexNeighbors(_tris[hit2->second.triangle].v[0]);
					else if (hit2->second.uv[1] > 0.9999f)
						addVertexNeighbors(_tris[hit2->second.triangle].v[2]);
					else {
						if ((atr = _adjs[hit2->second.triangle * 3 + 2]) != 3)
							neiSet.insert(atr >> 2);
					}
				}
				else if (hit2->second.uv[1] < 1e-5f){
					if (hit2->second.uv[0] > 0.9999f)
						addVertexNeighbors(_tris[hit2->second.triangle].v[1]);
					else {
						if ((atr = _adjs[hit2->second.triangle * 3]) != 3)
							neiSet.insert(atr >> 2);
					}
				}
				else if (hit2->second.uv[0] + hit2->second.uv[1] > 0.9999f){
					if((atr = _adjs[hit2->second.triangle * 3 + 1]) != 3)
						neiSet.insert(atr >> 2);
				}
				else
					;
				if (!neiSet.empty() && neiSet.find(hit->second.triangle) != neiSet.end())
					hit2 = hits.erase(hit2);
				else
					++hit2;
			}
			else
				break;
		}
	}
	return (int)hits.size();
}

bool materialTriangles::rayTriangleIntersection(const Vec3f &rayOrigin, const Vec3f &rayDirection, const int triangle, float &rayParam, float(&triParam)[2], Vec3f &intersect)
{
	Vec3f b, r, U, V, t[3];
	int *tr = triangleVertices(triangle);
	for (int i = 0; i < 3; ++i)
		getVertexCoordinate(tr[i], t[i]._v);
	b = rayOrigin - t[0];
	U = t[1] - t[0];
	V = t[2] - t[0];
	Mat3x3f m(-rayDirection, U, V);
	r = m.Robust_Solve_Linear_System(b);
	if (r._v[1]<-1e-4f || r._v[2]<-1e-4f || r._v[1] >1.0001f || r._v[2] >1.0001f || r._v[1] + r._v[2]>1.0001f) // allows for roundoff error
		return false;
	rayParam = r._v[0];
	triParam[0] = r._v[1];
	triParam[1] = r._v[2];
	intersect = t[0] + U*triParam[0] + V*triParam[1];
	return true;
}

void materialTriangles::getTriangleNormal(int triangle, float (&normal)[3], bool normalized)
{
	int *tr = &_tris[triangle].v[0];
	Vec3f v0,v1,t0((float(&)[3])_xyz[(tr[0]<<1)+tr[0]]),t1((float(&)[3])_xyz[(tr[1]<<1)+tr[1]]),t2((float(&)[3])_xyz[(tr[2]<<1)+tr[2]]);
	v0 = t1-t0;
	v1 = t2-t0;
	t0 = v0^v1;
	if(normalized)
		t0.normalize();
	normal[0]=t0.x(); normal[1]=t0.y(); normal[2]=t0.z();
}

void materialTriangles::getAreaNormal(const int triangle, const float(&uv)[2], const float radius, float(&normal)[3], bool normalized)
{
	Vec3f N, T, P;
	getBarycentricPosition(triangle, uv, P._v);
	std::set<int> trisDone;
	recurseTriangleNormals(triangle, trisDone, P._v, radius*radius, N._v);
	if (normalized)
		N.normalize();
	normal[0] = N._v[0]; normal[1] = N._v[1]; normal[2] = N._v[2];
}

void materialTriangles::recurseTriangleNormals(const int triangle, std::set<int> &trisDone, float (&center)[3], float radiusSq, float(&normalSum)[3])
{
	if (!trisDone.insert(triangle).second)
		return;
	Vec3f N;
	getTriangleNormal(triangle, N._v, false);
	normalSum[0] += N._v[0]; normalSum[1] += N._v[1]; normalSum[2] += N._v[2];
	for (int i = 0; i < 3; ++i) {
		unsigned int adj = _adjs[triangle * 3 + i];
		if (adj == 3)
			continue;
		getVertexCoordinate(_tris[adj>>2].v[((adj&3)+2)&3], N._v);
		if ((N - Vec3f(center)).length2() < radiusSq)
			recurseTriangleNormals(adj >> 2, trisDone, center, radiusSq, normalSum);
	}
}

void materialTriangles::getNearestHardEdge(float(&xyz)[3], int &triangle, int &edge, float &param, int materialLimit)
{	// Input xyz, then overwrites all 4 with the point on the nearest hard edge.
	// If materialLimit>-1 searches only triangle edges whose material==materialLimit.  If materialLimit==-1 searches all.
	if(!_adjacenciesComputed)
		findAdjacentTriangles();
	int *tr,i,j,n=(int)_tris.size();
	Vec3f pt(xyz), e0, e1;
	float minDsq=1e30f;
	unsigned int *adj;
	for(i=0; i<n; ++i)	{
		tr = &_tris[i].v[0];
		if (_tris[i].material<0 || (materialLimit>-1 && materialLimit != _tris[i].material))
			continue;
		adj = &_adjs[i*3];
		for (j = 0; j < 3; ++j) {
			if (adj[j]!=3 && _tris[i].material == _tris[adj[j] >> 2].material)  // not a hard edge
				continue;
			// candidate edge
			e0.set((float(&)[3])_xyz[tr[j] * 3]);
			e1.set((float(&)[3])_xyz[tr[(j + 1) % 3] * 3]);
			Vec3f v0 = pt - e0, v1 = pt - e1;
			float d = v0*v0, t;
			if (d > minDsq || v1.length2() > minDsq)
				continue;
			e1 -= e0;
			assert(e1*e1 > 0.0f);
			t = (v0*e1) / (e1*e1);
			if (t < 0.0f){
				if (minDsq > d) {
					minDsq = d;
					xyz[0] = e0._v[0];
					xyz[1] = e0._v[1];
					xyz[2] = e0._v[2];
					triangle = i;
					edge = j;
					param = 0.0;
				}
			}
			else if (t > 1.0f){
				if (minDsq > (d = v1*v1)) {
					minDsq = d;
					xyz[0] = (e0 + e1)._v[0];
					xyz[1] = (e0 + e1)._v[1];
					xyz[2] = (e0 + e1)._v[2];
					triangle = i;
					edge = j;
					param = 1.0;
				}
			}
			else {
				v1 = e0 + e1*t;
				v0 = v1 - pt;
				if ((d = v0*v0) < minDsq) {
					minDsq = d;
					xyz[0] = v1._v[0];
					xyz[1] = v1._v[1];
					xyz[2] = v1._v[2];
					triangle = i;
					edge = j;
					param = t;
				}
			}
		}
	}
}

void materialTriangles::interpolateEdgeTextures(int triangle, int edge, int newVert, float param)
{	// assumes triangle hasn't been changed yet by newVert

	assert(false);

	int *trVerts = &_tris[triangle].tex[0];
	float *txOut = getTexture(newVert),*txIn;
	txIn=getTexture(trVerts[edge]);
	txOut[0]=(1.0f-param)*txIn[0]; txOut[1]=(1.0f-param)*txIn[1];
	txIn=getTexture(trVerts[(edge+1)%3]);
	txOut[0]+=param*txIn[0]; txOut[1]+=param*txIn[1];
}

int materialTriangles::splitTriangleEdge(int triangle, int edge, const float parameter)
{	// Splits a triangle along edge(0-2) by parameter(0-1).
	// Creates 1 or 2 new triangles and one new vertex. Returns new vertex number of the input triangle split.
	// Does fix adjacency array so getNeighbors() works and texture interpolated.
	// Definitely should call getTriangleAdjacencies() ASAP.

	if (triangle == 26395)
		int junk = 0;

	assert(0.0f<=parameter && 1.0f>=parameter);
	int *trVerts = &_tris[triangle].v[0], *trTex = &_tris[triangle].tex[0];
	if (_tris[triangle].material<0)
		return -1;
	int tn,newVert= addVertices(1);
	int ve=trVerts[edge],ve1=trVerts[(edge+1)%3];
	float gv[3], tx[2];
	float *gvp=vertexCoordinate(ve), *txp=getTexture(trTex[edge]);
	gv[0]=(1.0f-parameter)*gvp[0]; gv[1]=(1.0f-parameter)*gvp[1]; gv[2]=(1.0f-parameter)*gvp[2];
	tx[0] = (1.0f - parameter) * txp[0]; tx[1] = (1.0f - parameter) * txp[1];
	gvp = vertexCoordinate(ve1);
	gv[0]+=parameter*gvp[0]; gv[1]+=parameter*gvp[1]; gv[2]+=parameter*gvp[2];
	txp = getTexture(trTex[(edge + 1) % 3]);
	tx[0] += parameter * txp[0]; tx[1] += parameter * txp[1];
	setVertexCoordinate(newVert,gv);
	int tx0 = addTexture();
	setTexture(tx0, tx);
	int v[3], tex[3];
	v[1] = trVerts[(edge + 1) % 3];
	tex[1] = trTex[(edge + 1) % 3];
	trVerts[(edge + 1) % 3] = newVert;
	trTex[(edge + 1) % 3] = tx0;
	v[0] = newVert;	v[2] = trVerts[(edge + 2) % 3];
	tex[0] = tx0;	tex[2] = trTex[(edge + 2) % 3];
	tn = addTriangle(v, _tris[triangle].material, tex);	// invalidates old _tris and _adjs pointers
	if(_adjs[triangle * 3 + edge] == 0x00000003)	{
		_adjs[tn*3] = 0x00000003;
		unsigned int adjTE = _adjs[triangle*3+((edge+1)%3)];
		_adjs[tn*3+1] = adjTE;
		if(adjTE != 3)
			_adjs[(adjTE>>2)*3 + (adjTE&3)] = (tn<<2)+1;
		_adjs[tn*3+2] = (triangle<<2)+((edge+1)%3);
		_adjs[triangle*3+((edge+1)%3)] = (tn<<2)+2;
		if((_vertexFace[v[1]]&0x3fffffff)==triangle)
			_vertexFace[v[1]] = (tn | 0x40000000);	// first vertex on a free edge
		_vertexFace[newVert] = (tn | 0x40000000);
		return newVert;
	}
	int tx1 = tx0;
	unsigned int* trAdjs = &(_adjs[triangle * 3]);
	int *trVertsA = &_tris[trAdjs[edge] >> 2].v[0], *trTexA = &_tris[trAdjs[edge]>>2].tex[0];
	unsigned int* trAdjsA = &(_adjs[(trAdjs[edge] >> 2) * 3]);
	int ea = trAdjs[edge] & 0x00000003, ta = trAdjs[edge] >> 2;
	if (trTexA[ea] != tex[1] || trTexA[(ea + 1) % 3] != trTex[edge]) {  // texture seam at edge
		txp = getTexture(trTexA[ea]);
		tx[0] = parameter * txp[0]; tx[1] = parameter * txp[1];
		txp = getTexture(trTexA[(ea + 1) % 3]);
		tx[0] += (1.0f - parameter) * txp[0]; tx[1] += (1.0f - parameter) * txp[1];
		tx1 = addTexture();
		setTexture(tx1, tx);
		if (_tris[triangle].material == _tris[ta].material) {
			int twoTx[2] = { tx0, tx1 };
			addOneMaterialTextureSeamVertex(newVert, twoTx);
		}
	}
	v[1] = trVertsA[(ea + 1) % 3];
	tex[1] = trTexA[(ea + 1) % 3];
	trVertsA[(ea + 1) % 3] = newVert;
	trTexA[(ea + 1) % 3] = tx1;
	v[0] = newVert;	v[2] = trVertsA[(ea + 2) % 3];
	tex[0] = tx1;	tex[2] = trTexA[(ea + 2) % 3];
	int tna = addTriangle(v, _tris[ta].material, tex);	// invalidates old _tris and _adjs pointers
	trAdjs = &(_adjs[triangle * 3]);
	trAdjsA = &(_adjs[(trAdjs[edge] >> 2) * 3]);
	// new adj assignments
	unsigned int ae1,aa1;
	trAdjs[edge] = (tna<<2);
	trAdjsA[ea] = (tn << 2);
	ae1 = trAdjs[(edge+1)%3];
	aa1 = trAdjsA[(ea + 1) % 3];
	trAdjs[(edge+1)%3] = (tn<<2) + 2;
	trAdjsA[(ea + 1) % 3] = (tna << 2) + 2;
	if(ae1 != 3)
		_adjs[(ae1>>2)*3 + (ae1&3)] = (tn << 2) + 1;
	if (aa1 != 3)
		_adjs[(aa1 >> 2) * 3 + (aa1 & 3)] = (tna << 2) + 1;
	// now set adjacencies for tn and tna
	trAdjs = &(_adjs[tn * 3]);
	trAdjsA = &(_adjs[tna * 3]);
	trAdjs[0] = (ta << 2) + ea;
	trAdjsA[0] = (triangle << 2) + edge;
	trAdjs[1] = ae1;
	trAdjsA[1] = aa1;
	trAdjs[2] = (triangle << 2) + ((edge + 1) % 3);
	trAdjsA[2] = (ta << 2) + ((ea + 1) % 3);
	// new vertexFace assignments
	_vertexFace[newVert] = triangle;
	if((_vertexFace[ve1]&0x3fffffff)==triangle) {
		if((_vertexFace[ve1]&0x40000000)>0)
			_vertexFace[ve1] = tn|0x40000000;
		else
			_vertexFace[ve1] = tn;
	}
	if((_vertexFace[ve]&0x3fffffff)==ta) { // va1
		if((_vertexFace[ve]&0x40000000)>0)
			_vertexFace[ve] = tna|0x40000000;
		else
			_vertexFace[ve] = tna;
	}
	return newVert;
}

int materialTriangles::addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2])
{	// creates 2 new triangles and one new vertex. Returns new vertex position number.
	// Input uvParameters[2] are parameters along vectors t1-t0 and t2-t0.
	// Does fix adjacency array so getNeighbors() works and texture & positions interpolated.
	// Definitely should getTriangleAdjacencies() ASAP.
	assert(uvParameters[0]>=0.0f && uvParameters[0]<=1.0f);
	assert(uvParameters[1]>=0.0f && uvParameters[1]<=1.0f);
	assert(uvParameters[0]+uvParameters[1]<=1.0001f);
	if (_tris[triangle].material < 0) {
		throw(std::logic_error("Trying to add a vertex into a deleted triangle."));
		return -1;
	}
	int* trVerts = triangleVertices(triangle);
	if(uvParameters[0]<0.0002f && uvParameters[1]<0.0002f)  // COURT check that these 12 lines create a deep point
		return trVerts[0];
	if(uvParameters[0]>0.9998f)
		return trVerts[1];
	if(uvParameters[1]>0.9998f)
		return trVerts[2];
	if(uvParameters[0]<0.0002f)
		return splitTriangleEdge(triangle,2,1.0f-uvParameters[1]);
	if(uvParameters[1]<0.0002f)
		return splitTriangleEdge(triangle,0,uvParameters[0]);
	if(uvParameters[0]+uvParameters[1]>0.9998f)
		return splitTriangleEdge(triangle,1,1.0f-uvParameters[0]);
	// now we know we will add a vertex and two triangles. These operations invalidate pointers due to possible reallocation.
	int* trTex = triangleTextures(triangle);
	int v[3],v0=trVerts[0],v1=trVerts[1],oldVert=trVerts[2], tx0 = trTex[0], tx1 = trTex[1], oldTx = trTex[2];
	int ret = addVertices(), rTx = addTexture();
	unsigned int *trAdjs = &_adjs[triangle*3];
	unsigned int a1=trAdjs[1],a2=trAdjs[2];
	float p,*pv = vertexCoordinate(v0);
	p = 1.0f - uvParameters[0] - uvParameters[1];
	float vec[3] = {p*pv[0],p*pv[1],p*pv[2]};
	pv = getTexture(trTex[0]);  // vertexTexture(v0);
	float tx[2]={p*pv[0],p*pv[1]};
	for(int i=0; i<2; ++i)	{
		pv = getTexture(trTex[i+1]);
		tx[0] += pv[0]*uvParameters[i];
		tx[1] += pv[1]*uvParameters[i];
		pv = vertexCoordinate(trVerts[i+1]);
		vec[0]+=pv[0]*uvParameters[i]; vec[1]+=pv[1]*uvParameters[i]; vec[2]+=pv[2]*uvParameters[i];
	}
	setVertexCoordinate(ret,vec);
	setTexture(rTx,tx);
	// assign vertices
	trVerts[2] = ret;
	trTex[2] = rTx;
	int t[3];
	v[0] = ret; v[1] = v1; v[2] = oldVert;
	t[0] = rTx; t[1] = trTex[1]; t[2] = oldTx;
	int t2, t1 = addTriangle(v, _tris[triangle].material, t);  // invalidates _tris and _adj pointers and iterators
	trTex = triangleTextures(triangle);
	v[0] = ret; v[2] = v0; v[1] = oldVert;
	t[0] = rTx; t[2] = trTex[0]; t[1] = oldTx;
	t2 = addTriangle(v, _tris[triangle].material, t);  // invalidates _tris and _adj pointers and iterators
	// assign adjs
	_adjs[triangle*3+1] = t1<<2;
	_adjs[triangle*3+2] = (t2<<2)+2;
	_adjs[t1*3] = (triangle<<2)+1;
	_adjs[t1*3+1] = a1;
	_adjs[t1*3+2] = t2<<2;
	_adjs[t2*3] = (t1<<2)+2;
	_adjs[t2*3+1] = a2;
	_adjs[t2*3+2] = (triangle<<2)+2;
	if(a1!=0x00000003)
		_adjs[(a1>>2)*3+(a1&0x00000003)] = (t1<<2)+1;
	if(a2!=0x00000003)
		_adjs[(a2>>2)*3+(a2&0x00000003)] = (t2<<2)+1;
	// assign vertexFace
	if((_vertexFace[oldVert]&0x3fffffff)==triangle)	{
		_vertexFace[oldVert]=t2;
		if(a2==0x00000003)
			_vertexFace[oldVert] |= 0x40000000;
	}
	if((_vertexFace[v1]&0x3fffffff)==triangle)	{
		_vertexFace[v1]=t1;
		if(a1==0x00000003)
			_vertexFace[v1] |= 0x40000000;
	}
	_vertexFace[ret]=triangle;
	// 	_vertexFace[v0] can stay unchanged
	return ret;
}

/* int materialTriangles::addTriangle(int(&vertices)[3], int material)
{  // now using texture seams making this simple call invalid.
	throw(std::logic_error("This addTriangle() call is obsolete.\n"));
	int retval = (int)_tris.size();
	matTriangle mt;
	mt.v[0] = vertices[0]; mt.v[1] = vertices[1]; mt.v[2] = vertices[2];
	mt.tex[0] = 0;  mt.tex[1] = 0; mt.tex[2] = 0;
	mt.material = material;
	_tris.push_back(mt);
	if(!_adjs.empty()) {
		_adjs.push_back(0x03);
		_adjs.push_back(0x03);
		_adjs.push_back(0x03); }
	_adjacenciesComputed = false;
	return retval;
} */

int materialTriangles::addTriangle(const int(&vertices)[3], const int material,  const int(&textures)[3])
{
	int retval = (int)_tris.size();
	matTriangle mt;
	mt.v[0] = vertices[0]; mt.v[1] = vertices[1]; mt.v[2] = vertices[2];
	mt.tex[0] = textures[0];  mt.tex[1] = textures[1]; mt.tex[2] = textures[2];
	mt.material = material;
	_tris.push_back(mt);
	if (!_adjs.empty()) {
		_adjs.push_back(0x03);
		_adjs.push_back(0x03);
		_adjs.push_back(0x03);
	}
	_adjacenciesComputed = false;
	return retval;
}

int materialTriangles::addVertices(int numberToAdd)
{
	int retval = (int)_xyz.size() / 3;
	for(unsigned int i=0; i<(unsigned)numberToAdd; ++i)
	{
		_xyz.push_back(0.0f);
		_xyz.push_back(0.0f);
		_xyz.push_back(0.0f);
		// in new version texture not present.  Force addTexture() instead
		if(!_vertexFace.empty())
			_vertexFace.push_back(0x80000000);
	}
	_adjacenciesComputed = false;
	return retval;
}

int materialTriangles::cloneVertex(int sourceVertex)
{  // makes a new vertex which is a copy of sourceVertex
	int retval = (int)_xyz.size() / 3;
	float V[3];
	getVertexCoordinate(sourceVertex, V);
	_xyz.push_back(V[0]);
	_xyz.push_back(V[1]);
	_xyz.push_back(V[2]);
	if (!_vertexFace.empty())
		_vertexFace.push_back(0x80000000);
	_adjacenciesComputed = false;
	return retval;
}

void materialTriangles::setVertexCoordinate(int vertex, const float (&newCoord)[3])
{
	float *v = &_xyz[(vertex<<1)+vertex];
	v[0]=newCoord[0];
	v[1]=newCoord[1];
	v[2]=newCoord[2];
}

void materialTriangles::cleanAndPack()
{	// warning - invalidates all triangle and vertex indices.
	_vertexFace.clear();
	_vertexFace.assign(_xyz.size()/3,0x80000000);
	_adjs.clear();
	int i,n=(int)_tris.size();
	std::vector<matTriangle> tmpTr;
	tmpTr.reserve(n);
	int *tr;
	for(i=0; i<n; ++i)	{
		if(_tris[i].material<0)	// deleted triangle
			continue;
		tmpTr.push_back(_tris[i]);
		tr = &_tris[i].v[0];
		for (int j = 0; j<3; ++j)
			_vertexFace[tr[j]] = 0;
	}
	_tris.assign(tmpTr.begin(), tmpTr.end());
	tmpTr.clear();
	n = (int)_vertexFace.size();
	int bot=0,top,newV=0;
	for(i=0; i<n; ++i)	{
		if(_vertexFace[i]>0)
			continue;
		_vertexFace[i] = newV++;
		top = (i<<1)+i;
		if(bot<top) {
			for(int j=0; j<3; ++j)
				_xyz[bot++]=_xyz[top++];
		}
		else
			bot += 3;
	}
	_xyz.resize(bot);
	bot = 0;
	for (i = 0; i<n; ++i)	{
		if (_vertexFace[i]>0x7fffffff)
			continue;
		top = i << 1;
		if (bot<top) {
			_uv[bot++] = _uv[top++];
			_uv[bot++] = _uv[top++];
		}
		else
			bot += 2;
	}
	_uv.resize(bot);
	n = (int)_tris.size();
	for (i = 0; i < n; ++i)
		for (int j = 0; j < 3; ++j)
			_tris[i].v[j] = _vertexFace[_tris[i].v[j]];
	_adjacenciesComputed=false;
	_vertexFace.clear();
}

void materialTriangles::cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap)
{	// warning - invalidates all triangle and vertex indices. Returns the new index numbers of the old vertices and triangles. -1 indicates deletion.
	std::vector<int> pp;
	newVertexMap.clear();
	newVertexMap.assign(_xyz.size()/3,-1);
	pp.assign(_xyz.size()/3,-1);
	int i,j,k,n=(int)_tris.size();
	for(i=0; i<n; ++i)	{
		if(_tris[i].material<0)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			newVertexMap[_tris[i].v[j]] = 1;
			pp[_tris[i].v[j]] = 1;
		}
	}
	n=(int)pp.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(pp[i]<0)
			continue;
		if(i>j)	{
			_xyz[(j<<1)+j]=_xyz[(i<<1)+i];
			_xyz[(j<<1)+j+1]=_xyz[(i<<1)+i+1];
			_xyz[(j<<1)+j+2]=_xyz[(i<<1)+i+2];
		}
		pp[i]=j;
		++j;
	}
	_xyz.resize((j<<1)+j);
	n = (int)newVertexMap.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(newVertexMap[i]<0)
			continue;
		_uv[(j<<1)+j]=_uv[(i<<1)+i];
		_uv[(j<<1)+j+1]=_uv[(i<<1)+i+1];
		newVertexMap[i]=j;
		++j;
	}
	_uv.resize(j<<1);
	n=(int)_tris.size();
	newTriangleMap.clear();
	newTriangleMap.assign(n,-1);
	k=0;
	for(i=0; i<n; ++i)	{
		if(_tris[i].material<0)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			assert(newVertexMap[_tris[i].v[j]]>-1);
			_tris[k].v[j] = newVertexMap[_tris[i].v[j]];
		}
		newTriangleMap[i] = k++;
	}
	_tris.resize(k);
	_adjacenciesComputed=false;
	_adjs.clear();
	_vertexFace.clear();
}

int materialTriangles::isManifoldConsistent()
{	// Closed manifold surface topology checker.  Returns # of topological handles or -1 if inconsistent
	typedef std::set<edge,edgeTest> edgeSet;
	edge E;
	edgeSet M;
	M.clear();
	int tnow[3];
	int i,j,numtris=(int)_tris.size();
	std::vector<int> posVec;
	posVec.assign(_xyz.size()/3,0);
	for(i=0; i<numtris; ++i)
	{
		if(_tris[i].material<0)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
			tnow[j] = _tris[i].v[j];
		for(j=0; j<3; j++)
		{
			if (tnow[j] >= (int)posVec.size())
				return -1;
			else
				posVec[tnow[j]] = 1;
			unsigned int tmp;
			if((tmp=(unsigned int)tnow[(j+1)%3])<(unsigned int)tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j];	}
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp;	}
			M.insert(E);
		}
	}
	for (i = 0; i<(int)posVec.size(); ++i)
		if(!posVec[i])
			return -1;
	int handles2 = (int)((_tris.size()) + posVec.size() - M.size());
	if(handles2&0x0001)
		return -1;
	return handles2>>1;
}

void materialTriangles::getMeanVertexNormal(int vertex, float(&normal)[3], int onlyMaterial)
{  // if onlyMaterial>-1 only use neighbor triangles with material == onlyMaterial
	if(!_adjacenciesComputed)	findAdjacentTriangles();
	std::vector<neighborNode> nei;
	getNeighbors(vertex,nei);
	if(nei.size()<2) {
		normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;
		return; }
	std::vector<neighborNode>::iterator nit=nei.begin();
	Vec3f last,now,p,mean(0.0f,0.0f,0.0f);
	getVertexCoordinate(vertex,p._v);
	int lastV;
	if(nei.front().triangle<0) {
		lastV = nei.front().vertex;
		++nit;
	}
	else
		lastV = nei.back().vertex;
	getVertexCoordinate(lastV, last._v);
	last -= p;
	while(nit!=nei.end()) {
		getVertexCoordinate(nit->vertex,now._v);
		now -= p;
		if (onlyMaterial<0 || (onlyMaterial > -1 && _tris[nit->triangle].material == onlyMaterial))
			mean += last^now;
		last = now;
		++nit;
	}
	mean.normalize();
	normal[0]=mean._v[0]; normal[1]=mean._v[1]; normal[2]=mean._v[2];
}

bool materialTriangles::deleteEdge(int triangle, int edge)
{ // always leaves triangle vertex[edge] behind and deletes vertex[edge+1]

	throw(std::logic_error("Needs software fix."));  // COURT needs to be fixed for new texture seam version

	if (_tris[triangle].material < 0)
		return false;
	unsigned int adj, *vf, *vf1, ae, ae1, ae2;
	ae1 = _adjs[triangle * 3 + ((edge + 1) % 3)];
	ae2 = _adjs[triangle * 3 + ((edge + 2) % 3)];
	if (ae1 != 3 && ae2 != 3 && triangleVertices(ae1 >> 2)[((ae1 & 3) + 2) % 3] == triangleVertices(ae2 >> 2)[((ae2 & 3) + 2) % 3])
		return false;  // would produce an illegal surface
	ae = _adjs[(triangle << 1) + triangle + edge];
	if (ae != 3) {
		assert(_tris[ae >> 2].material > -1);
		ae1 = _adjs[(ae >> 2) * 3 + (((ae & 3) + 1) % 3)];
		ae2 = _adjs[(ae >> 2) * 3 + (((ae & 3) + 2) % 3)];
		if (ae1 != 3 && ae2 != 3 && triangleVertices(ae1 >> 2)[((ae1 & 3) + 2) % 3] == triangleVertices(ae2 >> 2)[((ae2 & 3) + 2) % 3])
			return false;  // would produce an illegal surface
	}
	int ve, ve1;  // variables starting with v are vertices, t are triangles, vf are vertexFace, a are adjacencies
	ve = _tris[triangle].v[edge];
	ve1 = _tris[triangle].v[(edge + 1) % 3];
	vf = &_vertexFace[ve];
	vf1 = &_vertexFace[ve1];
	if (*vf & 0x40000000) {
		if (*vf1 & 0x40000000) {
			if (ae != 3)  // would create an illegal surface
				return false;
			else {
				assert((*vf & 0x3fffffff) == triangle);
				if ((*vf1 & 0x3fffffff) == triangle) {
					adj = _adjs[(triangle << 1) + triangle + ((edge + 2) % 3)];
					if (adj == 3) {  // single isolated triangle
						_vertexFace[_tris[triangle].v[0]] = 0x80000000;  // strand vertices
						_vertexFace[_tris[triangle].v[1]] = 0x80000000;
						_vertexFace[_tris[triangle].v[2]] = 0x80000000;
						_tris[triangle].material = -1;
						_adjs[triangle * 3] = 3;
						_adjs[triangle * 3 + 1] = 3;
						_adjs[triangle * 3 + 2] = 3;
						return true;
					}
					else
						*vf = (adj >> 2) | 0x40000000;
				}
				else
					*vf = *vf1;
			}
		}
		else {
			assert((*vf & 0x3fffffff) != triangle && ae != 3);
			if ((*vf & 0x3fffffff) == (ae >> 2)) {
				assert(ae1 == 3 && ae2 != 3);
				*vf = (ae2 >> 2) | 0x40000000;
			}
			// else leave as is
		}
	}
	else if (*vf1 & 0x40000000) {
		assert(ae != 3);
		if ((*vf1 & 0x3fffffff) == triangle) {
			adj = _adjs[(triangle << 1) + triangle + ((edge + 2) % 3)];
			assert(adj != 3);
			*vf = (adj >> 2) | 0x40000000;
		}
		else
			*vf = *vf1;
	}
	else {
		assert(ae != 3);
		if ((*vf & 0x3fffffff) == triangle || (*vf & 0x3fffffff) == (ae >> 2)) {
			adj = _adjs[(triangle << 1) + triangle + ((edge + 2) % 3)];
			assert(adj != 3);
			*vf = (adj >> 2);
		}
		// else leave as is
	}
	*vf1 = 0x80000000;  // strand
	// default behavior is remaining vertex is an average of the two. If you want something else (e.g. volume preservation) compute externally.
Vec3f vve, vve1;
vve.set((const float(&)[3])_xyz[ve * 3]);
vve1.set((const float(&)[3])_xyz[ve1 * 3]);
vve += vve1;
vve *= 0.5f;
setVertexCoordinate(ve, vve._v);
float uv0[2], uv1[2];
assert(false);
//	getVertexTexture(ve, uv0);
//	getVertexTexture(ve1, uv1);
uv0[0] += uv1[0];  uv0[1] += uv1[1];
uv0[0] *= 0.5f;  uv0[1] *= 0.5f;
//	setVertexTexture(ve, uv0);
auto delOneTri = [&](int tH, int eH) {
	_tris[tH].material = -1;  // mark deleted
	unsigned int* vf2;
	vf2 = &_vertexFace[_tris[tH].v[(eH + 2) % 3]];
	ae1 = _adjs[(tH << 1) + tH + ((eH + 1) % 3)];
	ae2 = _adjs[(tH << 1) + tH + ((eH + 2) % 3)];
	_adjs[(tH << 1) + tH] = 3;
	_adjs[(tH << 1) + tH + 1] = 3;
	_adjs[(tH << 1) + tH + 2] = 3;
	if (ae2 != 3)
		_adjs[(ae2 >> 2) * 3 + (ae2 & 3)] = ae1;
	if (ae1 != 3)
		_adjs[(ae1 >> 2) * 3 + (ae1 & 3)] = ae2;
	if (*vf2 & 0x40000000) {
		if ((*vf2 & 0x3fffffff) == tH) {
			assert(ae2 == 3);
			if (ae1 != 3)
				*vf2 = (ae1 >> 2) | 0x40000000;
			else
				*vf2 = 0x80000000; // strand deleted vertex
		}
		// else leave as is
	}
	else {  // part of closed ring
		if ((*vf2 & 0x3fffffff) == tH) {
			assert(ae1 != 3);
			*vf2 = (ae1 >> 2);
		}
	}
};
if (ae == 3)  // no opposite triangle case
delOneTri(triangle, edge);
else {
	delOneTri(triangle, edge);
	delOneTri(ae >> 2, ae & 3);
}
assert(!(*vf & 0x80000000));
int trNum = *vf & 0x3fffffff;
unsigned int stop;
int* tnow = &(_tris[trNum].v[0]);
assert(_tris[trNum].material > -1);
int j;
for (j = 0; j < 3; ++j)
	if (tnow[j] == ve || tnow[j] == ve1) {
		tnow[j] = ve;
		break;
	}
assert(j < 3);
// set triStart to the end adjacency code for counterclockwise traversal
if (*vf & 0x40000000)	// started on a free edge, will end on one
stop = 3;
else	// create adjacency code of starting edge
stop = (trNum << 2) + j;
adj = _adjs[(trNum << 1) + trNum + ((j + 2) % 3)];
while (adj != stop)
{
	trNum = adj >> 2;
	j = adj & 0x00000003;
	if (_tris[trNum].v[j] == ve1)
		_tris[trNum].v[j] = ve;
	else
		assert(_tris[trNum].v[j] = ve);
	adj = _adjs[(trNum << 1) + trNum + ((j + 2) % 3)];
}
return true;
}

bool materialTriangles::hasSelfIntersection(const bool isClosed, std::vector<std::pair<int, int> > &triangleIntersectPairs) {  // slow. should only be used to debug surface
	triangleIntersectPairs.clear();
	boundingBox<float> eb, tb;
	for (int n = numberOfTriangles(), i = 0; i < n; ++i) {
		if (triangleMaterial(i) < 0)
			continue;
		int* tr = triangleVertices(i);
		Vec3f e[2], t[3];
		int e0, e1;
		for (int j = 0; j < 3; ++j) {
			if (isClosed && tr[j] > tr[(j + 1) % 3])  // only need to do this edge once
				continue;
			e0 = tr[j];
			e1 = tr[(j + 1) % 3];
			getVertexCoordinate(e0, e[0]._v);
			getVertexCoordinate(e1, e[1]._v);
			eb.Empty_Box();
			eb.Enlarge_To_Include_Point(e[0]._v);
			eb.Enlarge_To_Include_Point(e[1]._v);
			e[1] -= e[0];
			for (int m, k = 0; k < n; ++k) {
				if (k == i || triangleMaterial(k) < 0)
					continue;
				int* t2 = triangleVertices(k);
				tb.Empty_Box();
				for (m = 0; m < 3; ++m){
					if (t2[m] == e0 || t2[m] == e1)
						break;
					getVertexCoordinate(t2[m], t[m]._v);
					tb.Enlarge_To_Include_Point(t[m]._v);
				}
				if (m < 3 || !eb.Intersection(tb))
					continue;
				t[1] -= t[0];
				t[2] -= t[0];
				Mat3x3f M;
				M.Initialize_With_Column_Vectors(t[1], t[2], -e[1]);
				Vec3f R = M.Robust_Solve_Linear_System(e[0] - t[0]);
				if (R.X <= 0.0f || R.Y <= 0.0f || R.Z <= 0.0f || R.X > 1.0f || R.Y >= 1.0f || R.Z >= 1.0f || R.X + R.Y >= 1.0f)
					continue;
				std::pair<int, int> ip;
				ip.first = k;
				ip.second = i;
				triangleIntersectPairs.push_back(ip);
			}
		}
	}
	if (triangleIntersectPairs.empty()) {
		std::cout << "materialTriangle object does not have self intersections\n";
		return false;
	}
	else {
		std::cout << "materialTriangle object has " << triangleIntersectPairs.size() << " self intersections-\n";
		return true;
	}
}

bool materialTriangles::topoCheck()
{
    std::vector<unsigned int> adjs;	// low 2 bits are the edge number of the adjacent triangle.
	adjs.assign(_adjs.begin(),_adjs.end());
	std::vector<unsigned int> vertexFace;
	vertexFace.assign(_vertexFace.begin(),_vertexFace.end());
	_adjacenciesComputed=false; // force computation
	findAdjacentTriangles();
	int i,n;
	n = (int)adjs.size();
	for(i=0; i<n; ++i) {
		if((i%3==0) && _tris[i/3].material<0) {
			i+=2;
			continue; }
		if(adjs[i]!=_adjs[i])
			return false;
	}
	std::vector<materialTriangles::neighborNode> nei;
	std::vector<materialTriangles::neighborNode>::iterator nit;
	n = (int)vertexFace.size();
	for(i=0; i<n; ++i) {
		unsigned int vo,vf = vertexFace[i];
		vo = _vertexFace[i];
		if(vo>0xfffffffe || vf>0xfffffffe) {
			if(vo!=vf) // should be stranded vertex
				return false;
			else
				continue;
		}
		if((vo&0x40000000)>0) {
			if(vo!=vf)
				return false;
		}
		else {
			getNeighbors(i,nei);
			for(nit=nei.begin(); nit!=nei.end(); ++nit) {
				if(nit->triangle==vf)
					break;
			}
			if(nit==nei.end()) { // could be stranded vertex
				int j,k;
				for(j=0; j<(int)_tris.size(); ++j) {
					if(_tris[j].material<0)
						continue; 
					for (k = 0; k < 3; ++k) {
						if (_tris[j].v[k] == i)
							break;
					}
					if (k < 3)
						break;
				}
				if(j<n) // not a stranded vertex
					return false;
			}
		}
	}
	return true;
}

bool materialTriangles::inside(const float(&xyz)[3])
{  // if a closed manifold surface, returns if xyz is inside
	std::vector<unsigned char> trDone;
	int nI=0,n = (int)_tris.size();
	int *t;
	trDone.assign(n, 0x00);
	float *p, *q, *r,U[3],V[3],X[2],u,v,c;
	for (int i = 0; i < n; ++i) {
		if (trDone[i]>0x00)
			continue;
		t = &_tris[i].v[0];
		p = vertexCoordinate(t[0]);
		q = vertexCoordinate(t[1]);
		r = vertexCoordinate(t[2]);
		if ((xyz[0]<p[0] && xyz[0]<q[0] && xyz[0]<r[0]) || (xyz[0]>p[0] && xyz[0]>q[0] && xyz[0]>r[0]))
			continue;
		if ((xyz[1]<p[1] && xyz[1]<q[1] && xyz[1]<r[1]) || (xyz[1]>p[1] && xyz[1]>q[1] && xyz[1]>r[1]))
			continue;
		U[0] = q[0] - p[0]; U[1] = q[1] - p[1]; U[2] = q[2] - p[2];
		V[0] = r[0] - p[0]; V[1] = r[1] - p[1]; V[2] = r[2] - p[2];
		X[0] = xyz[0] - p[0]; X[1] = xyz[1] - p[1];
		if (fabs(V[1]) > fabs(V[0])) {
			c = V[0] / V[1];
			u = X[0] - c*X[1];
			c = U[0] - c*U[1];
			if (fabs(c) < 1e-16f)
				continue;
			u /= c;
			v = (X[1] - u*U[1]) / V[1];
		}
		else {
			c = V[1] / V[0];
			u = X[1] - c*X[0];
			c = U[1] - c*U[0];
			if (fabs(c) < 1e-16f)
				continue;
			u /= c;
			v = (X[0] - u*U[0]) / V[0];
		}
		if (u<0.0f || v<0.0f || u>1.0f || v>1.0f || u + v > 1.0f)
			continue;
		if (v == 0.0f) {
			if (!_adjacenciesComputed)
				findAdjacentTriangles();
			trDone[_adjs[i * 3] >> 2] = 0xff;
		}
		if (u == 0.0f) {
			if (!_adjacenciesComputed)
				findAdjacentTriangles();
			trDone[_adjs[i * 3 + 2] >> 2] = 0xff;
		}
		if (u + v == 1.0f) {
			if (!_adjacenciesComputed)
				findAdjacentTriangles();
			trDone[_adjs[i * 3 + 1] >> 2] = 0xff;
		}
		c = p[2] + u*U[2] + v*V[2];
		if (c > xyz[2])
			++nI;
	};
	return (nI&0x01)>0;
}

void materialTriangles::partitionTriangleMaterials()
{  // returns the indices into triangle array of the next element beyond each key first material
	std::vector<matTriangle>::iterator prev,sit = _tris.begin(), eit = _tris.end();
	int matNow = 0;
	prev = sit;
	while (true) {
		sit = std::stable_partition(sit, eit, [&matNow](matTriangle &mt){ return mt.material < matNow; });
		if (sit != eit) {
			prev = sit;
			++matNow; }
		else
			break;
	}
	_adjacenciesComputed = false;
}

/* int materialTriangles::getMaterialEnd(const int material)
{ // after partitionMaterials(), input material and will return the index of one past the last triangle of this material. Return 0 means no material found. Return -1 means material requested exceeds maximum material.
	if (material>_matEnds.rbegin()->first)
		return -1L;
	auto mit = _matEnds.find(material);
	if (mit == _matEnds.end())
		return 0;
	return mit->second;
} */

void materialTriangles::closestPoint(const float(&xyz)[3], int &triangle, float(&uv)[2], int onlyMaterial)  // closest barycentric position to point xyz
{
	Vec3f P, Q, R;
	float t, minT, uvNow[2], dsq, minDsq = 1e32f;
	unsigned int bestE;
	P.set(xyz);
	// examine only unique edges
	unsigned int edge,adj;
	for (int n = (int)_tris.size(), j, i = 0; i < n; ++i) {
		if (_tris[i].material < 0 || (onlyMaterial > -1 && _tris[i].material != onlyMaterial))
			continue;
		for (j = 0; j < 3; ++j) {
			edge = (i << 2) + j;
			adj = _adjs[(i << 1) + i + j];
			if (adj != 3) {
				if (adj < edge && (onlyMaterial < 0 || _tris[adj >> 2].material == onlyMaterial)) // only do this edge once
					continue;
			}
			Q.set((const float(&)[3])_xyz[_tris[i].v[j] * 3]);
			R.set((const float(&)[3])_xyz[_tris[i].v[(j + 1) % 3] * 3]);
			R -= Q;
			dsq = R*R;
			if (dsq < 1e-16f)
				continue;
			t = ((P - Q)*R) / dsq;
			if (t > 1.0f)
				t = 1.0f;
			if (t < 0.0f)
				t = 0.0f;
			Q += R*t;
			dsq = (P - Q).length2();
			if (dsq < minDsq) {
				minDsq = dsq;
				bestE = edge;
				minT = t;
			}
		}
	}
	triangle = bestE >> 2;
	getBarycentricProjection(triangle, xyz, uv);
	bool midTri = false;
	if (uv[0] >= 0.0f && uv[1] >= 0.0f && uv[0] + uv[1] <= 1.0f) {
		getBarycentricPosition(triangle, uv, R._v);
		if ((dsq = (P - R).length2()) < minDsq) {
			minDsq = dsq;
			midTri = true;
		}
	}
	edge = _adjs[triangle * 3 + (bestE & 3)];
	if (edge != 3 && (onlyMaterial < 0 || _tris[edge >> 2].material == onlyMaterial)) {
		getBarycentricProjection(edge >> 2, xyz, uvNow);
		if (uvNow[0] >= 0.0f && uvNow[1] >= 0.0f && uvNow[0] + uvNow[1] <= 1.0f) {
			getBarycentricPosition(edge >> 2, uvNow, R._v);
			if ((P - R).length2() < minDsq) {
				midTri = true;
				triangle = edge >> 2;
				uv[0] = uvNow[0];
				uv[1] = uvNow[1];
			}
		}
	}
	if (midTri)
		return;
	if ((bestE & 3) < 1) {
		uv[0] = minT;
		uv[1] = 0.0f;
	}
	else if ((bestE & 3) < 2) {
		uv[0] = 1.0f - minT;
		uv[1] = minT;
	}
	else {
		uv[1] = 1.0f - minT;
		uv[0] = 0.0f;
	}
}

bool materialTriangles::textureFind(const float(&txIn)[2], const int materialIn, int &triangle, float(&uv)[2])
{  // finds the triangle in mt and its uv parametric location closest to texture txIn and material materialIn.
	// if direct hit, returns true. If no direct hit will find closest and return false;
	float dist2 = 1e32f, uvBest[2], uvNow[2], txNow[2];  // tBox[4], *tx,
	int *tv, bestTri;
	bool doDist2;
	for (int j, n = numberOfTriangles(), i = 0; i < n; ++i) {
		if (materialIn>-1 && _tris[i].material != materialIn)
			continue;
		tv = _tris[i].v;
		// now that we decided to look for closest texture, consider nuking bounding box
		/*		tBox[0] = tBox[2] = 1e32f;
		tBox[1] = tBox[3] = -1e32f;
		for (j = 0; j < 3; ++j) {
		tx = &mt->_uv[tv[j]<<1];
		if (tx[0] < tBox[0])
		tBox[0] = tx[0];
		if (tx[0] > tBox[1])
		tBox[1] = tx[0];
		if (tx[1] < tBox[2])
		tBox[2] = tx[1];
		if (tx[1] > tBox[3])
		tBox[3] = tx[1];
		}
		if (txIn[0]<tBox[0] || txIn[0]>tBox[1] || txIn[1]<tBox[2] || txIn[1]>tBox[3])
		continue; */
		float det, a0, a1, b0, b1, c0, c1, *trx[3];
		for (j = 0; j < 3; ++j)
			trx[j] = &_uv[tv[j] << 1];
		a0 = trx[1][0] - trx[0][0];
		a1 = trx[1][1] - trx[0][1];
		b0 = trx[2][0] - trx[0][0];
		b1 = trx[2][1] - trx[0][1];
		c0 = txIn[0] - trx[0][0];
		c1 = txIn[1] - trx[0][1];
		det = a0*b1 - a1*b0;
		if (fabs(det) < 1e-16f)
			continue;
		uv[0] = (c0*b1 - c1*b0) / det;
		uv[1] = (a0*c1 - a1*c0) / det;
		if (uv[0] > -1e-5f && uv[1] > -1e-5f && uv[0] < 1.0001f && uv[1] < 1.0001f && uv[0] + uv[1] < 1.0001f) {
			triangle = i;
			if (uv[0] < 0.0f)
				uv[0] = 0.0f;
			if (uv[1] < 0.0f)
				uv[1] = 0.0f;
			if (uv[0] > 1.0f)
				uv[0] = 1.0f;
			if (uv[1] > 1.0f)
				uv[1] = 1.0f;
			if (uv[0] + uv[1] > 1.0f) {
				a0 = uv[0] + uv[1];
				uv[0] /= a0;
				uv[1] /= a0;
			}
			return true;
		}
		doDist2 = false;
		if (uv[0] < 0.0f) {
			doDist2 = true;
			uvNow[0] = 0.0f;
			assert(b0*b0 + b1*b1 > 1e-16f);
			uvNow[1] = (c0*b0 + c1*b1) / (b0*b0 + b1*b1);
			if (uvNow[1] < 0.0f)
				uvNow[1] = 0.0f;
			else if (uvNow[1] > 1.0f)
				uvNow[1] = 1.0f;
			else
				;
		}
		if (uv[1] < 0.0f) {
			doDist2 = true;
			uvNow[1] = 0.0f;
			assert(a0*a0 + a1*a1 > 1e-16f);
			uvNow[0] = (c0*a0 + c1*a1) / (a0*a0 + a1*a1);
			if (uvNow[0] < 0.0f)
				uvNow[0] = 0.0f;
			else if (uvNow[0] > 1.0f)
				uvNow[0] = 1.0f;
			else
				;
		}
		if (!doDist2 && uv[0] + uv[1] > 1.0f) {
			doDist2 = true;
			b0 = trx[2][0] - trx[1][0];
			b1 = trx[2][1] - trx[1][1];
			c0 = txIn[0] - trx[1][0];
			c1 = txIn[1] - trx[1][1];
			assert(b0*b0 + b1*b1 > 1e-16f);
			uvNow[1] = (c0*b0 + c1*b1) / (b0*b0 + b1*b1);
			if (uvNow[1] < 0.0f)
				uvNow[1] = 0.0f;
			else if (uvNow[1] > 1.0f)
				uvNow[1] = 1.0f;
			else
				;
			uvNow[0] = 1.0f - uvNow[1];
		}
		assert(doDist2);
		getBarycentricTexture(i, uvNow, txNow);
		a0 = txIn[0] - txNow[0];
		a1 = txIn[1] - txNow[1];
		b1 = a0*a0 + a1*a1;
		if (b1 < dist2) {
			dist2 = b1;
			uvBest[0] = uvNow[0];
			uvBest[1] = uvNow[1];
			bestTri = i;
		}
	}
	triangle = bestTri;
	uv[0] = uvBest[0];
	uv[1] = uvBest[1];
	return false;
}

void materialTriangles::clear()
{
	_tris.clear();
//	_matEnds.clear();
	_xyz.clear();
	_uv.clear();
	_adjacenciesComputed= false;
	_adjs.clear();
	_vertexFace.clear();
	_name.assign("");
}

bool materialTriangles::checkTopology()
{  // debugging routine to check veracity of current topological arrays
	std::vector<unsigned int> adjs, vf;	// low 2 bits are the edge number of the adjacent triangle.
	adjs.assign(_adjs.begin(), _adjs.end());
	vf.assign(_vertexFace.begin(), _vertexFace.end());
	if (_vertexFace.size() != _xyz.size() / 3)
		return false;
	if (_adjs.size() != _tris.size()*3)
		return false;
	findAdjacentTriangles(true);
	if(adjs.size() != _adjs.size())
		return false;
	if(vf.size() != _vertexFace.size())
		return false;
	for (int n = (int)adjs.size(), i = 0; i < n; ++i) {
		if (adjs[i] != _adjs[i]) {
			return false;
		}
	}
	std::vector<materialTriangles::neighborNode> nei;
	std::vector<materialTriangles::neighborNode>::iterator nit;
	for (int n = (int)vf.size(), i = 0; i < n; ++i) {
		if (vf[i] != _vertexFace[i]) {
			if ((vf[i] & 0x80000000) != (_vertexFace[i] & 0x80000000))
				return false;
			if ((vf[i] & 0x80000000) && (_vertexFace[i] != vf[i]))
				return false;
			if ((vf[i] & 0x40000000) != (_vertexFace[i] & 0x40000000))
				return false;
			if ((vf[i] & 0x40000000) && (_vertexFace[i] != vf[i]))
				return false;
			if ((_vertexFace[i] & 0xc0000000) < 1) {
				getNeighbors(i, nei);
				assert(nei.front().triangle > -1);
				for (nit = nei.begin(); nit != nei.end(); ++nit) {
					if (nit->triangle == vf[i])
						break;
				}
				if (nit == nei.end())
					return false;
			}
		}
	}
	return true;
}

void materialTriangles::correctLocalNeighborArrays(std::vector<int> &changedTriangles)
{  // does local patching of adjacency and vertexFace arrays for input changedTriangles and the triangles adjacent to them.
	// assumes adjacencies have been computed before triangle vertices were changed.
	std::map<std::pair<int, int>, unsigned int> M;
	std::vector<unsigned int> TE;
	TE.reserve(changedTriangles.size() * 3);
	std::sort(changedTriangles.begin(), changedTriangles.end());
	std::set<int> verts;
	for (auto ct : changedTriangles){
		if (_tris[ct].material < 0)
			continue;
		for (int i = 0; i < 3; ++i){
			int v = _tris[ct].v[i];
			_vertexFace[v] = ct;
			verts.insert(v);
			unsigned int adj = _adjs[ct * 3 + i];
			if (adj == 3)
				continue;
			M.insert(std::make_pair(std::make_pair(_tris[ct].v[(i + 1) % 3], _tris[ct].v[i]), (ct << 2) + i));
			if (!std::binary_search(changedTriangles.begin(), changedTriangles.end(), (int)(adj>>2)))
				TE.push_back(adj);
			else
				_adjs[ct * 3 + i] = 3;
		}
	}
	for (auto te : TE){
		auto tep = std::make_pair(_tris[te >> 2].v[te & 3], _tris[te >> 2].v[((te & 3)+1)%3]);
		auto mit = M.find(tep);
		assert(mit != M.end());
		_adjs[(te >> 2) * 3 + (te & 3)] = mit->second;
		_adjs[(mit->second >> 2) * 3 + (mit->second & 3)] = te;
		M.erase(mit);
	}
	for (auto ct : changedTriangles){
		if (_tris[ct].material < 0)
			continue;
		for (int i = 0; i < 3; ++i){
			auto tep = std::make_pair(_tris[ct].v[i], _tris[ct].v[(i + 1) % 3]);
			auto mit = M.find(tep);
			if (mit != M.end()){
				_adjs[ct * 3 + i] = mit->second;
				_adjs[(mit->second >> 2) * 3 + (mit->second & 3)] = (ct<<2) + i;
				M.erase(mit);
			}
		}
	}
	assert(M.empty());
	for (auto v : verts){
		unsigned int teStart, te, te2;
		for (int i = 0; i < 3; ++i){
			if (_tris[_vertexFace[v]].v[i] == v){
				teStart = (_vertexFace[v] << 2) + i;
				break;
			}
		}
		te = teStart;
		do{
			te2 = _adjs[(te >> 2) * 3 + (te & 3)];
			if (te2 == 3){
				_vertexFace[v] = (te >> 2) | 0x40000000;
				break;
			}
			te = (te2&0xfffffffc) + (((te2&3)+1)%3);
		} while (te != teStart);
	}
}

void materialTriangles::collectCreateTextureSeams() {
	// only called on load.  Finds same material texture seam txVerts. Gets or creates when necessary intermaterial texture vertices,
	// giving texture priority to material 2 then 7. Material 1 is unimportant as these triangles colored blue and are unselectable.
	findAdjacentTriangles(true, false);
	struct txVert {
		int material;
		int txId;
	}tv;
	std::unordered_multimap<int, txVert> multiMatVerts;  // index is position index and material of positions with multiple materials, second is texture index
	auto uniqueMtxAdd = [&](int& v, txVert& tx) {
		auto pr = multiMatVerts.equal_range(v);
		if (pr.first == multiMatVerts.end())
			multiMatVerts.insert(std::make_pair(v, tx));
		else {
			auto it = pr.first;
			while (it != pr.second) {
				if (it->second.material == tx.material)
					break;
				++it;
			}
			if (it == pr.second)
				multiMatVerts.insert(std::make_pair(v, tx));
		}
	};
	auto uniqueAddTx = [&](int& v, int& tx) {
		auto pr = _oneMaterialSeams.equal_range(v);
		if (pr.first == _oneMaterialSeams.end())
			_oneMaterialSeams.insert(std::make_pair(v, tx));
		else {
			auto it = pr.first;
			while (it != pr.second) {
				if (it->second == tx)
					break;
				++it;
			}
			if (it == pr.second)
				_oneMaterialSeams.insert(std::make_pair(v, tx));
		}
	};
	for (int n = numberOfTriangles(), i = 0; i < n; ++i) {
		int mati = triangleMaterial(i);
		unsigned int* adjs = triAdjs(i);
		int* tr = triangleVertices(i);
		int* tx = triangleTextures(i);
		for (int j = 0; j < 3; ++j) {
			int matj = triangleMaterial(adjs[j] >> 2);
			int* adjTx = triangleTextures(adjs[j] >> 2);
			if (mati > matj) {  // only need to do this edge once
				tv.material = mati;
				tv.txId = tx[j];
				uniqueMtxAdd(tr[j], tv);
				tv.material = matj;
				tv.txId = adjTx[((adjs[j] & 3) + 1) % 3];
				uniqueMtxAdd(tr[j], tv);
			}
			else if (mati == matj) {  // only need to do one vertex since other side will do the other
				int vtxA = adjTx[((adjs[j] & 3) + 1) % 3];
				if (tx[j] != vtxA) {
					uniqueAddTx(tr[j], tx[j]);
					uniqueAddTx(tr[j], vtxA);
				}
			}
			else
				;
		}
	}
	auto mvit = multiMatVerts.begin();
	while (mvit != multiMatVerts.end()) {
		auto start = mvit;
		auto source = multiMatVerts.end();
		while (mvit != multiMatVerts.end() && mvit->first == start->first) {
			if (mvit->second.material == 2)
				source = mvit;
			else if (mvit->second.material > 6 && source == multiMatVerts.end())
				source = mvit;
			else
				;
			++mvit;
		}
		assert(source != multiMatVerts.end());
		std::vector<neighborNode> nei;
		while (start != mvit) {
			if (start != source && start->second.txId == source->second.txId) {
				start->second.txId = addTexture();
				float tex[2];
				getTexture(source->second.txId, tex);
				setTexture(start->second.txId, tex);
				if (nei.empty())
					getNeighbors(source->first, nei);
				for (auto& n : nei) {
					if (n.triangle > -1 && triangleMaterial(n.triangle) == start->second.material) {
						int* trv = triangleVertices(n.triangle);
						int* ttx = triangleTextures(n.triangle);
						for (int j = 0; j < 3; ++j) {
							if (trv[j] == source->first) {
								ttx[j] = start->second.txId;
								break;
							}
						}
					}
				}
			}
			++start;
		}
	}
}

void materialTriangles::addOneMaterialTextureSeamVertex(int vertex, int(&textures)[2]) {

	if (vertex == 8151)
		std::cout << "In addTexSeam at 8151";
	std::vector<int> newTx;
	newTx.reserve(2);
	auto pr = _oneMaterialSeams.equal_range(vertex);
	if (pr.first == _oneMaterialSeams.end()) {
		for (int i = 0; i < 2; ++i)
			_oneMaterialSeams.insert(std::make_pair(vertex, textures[i]));
	}
	else {
		for (int i = 0; i < 2; ++i) {
			auto it = pr.first;
			while (it != pr.second) {
				if (it->second == textures[i])
					break;
				++it;
			}
			if (it == pr.second)
				newTx.push_back(textures[i]);
		}
	}
	for(auto &nt : newTx)
		_oneMaterialSeams.insert(std::make_pair(vertex, nt));
}

float materialTriangles::getDiameter() {
	findAdjacentTriangles(false, false);
	boundingBox<float> bb;
	bb.Empty_Box();
	Vec3f v, max;
	for (int n = (int)_vertexFace.size(), i = 0; i < n; ++i) {
		if (_vertexFace[i] == 0x80000000)
			continue;
		getVertexCoordinate(i, v._v);
		bb.Enlarge_To_Include_Point(v._v);
	}
	bb.Minimum_Corner(v._v);
	bb.Maximum_Corner(max._v);
	return (max - v).length();
}
