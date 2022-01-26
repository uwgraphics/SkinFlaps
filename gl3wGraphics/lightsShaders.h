// file: lightsShaders.h
// Author: CourtCutting, MD
// Date: January 31,2012
// Purpose: This class does basic lighting and coordinated shaders.  It expects
//	vertex attributes vVertex(vec4), vNormal(vec3), and vTexture(vec2) as inputs.
//	It also expects uniforms mvpMatrix(mat4), mvMatrix(mat4), and mNormal(mat3)
//	which are the modelview-projection, modelview, and inverse modelview rotation
//	matrices respectively. If a 2Dtexture is associated with the object, that is
//	also handled here. It is expected in the future that this class will be jazzed
//	up significantly by a shader guru, but it at least provides a basic level.

#ifndef __lightsShaders_h__
#define __lightsShaders_h__

#include <vector>
#include <map>
#include <string>
#include <GL/gl3w.h>

// forward declarations
class GLmatrices;

class lightsShaders
{
public:
	GLuint getTextureProgramNumber() {return _textureProgram;}
	GLuint getLineProgramNumber() {return _lineProgram;}
	void useGlslProgram(GLuint programNumber);  // careful - no error checking for validity
	GLuint getOrCreateLineProgram();
	GLuint getOrCreateColorProgram();
	void setColor(GLfloat *color);
	void setMaterial(int material);  // added 3/7/15
	void createTextureProgram();
	bool createCustomProgram(GLuint &program, const char *vertexShader, const char *fragmentShader, std::vector<std::string> &attributes);
	bool createNormalTangentProgram(); // creates geometry shader program for computation of normals and tangents
	bool loadCustomProgram(GLuint &program, const char *vertexShaderFile, const char *fragmentShaderFile, std::vector<std::string> &attributes);
	void setModelMatrix(GLfloat *model);
	void setProgramUniforms(GLuint program);
	void setGLmatrices(GLmatrices *glM) {_glM=glM;}
	void clear();    // clears graphics card of programs and uniforms
	lightsShaders();
	~lightsShaders();

private:
	static bool createProgramWithAttributes(GLuint &program, const char *vertexShader, const char *fragmentShader, std::vector<std::string> &attributes);
	GLmatrices *_glM;
	static GLuint _textureProgram;
	static GLuint _colorProgram;
	static GLuint _lineProgram;
	static GLuint _normalTangentProgram;
//	GLuint _textureBufferObjects[2],_texBOBuffers[2];
	GLuint _currentProgram;
	GLfloat _objectColor[4];
	GLint _material;
	GLfloat _vEyeLight[3];
	GLfloat _vAmbientColor[4];
	GLfloat _vDiffuseColor[4];
	GLfloat _vSpecularColor[4];
	GLfloat _MVP[16], _modelMat[16], _normMat[9];
	struct progUniforms {
		bool	notDoneOnce;		// Next lines should only be loaded once
		GLint	locAmbient;			// The location of the ambient color
		GLint   locDiffuse;			// The location of the diffuse color
		GLint   locSpecular;		// The location of the specular color
		GLint	locLight;			// The location of the Light in eye coordinates
		GLint   locTexture0;		// Textured object with one texture
		GLint   locTexture1;		// Second texture
		GLint   locTexture2;		// Third texture
		GLint   locTexture3;		// Forth texture
		// next lines are frequently changed
		GLint	locMVP;				// The location of the ModelViewProjection matrix uniform
		GLint	locMV;				// The location of the ModelView matrix uniform
		GLint	locPM;				// The location of the Projection matrix uniform
		GLint	locNM;				// The location of the Normal matrix uniform
		GLint	locObjColor;		// Only for colored, non-textured objects
		GLint	locMaterial;		// Material ID specifying which texturing algorithm to apply. Added 3/7/2015.
	}_texUni, _lineUni, _colorUni;
	std::map<GLuint, progUniforms> _programUniforms;  // uniform locations associated with each custom program
//	struct lineUniforms{
//		GLint	locMVP;				// The location of the ModelViewProjection matrix uniform
//		GLint	locMV;				// The location of the ModelView matrix uniform
//		GLint	color;				// Only for colored, non-textured objects
//	}_lineUni;

/*    const char *__GTVertexShaderColoredLine;
    const char *__GTFragmentShaderColoredLine;
    const char *__GTVertexShaderColoredPhong;
    const char *__GTFragmentShaderColoredPhong;
    const char *__GTVertexShaderDefault;
    const char *__GTFragmentShaderDefault;
    std::vector< const char* > cleanup_shaders; */
    
};


#endif	// __lightsShaders_h__
