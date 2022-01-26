// textures.h
// Author: Court Cutting with help from Richard Wright and others
// Date: January 22, 2012
// Purpose: Texture and bump map loader
#ifndef __textures_h__
#define __textures_h__

#include <map>
#include <GL/gl3w.h>			// OpenGL Extension "autoloader"

class textures
{
public:
	int loadTexture(int txId, const char *fileName);	// 0 means load failure, otherwise is txId
	bool textureExists(int txId) { return _textures.find(txId) != _textures.end(); }
	int textureExists(std::string &textureName); // 0 return is no texture found, otherwise returns txId
	GLuint getOGLtextureNumber(int txId) { return _textures[txId].texture; }
	bool getTextureSize(const int txId, int& width, int& height);
	void clear();	// clears textures from graphics card
	textures();
	~textures();

private:
	struct tex{
		std::string name;
		GLuint texture;
		int width;
		int height;
	};
	std::map<int,tex> _textures;
	bool LoadTGATexture(const char *szFileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode, int& width, int& height);
	bool loadBMPTexture(const char *fileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode, int& width, int& height);
	bool loadJpgTexture(const char* fileName, int& width, int& height);
	GLbyte* gltReadTGABits(const char *szFileName, GLint *iWidth, GLint *iHeight, GLint *iComponents, GLenum *eFormat);
	// Define targa header. This is only used locally.
	#pragma pack(push)
	#pragma pack(1)
	typedef struct
	{
		GLbyte	identsize;              // Size of ID field that follows header (0)
		GLbyte	colorMapType;           // 0 = None, 1 = paletted
		GLbyte	imageType;              // 0 = none, 1 = indexed, 2 = rgb, 3 = grey, +8=rle
		unsigned short	colorMapStart;          // First colour map entry
		unsigned short	colorMapLength;         // Number of colors
		unsigned char 	colorMapBits;   // bits per palette entry
		unsigned short	xstart;                 // image x origin
		unsigned short	ystart;                 // image y origin
		unsigned short	width;                  // width in pixels
		unsigned short	height;                 // height in pixels
		GLbyte	bits;                   // bits per pixel (8 16, 24, 32)
		GLbyte	descriptor;             // image descriptor
	} TGAHEADER;
	#pragma pack(8)
	struct RGB { 
		GLbyte blue;
		GLbyte green;
		GLbyte red;
		GLbyte alpha;
	};
	#pragma pack(pop)
};
#endif	// __textures_h__
