// textures.cpp
// Author: Court Cutting with help from Richard Wright and others
// Date: January 22, 2012
// Purpose: Texture and bump map loader

#include <fstream>
#ifdef linux
#include <stdlib.h>
#endif
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "stb_image.h"  // above defines needed

#include "Bitmap.h"
#include "textures.h"

textures::textures()
{
	_textures.clear();
}

textures::~textures() {
    clear();
}

void textures::clear() {	// clears textures from graphics card
    std::map<int,tex>::iterator tit;
    for(tit=_textures.begin(); tit!=_textures.end(); ++tit)	{
        glDeleteTextures(1, &(tit->second.texture));
    }
    _textures.clear();
}

int textures::textureExists(std::string &textureName)
{	// 0 return is no texture found
	std::map<int,tex>::iterator tit;
	for(tit=_textures.begin(); tit!=_textures.end(); ++tit)	{
		if(tit->second.name==textureName)
			return tit->first;
	}
	return 0L;
}

int textures::loadTexture(int txId, const char *fileName)	{
	int existTx = textureExists(std::string(fileName));
	if (existTx > 0) {
		if (existTx == txId)
			return txId;
		else
			return 0;
	}
	std::pair<std::map<int, tex>::iterator,bool> pr = _textures.insert(std::make_pair(txId,tex()));
	if (!pr.second)
		return 0;
	pr.first->second.name = fileName;
	glGenTextures(1, &(pr.first->second.texture));
	glBindTexture(GL_TEXTURE_2D, pr.first->second.texture);
	bool ret;
	if (pr.first->second.name.size() - pr.first->second.name.rfind(".bmp") == 4)
		ret = loadBMPTexture(pr.first->second.name.c_str(), GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_TEXTURE_WRAP_S, pr.first->second.width, pr.first->second.height);
    else if (pr.first->second.name.size() - pr.first->second.name.rfind(".jpg") == 4)
        ret = loadJpgTexture(pr.first->second.name.c_str(), pr.first->second.width, pr.first->second.height);
    else if (pr.first->second.name.size() - pr.first->second.name.rfind(".tga") == 4)
		ret = LoadTGATexture(pr.first->second.name.c_str(), GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_TEXTURE_WRAP_S, pr.first->second.width, pr.first->second.height);
	else
		ret = false;
	if(!ret)	{
		glDeleteTextures(1, &(pr.first->second.texture));
		_textures.erase(pr.first);
		return 0xffffffff;
	}
	else
		return pr.first->first;
}

bool textures::getTextureSize(const int txId, int &width, int &height) {
    auto tid = _textures.find(txId);
    if (tid == _textures.end())
        return false;
    width = tid->second.width;
    height = tid->second.height;
    return true;
}

bool textures::loadJpgTexture(const char* fileName, int &width, int &height) {
    int channels;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* img = stbi_load(fileName, &width, &height, &channels, 0);
    if (img == NULL) {
        printf("Error in loading the image\n");
        return false;
    }
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_S);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_TEXTURE_WRAP_S);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, channels, width, height, 0,
        GL_RGB, GL_UNSIGNED_BYTE, img);
    glGenerateMipmap(GL_TEXTURE_2D);
    stbi_image_free(img);
    return true;
}

bool textures::loadBMPTexture(const char *fileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode, int& width, int& height)
{
	Bitmap bmp;
	// Read the texture bits
	bmp.loadBMP(fileName);
	if(bmp.data == NULL) 
		return false;
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
    
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, bmp.width, bmp.height, 0,
		GL_RGB, GL_UNSIGNED_BYTE, bmp.data);
    delete[] bmp.data ;
	bmp.data = NULL;
    if(minFilter == GL_LINEAR_MIPMAP_LINEAR || 
       minFilter == GL_LINEAR_MIPMAP_NEAREST ||
       minFilter == GL_NEAREST_MIPMAP_LINEAR ||
       minFilter == GL_NEAREST_MIPMAP_NEAREST)
        glGenerateMipmap(GL_TEXTURE_2D);
	return true;
}

// Load a TGA as a 2D Texture. Completely initialize the state
bool textures::LoadTGATexture(const char *szFileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode, int& width, int& height)
{
	GLbyte *pBits;
	int nWidth, nHeight, nComponents;
	GLenum eFormat;
	
	// Read the texture bits
	pBits = gltReadTGABits(szFileName, &nWidth, &nHeight, &nComponents, &eFormat);
	if(pBits == NULL) 
		return false;
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
    
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, nComponents, nWidth, nHeight, 0,
				 eFormat, GL_UNSIGNED_BYTE, pBits);
	
    free(pBits);
    
    if(minFilter == GL_LINEAR_MIPMAP_LINEAR || 
       minFilter == GL_LINEAR_MIPMAP_NEAREST ||
       minFilter == GL_NEAREST_MIPMAP_LINEAR ||
       minFilter == GL_NEAREST_MIPMAP_NEAREST)
        glGenerateMipmap(GL_TEXTURE_2D);
    
	return true;
}

////////////////////////////////////////////////////////////////////
// Allocate memory and load targa bits. Returns pointer to new buffer,
// height, and width of texture, and the OpenGL format of data.
// Call free() on buffer when finished!
// This only works on pretty vanilla targas... 8, 24, or 32 bit color
// only, no palettes, no RLE encoding.
GLbyte* textures::gltReadTGABits(const char *szFileName, GLint *iWidth, GLint *iHeight, GLint *iComponents, GLenum *eFormat)
{
    FILE *pFile;			// File pointer
    TGAHEADER tgaHeader;		// TGA file header
    unsigned int lImageSize;		// Size in bytes of image
    short sDepth;			// Pixel depth;
    GLbyte	*pBits = NULL;          // Pointer to bits
    // Default/Failed values
    *iWidth = 0;
    *iHeight = 0;
    *eFormat = GL_RGB;
    *iComponents = GL_RGB;
    // Attempt to open the file
    pFile = fopen(szFileName, "rb");
    if(pFile == NULL)
        return NULL;
    // Read in header (binary)
    fread(&tgaHeader, 18/* sizeof(TGAHEADER)*/, 1, pFile);
    // Do byte swap for big vs little endian
#ifdef __APPLE__
    LITTLE_ENDIAN_WORD(&tgaHeader.colorMapStart);
    LITTLE_ENDIAN_WORD(&tgaHeader.colorMapLength);
    LITTLE_ENDIAN_WORD(&tgaHeader.xstart);
    LITTLE_ENDIAN_WORD(&tgaHeader.ystart);
    LITTLE_ENDIAN_WORD(&tgaHeader.width);
    LITTLE_ENDIAN_WORD(&tgaHeader.height);
#endif
    // Get width, height, and depth of texture
    *iWidth = tgaHeader.width;
    *iHeight = tgaHeader.height;
    sDepth = tgaHeader.bits / 8;
    // Put some validity checks here. Very simply, I only understand
    // or care about 8, 24, or 32 bit targa's.
    if(tgaHeader.bits != 8 && tgaHeader.bits != 24 && tgaHeader.bits != 32)
        return NULL;
    // Calculate size of image buffer
    lImageSize = tgaHeader.width * tgaHeader.height * sDepth;
    // Allocate memory and check for success
    pBits = (GLbyte*)malloc(lImageSize * sizeof(GLbyte));
    if(pBits == NULL)
        return NULL;
    // Read in the bits
    // Check for read error. This should catch RLE or other 
    // weird formats that I don't want to recognize
    if(fread(pBits, lImageSize, 1, pFile) != 1)
		{
        free(pBits);
        return NULL;
		}
    // Set OpenGL format expected
    switch(sDepth)
		{
#ifndef OPENGL_ES
        case 3:     // Most likely case
            *eFormat = GL_BGR;
            *iComponents = GL_RGB;
            break;
#endif
        case 4:
            *eFormat = GL_BGRA;
            *iComponents = GL_RGBA;
            break;
        case 1:
//            *eFormat = GL_LUMINANCE;
 //           *iComponents = GL_LUMINANCE;
            break;
        default:        // RGB
            // If on the iPhone, TGA's are BGR, and the iPhone does not 
            // support BGR without alpha, but it does support RGB,
            // so a simple swizzle of the red and blue bytes will suffice.
            // For faster iPhone loads however, save your TGA's with an Alpha!
        break;
		}
    // Done with File
    fclose(pFile);
    // Return pointer to image data
    return pBits;
}

