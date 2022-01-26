#ifndef _BITMAP_H
#define _BITMAP_H
//File: Bitmap.h
//Written by:     Mark Bernard
//on GameDev.net: Captain Jester
//e-mail: mark.bernard@rogers.com
//Please feel free to use and abuse this code as much
//as you like.  But, please give me some credit for
//starting you off on the right track.
//
//The file Bitmap.cpp goes along with this file
//
#include <iostream>
#include <cstdio>
#include <string>
using namespace std;

const short BITMAP_MAGIC_NUMBER=19778;
const int RGB_BYTE_SIZE=3;

#pragma pack(push,bitmap_data,1)

typedef struct tagRGBQuad {
	char rgbBlue;
	char rgbGreen;
	char rgbRed;
	char rgbReserved;
} RGBQuad;

typedef struct tagBitmapFileHeader {
	unsigned short bfType;
	unsigned int bfSize;
	unsigned short bfReserved1;
	unsigned short bfReserved2;
	unsigned int bfOffBits;
} BitmapFileHeader;

typedef struct tagBitmapInfoHeader {
	unsigned int biSize;
	int biWidth;
	int biHeight;
	unsigned short biPlanes;
	unsigned short biBitCount;
	unsigned int biCompression;
	unsigned int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	unsigned int biClrUsed;
	unsigned int biClrImportant;
} BitmapInfoHeader;

#pragma pack(pop,bitmap_data)

class Bitmap {
public:
    //variables
    RGBQuad *colours;
    char *data;
    bool loaded;
    int width,height;
    unsigned short bpp;
    string error;
    //methods
    Bitmap(void);
    ~Bitmap();
    bool loadBMP(const char *);
private:
    //variables
    BitmapFileHeader bmfh;
    BitmapInfoHeader bmih;
    int byteWidth;            //the width in bytes of the image
    int padWidth;             //the width in bytes of the added image
    unsigned int dataSize;                //size of the data in the file
    //methods
    void reset(void);
    bool convert24(char *);		//convert to 24bit RGB bottom up data
    bool convert8(char *);		//convert to 24bit RGB bottom up data
};

#endif //_BITMAP_H

