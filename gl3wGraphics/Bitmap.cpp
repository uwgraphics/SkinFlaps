//File: Bitmap.cpp
//Written by:     Mark Bernard
//on GameDev.net: Captain Jester
//e-mail: mark.bernard@rogers.com
//Please feel free to use and abuse this code as much
//as you like.  But, please give me some credit for
//starting you off on the right track.
//
//The file Bitmap.h goes along with this file
//
#include "Bitmap.h"

//basic constructor
Bitmap::Bitmap(){
    reset();
}

//destructor
Bitmap::~Bitmap(){
    if(colours!=0) {
        delete[] colours;
    }
    if(data!=0) {
        delete[] data;
    }
}

//load a bitmap from a file and represent it correctly
//in memory
bool Bitmap::loadBMP(const char *file) {
    FILE *in;                  //file stream for reading
    char *tempData;       //temp storage for image data
    int numColours;            //total available colours

    //bitmap is not loaded yet
    loaded=false;
    //make sure memory is not lost
    if(colours!=0) {
        delete[] colours;
    }
    if(data!=0) {
        delete[] data;
    }

    //open the file for reading in binary mode
    in=fopen(file,"rb");

    //if the file does not exist return in error
    if(in==NULL) {
        error="File not found";
        return false;
    }

    //read in the entire BITMAPFILEHEADER
    fread(&bmfh,sizeof(BitmapFileHeader),1,in);
	cout << "sizeof(BitmapFileHeader)=" << sizeof(BitmapFileHeader) << endl;
    //check for the magic number that says this is a bitmap
    if(bmfh.bfType!=BITMAP_MAGIC_NUMBER) {
        error="File is not in DIB format";
        fclose(in);
        return false;
    }

    //read in the entire BITMAPINFOHEADER
    fread(&bmih,sizeof(BitmapInfoHeader),1,in);
	cout << "sizeof(BitmapInfoHeader)=" << sizeof(BitmapInfoHeader) << endl;

    //save the width, height and bits per pixel for external use
    width=bmih.biWidth;
    height=bmih.biHeight;
    bpp=bmih.biBitCount;
	cout << "biBitCount      =" << bmih.biBitCount << endl;
	cout << "biClrImportant  =" << bmih.biClrImportant << endl;
	cout << "biClrUsed       =" << bmih.biClrUsed << endl;
	cout << "biCompression   =" << bmih.biCompression << endl;
	cout << "biHeight        =" << bmih.biHeight << endl;
	cout << "biPlanes        =" << bmih.biPlanes << endl;
	cout << "biSize          =" << bmih.biSize << endl;
	cout << "biSizeImage     =" << bmih.biSizeImage << endl;
	cout << "biWidth         =" << bmih.biWidth << endl;
	cout << "biXPelsPerMeter =" << bmih.biXPelsPerMeter << endl;
	cout << "biYPelsPerMeter =" << bmih.biYPelsPerMeter << endl;

    //calculate the size of the image data with padding
    dataSize=(width*height*(unsigned int)(bmih.biBitCount/8.0));

    //calculate the number of available colours
    numColours=1<<bmih.biBitCount;

    //if the bitmap is not 8 bits per pixel or more
    //return in error
    if(bpp<8) {
        error="File is not 8 or 24 bits per pixel";
        fclose(in);
        return false;
    }

    //load the palette for 8 bits per pixel
    if(bpp==8) {
    	colours=new RGBQuad[numColours];
    	fread(colours,sizeof(RGBQuad),numColours,in);
    }

    //set up the temporary buffer for the image data
    tempData=new char[dataSize];

    //exit if there is not enough memory
    if(tempData==NULL) {
        error="Not enough memory to allocate a temporary buffer";
        fclose(in);
        return false;
    }

    //read in the entire image
    fread(tempData,sizeof(char),dataSize,in);

    //close the file now that we have all the info
    fclose(in);

    //calculate the witdh of the final image in bytes
    byteWidth=padWidth=(int)((float)width*(float)bpp/8.0);

    //adjust the width for padding as necessary
    while(padWidth%4!=0) {
        padWidth++;
    }

    //change format from GBR to RGB
    if(bpp==8) {
    	loaded=convert8(tempData);
   	}
    else if(bpp==24) {
    	loaded=convert24(tempData);
   	}

    //clean up memory
    delete[] tempData;

    //bitmap is now loaded
    error="Bitmap loaded";

    //return success
    return loaded;
}

//function to set the inital values
void Bitmap::reset(void) {
    loaded=false;
    colours=0;
    data=0;
    error="";
}

bool Bitmap::convert24(char* tempData) {
	int offset,diff;

	diff=width*height*RGB_BYTE_SIZE;
    //allocate the buffer for the final image data
    data=new char[diff];

    //exit if there is not enough memory
    if(data==NULL) {
        error="Not enough memory to allocate an image buffer";
        delete[] data;
        return false;
    }

    if(height>0) {
        offset=padWidth-byteWidth;
        //count backwards so you start at the front of the image
        for(unsigned int i=0;i<dataSize;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+i+2)=*(tempData+i);
            *(data+i+1)=*(tempData+i+1);
            *(data+i)=*(tempData+i+2);
        }
    }

    //image parser for a forward image
    else {
        offset=padWidth-byteWidth;
        int j=dataSize-3;
        //count backwards so you start at the front of the image
		//here you can start from the back of the file or the front,
		//after the header  The only problem is that some programs
		//will pad not only the data, but also the file size to
		//be divisible by 4 bytes.
        for(unsigned int i=0;i<dataSize;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+j+2)=*(tempData+i);
            *(data+j+1)=*(tempData+i+1);
            *(data+j)=*(tempData+i+2);
            j-=3;
        }
    }

    return true;
}

bool Bitmap::convert8(char* tempData) {
	int offset,diff;

	diff=width*height*RGB_BYTE_SIZE;
    //allocate the buffer for the final image data
    data=new char[diff];

    //exit if there is not enough memory
    if(data==NULL) {
        error="Not enough memory to allocate an image buffer";
        delete[] data;
        return false;
    }

    if(height>0) {
        offset=padWidth-byteWidth;
        int j=0;
        //count backwards so you start at the front of the image
        for(unsigned int i=0;i<dataSize*RGB_BYTE_SIZE;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+i)=colours[*(tempData+j)].rgbRed;
            *(data+i+1)=colours[*(tempData+j)].rgbGreen;
            *(data+i+2)=colours[*(tempData+j)].rgbBlue;
            j++;
        }
    }

    //image parser for a forward image
    else {
        offset=padWidth-byteWidth;
        int j=dataSize-1;
        //count backwards so you start at the front of the image
        for(unsigned int i=0;i<dataSize*RGB_BYTE_SIZE;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+i)=colours[*(tempData+j)].rgbRed;
            *(data+i+1)=colours[*(tempData+j)].rgbGreen;
            *(data+i+2)=colours[*(tempData+j)].rgbBlue;
            j--;
        }
    }

    return true;
}
