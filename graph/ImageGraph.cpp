#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "FreeImage.h"
#include "MyFloatType.h"

#define BPP 24

// #define DIAGS

// Cones
// int const StartD = 14;
// int const D = 50;
// char FileLeft[] = "ConesL_S14_D50.png";
// char FileRight[] = "ConesR_S14_D50.png";

// Kopi1
// int const StartD = 12+5; // =17
// int const D = 23+12+5+5; // =45
// char FileLeft[] = "Kopi1L_S17_D45.png";
// char FileRight[] = "Kopi1R_S17_D45.png";

// Kopi2
// int const StartD = 25+5; // =30
// int const D = 45+25+5+5; // =80
// char FileLeft[] = "Kopi2L_S30_D80.png";
// char FileRight[] = "Kopi2R_S30_D80.png";

// Tsukuba
// int const StartD = 17;
// int const D = 17;
// char FileLeft[] = "tsukubaL_S17_D17.png";
// char FileRight[] = "tsukubaR_S17_D17.png";

// Leaves
int const StartD = 15;
int const D = 30;
char FileLeft[] = "leaves2L.png";
char FileRight[] = "leaves2R.png";

FloatType CostScale = 1.0/10.0;
FloatType const MinValue = 1e-3;
FloatType const Gamma1 = 0.05; // continuity enforcement, smaller => more continuous
FloatType const Gamma2 = 0.05; // continuity enforcement, smaller => more continuous
#ifdef DIAGS
FloatType const Gamma3 = 0.7; // continuity enforcement, smaller => more continuous
FloatType const Gamma4 = 0.7; // continuity enforcement, smaller => more continuous
#endif

using namespace std;

/** Generic image loader
	@param lpszPathName Pointer to the full file name
	@param flag Optional load flag constant
	@return Returns the loaded dib if successful, returns NULL otherwise
*/
FIBITMAP* GenericLoader(const char* lpszPathName, int flag);

struct MyRGB
{
	FloatType R;
	FloatType G;
	FloatType B;
};

namespace GenerateUtils
{
	void MyPut(int value)
	{
		fwrite(&(value), sizeof(value), 1, stdout);
	}
	void MyPut(FloatType value)
	{
		fwrite(&(value), sizeof(value), 1, stdout);
	}
	void ComputeDataCosts(
		vector< vector< MyRGB > > const & ImageLeft ,
		vector< vector< MyRGB > > const & ImageRight ,
		vector< vector< vector< FloatType > > > & DataCosts);
};
using namespace GenerateUtils;







//  __       __            __              ___  ___
// |  \     /  \          |  \            /   \|   \
// | $$\   /  $$  ______   \$$ _______   /  $$$ \$$$\
// | $$$\ /  $$$ |      \ |  \|       \ |  $$     \$$\
// | $$$$\  $$$$  \$$$$$$\| $$| $$$$$$$\| $$      | $$
// | $$\$$ $$ $$ /      $$| $$| $$  | $$| $$      | $$
// | $$ \$$$| $$|  $$$$$$$| $$| $$  | $$ \$$\_  _/  $$
// | $$  \$ | $$ \$$    $$| $$| $$  | $$  \$$ \|   $$
//  \$$      \$$  \$$$$$$$ \$$ \$$   \$$   \$$$ \$$$




int main ( int argc , char ** argv )
{
	FreeImage_Initialise();
	// LEFT means Left Eye
	FIBITMAP * bitmap1 = GenericLoader(FileLeft,0);
	if( ! bitmap1 ) cerr << "no bitmap1" << endl;
	// RIGHT means Right Eye
	FIBITMAP * bitmap2 = GenericLoader(FileRight,0);
	if( ! bitmap2 ) cerr << "no bitmap2" << endl;
	RGBQUAD color ;

	int width = FreeImage_GetWidth(bitmap1);
	int height = FreeImage_GetHeight(bitmap1);
	{
		int width2 = FreeImage_GetWidth(bitmap2);
		int height2 = FreeImage_GetHeight(bitmap2);
		if (width != width2 || height != height2)
		{
			cerr << "error: different size images." << endl ;
			exit(1);
		}
	}

	vector< vector< MyRGB > > ImageLeft;
	vector< vector< MyRGB > > ImageRight;

	// (0,0) is bottom left of image
	ImageLeft.resize(width);
	ImageRight.resize(width);
	for(int i=0;i<width;i++)
	{
		ImageLeft[i].resize(height);
		ImageRight[i].resize(height);
		for(int j=0;j<height;j++)
		{
			FreeImage_GetPixelColor( bitmap1 , i , j , &color );
			ImageLeft[i][j].R = color.rgbRed;
			ImageLeft[i][j].G = color.rgbGreen;
			ImageLeft[i][j].B = color.rgbBlue;
			FreeImage_GetPixelColor( bitmap2 , i , j , &color );
			ImageRight[i][j].R = color.rgbRed;
			ImageRight[i][j].G = color.rgbGreen;
			ImageRight[i][j].B = color.rgbBlue;
		}
	}

	FreeImage_Unload(bitmap1);
	FreeImage_Unload(bitmap2);
	FreeImage_DeInitialise();

	// Do some truncation:
	if(StartD>0)
		ImageLeft.erase(ImageLeft.begin(),ImageLeft.begin()+StartD);
	else
		ImageRight.erase(ImageRight.begin(),ImageRight.begin()-StartD);

	vector< vector< vector< FloatType > > > DataCosts;
	ComputeDataCosts(ImageLeft,ImageRight,DataCosts);



	width = DataCosts.size();
	cerr << "width = " << width << endl;
	cerr << "height = " << height << endl;


	// number of variable nodes
	MyPut(width*height);

	// cardinality of each variable
	for(int i=0;i<width*height;i++)
	{
		MyPut(D);
	}

	// number of function nodes
	#ifndef DIAGS
	MyPut( (int)( height*(width-1) + width*(height-1) + width*height ) );
	#else
	MyPut( (int)( 2*(height-1)*(width-1) + height*(width-1) + width*(height-1) + width*height ) );
	#endif

	#ifdef DIAGS
	// variable nodes are enumerated left to right, then down
	for(int i=0;i<width-1;i++)
	{
		for(int j=0;j<height-1;j++)
		{
			// function node degree
			MyPut((int)2);
			// connected variable nodes
			MyPut((int)(i + j*width));
			MyPut((int)(i+1 + (j+1)*width));
			// which values vector
			MyPut((int)(width*height+1));
		}
	}

	for(int i=0;i<width-1;i++)
	{
		for(int j=0;j<height-1;j++)
		{
			// function node degree
			MyPut((int)2);
			// connected variable nodes
			MyPut((int)(i+1 + j*width));
			MyPut((int)(i + (j+1)*width));
			// which values vector
			MyPut((int)(width*height+1));
		}
	}
	#endif

	for(int i=0;i<width-1;i++)
	{
		for(int j=0;j<height;j++)
		{
			// function node degree
			MyPut((int)2);
			// connected variable nodes
			MyPut((int)(i + j*width));
			MyPut((int)(i+1 + j*width));
			// which values vector
			MyPut((int)(width*height));
		}
	}

	for(int i=0;i<height-1;i++)
	{
		for(int j=0;j<width;j++)
		{
			// function node degree
			MyPut((int)2);
			// connected variable nodes
			MyPut((int)(j + i*width));
			MyPut((int)(j+width + i*width));
			// which values vector
			MyPut((int)(width*height));
		}
	}

	for(int i=0;i<width*height;i++)
	{
		// function node degree
		MyPut((int)1);
		// connected variable nodes
		MyPut((int)i);
		// which values vector
		MyPut((int)i);
	}






	// how many values vectors
	#ifdef DIAGS
	MyPut(width*height+2);
	#else
	MyPut(width*height+1);
	#endif

	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			MyPut(D);
			// function kernel values
			for(int k=0;k<D;k++)
			{
				MyPut((FloatType)( DataCosts[j][i][k] ) ) ;
			}
		}
	}

	MyPut(D*D);
	// function kernel values
	for(int j=0;j<D;j++)
	{
		for(int k=0;k<D;k++)
		{
			if(j==k) MyPut((FloatType)1.0) ;
			else if(abs(j-k)==1) MyPut((FloatType)Gamma2) ;
			else MyPut((FloatType)Gamma1);
		}
	}

	#ifdef DIAGS
	MyPut(D*D);
	// function kernel values
	for(int j=0;j<D;j++)
	{
		for(int k=0;k<D;k++)
		{
			if(j==k) MyPut((FloatType)1.0) ;
			else if(abs(j-k)==1) MyPut((FloatType)Gamma4) ;
			else MyPut((FloatType)Gamma3);
		}
	}
	#endif

	MyPut(width);
	MyPut(height);
	MyPut(D);
	return 0;
}




//   ______                                                 __                _______              __                 ______                         __                  ___  ___
//  /      \                                               |  \              |       \            |  \               /      \                       |  \                /   \|   \
// |  $$$$$$\  ______   ______ ____    ______   __    __  _| $$_     ______  | $$$$$$$\  ______  _| $$_     ______  |  $$$$$$\  ______    _______  _| $$_     _______  /  $$$ \$$$\
// | $$   \$$ /      \ |      \    \  /      \ |  \  |  \|   $$ \   /      \ | $$  | $$ |      \|   $$ \   |      \ | $$   \$$ /      \  /       \|   $$ \   /       \|  $$     \$$\
// | $$      |  $$$$$$\| $$$$$$\$$$$\|  $$$$$$\| $$  | $$ \$$$$$$  |  $$$$$$\| $$  | $$  \$$$$$$\\$$$$$$    \$$$$$$\| $$      |  $$$$$$\|  $$$$$$$ \$$$$$$  |  $$$$$$$| $$      | $$
// | $$   __ | $$  | $$| $$ | $$ | $$| $$  | $$| $$  | $$  | $$ __ | $$    $$| $$  | $$ /      $$ | $$ __  /      $$| $$   __ | $$  | $$ \$$    \   | $$ __  \$$    \ | $$      | $$
// | $$__/  \| $$__/ $$| $$ | $$ | $$| $$__/ $$| $$__/ $$  | $$|  \| $$$$$$$$| $$__/ $$|  $$$$$$$ | $$|  \|  $$$$$$$| $$__/  \| $$__/ $$ _\$$$$$$\  | $$|  \ _\$$$$$$\ \$$\_  _/  $$
//  \$$    $$ \$$    $$| $$ | $$ | $$| $$    $$ \$$    $$   \$$  $$ \$$     \| $$    $$ \$$    $$  \$$  $$ \$$    $$ \$$    $$ \$$    $$|       $$   \$$  $$|       $$  \$$ \|   $$
//   \$$$$$$   \$$$$$$  \$$  \$$  \$$| $$$$$$$   \$$$$$$     \$$$$   \$$$$$$$ \$$$$$$$   \$$$$$$$   \$$$$   \$$$$$$$  \$$$$$$   \$$$$$$  \$$$$$$$     \$$$$  \$$$$$$$    \$$$ \$$$
//                                   | $$
//                                   | $$
//                                    \$$


namespace GenerateUtils
{
	// Ignore Left image columns left of StartD
	// pixel (StartD,0) of LeftImg should correspond with a pixel
	//     between RightImg (0,0) and (D-1,0)
	// 
	// To choose D and StartD - pixel (X+StartD,Y) from LeftImg should correspond
	//     with px in RightImg between (X+0,Y) and (X+D,Y)
	//   - Visually find biggest shift to the right and to the left, relative to LeftImg
	//   - left relative shift should be closer, right relative shift farther objects
	//   - choose StartD = LeftShift + SmallExcess (could be negative if LeftShift<0)
	//   - choose D = LeftShift + RightShift + SmallExcess (should certainly be positive)
	void ComputeDataCosts(
		vector< vector< MyRGB > > const & ImageLeft ,
		vector< vector< MyRGB > > const & ImageRight ,
		vector< vector< vector< FloatType > > > & DataCosts)
	{
		DataCosts.resize( min( ImageRight.size()-D+1 , ImageLeft.size() ) );
		for(int i=0;i<DataCosts.size();i++)
		{
			DataCosts[i].resize(ImageLeft[0].size());
			for(int j=0;j<DataCosts[i].size();j++)
			{
				DataCosts[i][j].assign(D,0.0);
				for(int k=0;k<D;k++)
				{
					MyRGB const & PixL = ImageLeft [i][j];
					MyRGB const & PixR = ImageRight[  i + k   ][j];
					DataCosts[i][j][k] = exp( -CostScale * ( abs(PixR.R-PixL.R)+abs(PixR.G-PixL.G)+abs(PixR.B-PixL.B) ) );
				}
			}
		}
		for(int i=0;i<DataCosts.size();i++)
		{
			for(int j=0;j<DataCosts[i].size();j++)
			{
				FloatType MaxVal = DataCosts[i][j][0];
				for(int k=1;k<D;k++)
				{
					if(DataCosts[i][j][k] > MaxVal) MaxVal = DataCosts[i][j][k];
				}
				for(int k=1;k<D;k++)
				{
					DataCosts[i][j][k] = DataCosts[i][j][k] / MaxVal;
					if(DataCosts[i][j][k] < MinValue) DataCosts[i][j][k] = MinValue;
				}
			}
		}
	}
};





//   ______                                           __            __                                 __                        ___  ___
//  /      \                                         |  \          |  \                               |  \                      /   \|   \
// |  $$$$$$\  ______   _______    ______    ______   \$$  _______ | $$       ______    ______    ____| $$  ______    ______   /  $$$ \$$$\
// | $$ __\$$ /      \ |       \  /      \  /      \ |  \ /       \| $$      /      \  |      \  /      $$ /      \  /      \ |  $$     \$$\
// | $$|    \|  $$$$$$\| $$$$$$$\|  $$$$$$\|  $$$$$$\| $$|  $$$$$$$| $$     |  $$$$$$\  \$$$$$$\|  $$$$$$$|  $$$$$$\|  $$$$$$\| $$      | $$
// | $$ \$$$$| $$    $$| $$  | $$| $$    $$| $$   \$$| $$| $$      | $$     | $$  | $$ /      $$| $$  | $$| $$    $$| $$   \$$| $$      | $$
// | $$__| $$| $$$$$$$$| $$  | $$| $$$$$$$$| $$      | $$| $$_____ | $$_____| $$__/ $$|  $$$$$$$| $$__| $$| $$$$$$$$| $$       \$$\_  _/  $$
//  \$$    $$ \$$     \| $$  | $$ \$$     \| $$      | $$ \$$     \| $$     \\$$    $$ \$$    $$ \$$    $$ \$$     \| $$        \$$ \|   $$
//   \$$$$$$   \$$$$$$$ \$$   \$$  \$$$$$$$ \$$       \$$  \$$$$$$$ \$$$$$$$$ \$$$$$$   \$$$$$$$  \$$$$$$$  \$$$$$$$ \$$         \$$$ \$$$




/** Generic image loader
	@param lpszPathName Pointer to the full file name
	@param flag Optional load flag constant
	@return Returns the loaded dib if successful, returns NULL otherwise
*/
FIBITMAP* GenericLoader(const char* lpszPathName, int flag) {
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;

	// check the file signature and deduce its format
	// (the second argument is currently not used by FreeImage)
	fif = FreeImage_GetFileType(lpszPathName, 0);
	if(fif == FIF_UNKNOWN) {
		// no signature ?
		// try to guess the file format from the file extension
		fif = FreeImage_GetFIFFromFilename(lpszPathName);
	}
	// check that the plugin has reading capabilities ...
	if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
		// ok, let's load the file
		FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
		// unless a bad file format, we are done !
		return dib;
	}
	return NULL;
}







