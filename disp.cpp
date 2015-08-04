/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"


int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
 -- NOTE: this function is optional and not part of the API, but you may want to use it within the display function.
*/
	//======= Boundary Check =======//
	if (framebuffer == NULL) {
		return GZ_FAILURE;
	}
	if (width <= 0 || height <= 0 || width > MAXXRES || height > MAXYRES) {
		return GZ_FAILURE;
	}
	*framebuffer = (char*)malloc( width * height * sizeof(GzPixel) );  
	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	//======= Boundary Check =======//
	if (xRes <= 0 || yRes <= 0 || xRes > MAXXRES || yRes > MAXYRES) {
		return GZ_FAILURE;
	}
	
	*display = (GzDisplay*)malloc(sizeof(GzDisplay));
	(*display)->xres = xRes;
	(*display)->yres = yRes;
	//GzNewFrameBuffer( (char*)(*display)->fbuf, xRes, yRes );
	(*display)->fbuf = (GzPixel*)malloc(xRes*xRes*sizeof(GzPixel));

	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* clean up, free memory */
	free( (display->fbuf) );
	free( display );
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* pass back values for a display */
	*xRes = display->xres;
	*yRes = display->yres;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */
	
	for(int i=0; i<(display->xres)*(display->yres); i++){
		display->fbuf[i].blue  = 2000;//3500;//3080;//2500;
		display->fbuf[i].green = 2000;//3500;//3080;//1500;
		display->fbuf[i].red   = 2000;//3500;//3080;//1500;
		display->fbuf[i].z = INT_MAX;
	}
	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */
	//============= Boundary Check =============//
	if( a >= 4095 ) a = 4095;
	if( b >= 4095 ) b = 4095;
	if( g >= 4095 ) g = 4095;
	if( r >= 4095 ) r = 4095;
	if( z >= INT_MAX ) z = INT_MAX;
	if( a <= 0 ) a = 0;
	if( b <= 0 ) b = 0;
	if( g <= 0 ) g = 0;
	if( r <= 0 ) r = 0;
	if( z <= 0 ) z = 0;
	if( i > display->xres ) return GZ_FAILURE;
	if( j > display->yres ) return GZ_FAILURE;
	if( i < 0 ) return GZ_FAILURE;
	if( j < 0 ) return GZ_FAILURE;
	//==========================================//
	display->fbuf[j*(display->xres)+i].alpha = a;
	display->fbuf[j*(display->xres)+i].blue = b;
	display->fbuf[j*(display->xres)+i].green = g;
	display->fbuf[j*(display->xres)+i].red = r;
	display->fbuf[j*(display->xres)+i].z = z;
	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* pass back pixel value in the display */
	*a = display->fbuf[j*(display->xres)+i].alpha;
	*b = display->fbuf[j*(display->xres)+i].blue;
	*g = display->fbuf[j*(display->xres)+i].green;
	*r = display->fbuf[j*(display->xres)+i].red;
	*z = display->fbuf[j*(display->xres)+i].z;
	
	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
/* write pixels to ppm file -- "P6 %d %d 255\r" */
	fprintf(outfile, "P6 %d %d 255\r", (display->xres), (display->yres));
	for(int i=0; i<(display->xres)*(display->yres);i++){
		fprintf(outfile, "%c%c%c", (char)((display->fbuf[i].red)>>4), (char)((display->fbuf[i].green)>>4), (char)((display->fbuf[i].blue)>>4));	
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{
/* write pixels to framebuffer: 
	- Put the pixels into the frame buffer
	- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
	- Not red, green, and blue !!!
*/
	
	for(int i=0; i<(display->xres)*(display->yres); i++){
		framebuffer[3*i] = (display->fbuf[i].blue)>>4;
		framebuffer[3*i+1] = (display->fbuf[i].green)>>4;
		framebuffer[3*i+2] = (display->fbuf[i].red)>>4;
	}

	return GZ_SUCCESS;
}