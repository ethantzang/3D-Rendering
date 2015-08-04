/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	"math.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
    }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */
	//Bounding check
	if (u < 0) u = 0;//1 - ceil(u) + u;//u = 0;
	if (u > 1) u = 1;//u - floor(u);//1;
	if (v < 0) v = 0;//1 - ceil(v) + v;//0;
	if (v > 1) v = 1;//v - floor(v);//1;
	//Scale to image texture size
	float imgU = u*(xs-1);
	float imgV = v*(ys-1);
	//Determine ABCD coord
	int tmpX = floor(imgU);
	int tmpY = floor(imgV);
	int A = (tmpY * xs) + tmpX;
	tmpX =  ceil(imgU);
	tmpY = floor(imgV);
	int B = (tmpY * xs) + tmpX;
	tmpX = ceil(imgU);
	tmpY = ceil(imgV);
	int C = (tmpY * xs) + tmpX;
	tmpX = floor(imgU);
	tmpY =  ceil(imgV);
	int D = (tmpY * xs) + tmpX;
	//
	float s = imgU - floor(imgU);
	float t = imgV - floor(imgV);
	// Interpolate
	color[RED]	 = s*t*image[C][RED]   + (1-s)*t*image[D][RED]   + s*(1-t)*image[B][RED]   + (1-s)*(1-t)*image[A][RED];
	color[GREEN] = s*t*image[C][GREEN] + (1-s)*t*image[D][GREEN] + s*(1-t)*image[B][GREEN] + (1-s)*(1-t)*image[A][GREEN];
	color[BLUE]  = s*t*image[C][BLUE]  + (1-s)*t*image[D][BLUE]  + s*(1-t)*image[B][BLUE]  + (1-s)*(1-t)*image[A][BLUE];
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	/*Checkerboard
	//GzColor color1, color2;
	int size = 8;
	int stride = 2;
	int index = 0;
	//Bounding Check
	if (u < 0) u = 0;//1 - ceil(u) + u;//u = 0;
	if (u > 1) u = 1;//u - floor(u);//1;
	if (v < 0) v = 0;//1 - ceil(v) + v;//0;
	if (v > 1) v = 1;//v - floor(v);//1;
	//Scale to image texture size
	float checkerBoardU = u*(size-1);
	float checkerBoardV = v*(size-1);
	int tmpU = checkerBoardU/stride;
	int tmpV = checkerBoardV/stride;
	if( (tmpU+tmpV)%2 == 0 ){
		color[RED]   = 200;
		color[GREEN] = 200;
		color[BLUE]  = 200;
	}
	else{
		color[RED]   = 50;
		color[GREEN] = 50;
		color[BLUE]  = 50;
	}*/
	/*Checkerboard
	if (u < 0) u = 0;//1 - ceil(u) + u;//u = 0;
	if (u > 1) u = 1;//u - floor(u);//1;
	if (v < 0) v = 0;//1 - ceil(v) + v;//0;
	if (v > 1) v = 1;//v - floor(v);//1;
	bool x=true,y=true;
	u = u*100;
	u = (int)u%10;
	if(u>=0 && u<5) x=true;
	else x = false;
	v= v*100;
	v = (int)v%10;
	if(v>=0 && v<5) y=true;
	else y = false;

	if(x^y){
		color[RED] = 100.0;
		color[GREEN] = 0;
		color[BLUE] = 0;
	}else{
		color[RED] = 0;
		color[GREEN] = 0;
		color[BLUE] = 255.0;
	}*/
	float xc = -0.5;
	float yc = -0.5;
	float zoom = 5.0;
	int maxIteration = 30;

	float xmin = xc - 0.5*zoom;
	float ymin = yc - 0.5*zoom;
	float x = xmin + zoom*u;
	float y = ymin + zoom*v;

	float new_a = x;
	float new_b = y;
	float old_a, old_b, ab,h;
	int i;
	for( i=0; i < maxIteration; i++ ){
		old_a = new_a * new_a;
		old_b = new_b * new_b;
		ab = 2 * new_a * new_b;
		if ((old_a + old_b) > 255) break;
		new_a = old_a - old_b + x;
		new_b = ab + y;
	}

	float remap = i - log(log(old_a+old_b)/2.0) / log(2.0);
	remap = 1 - remap;

	color[RED] = remap;
	color[GREEN] = remap;
	color[BLUE] = remap;

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

