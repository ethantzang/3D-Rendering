/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

//--------------------------------------------
#ifndef GZ_EDGE
typedef	struct {	
  GzCoord	Spoint;
  GzCoord	Epoint;
  GzCoord	currentPoint;
  GzCoord	slope;
} GzEdge;
#define GZ_EDGE
#endif;

int RunScanLine(GzRender *render, int numParts, GzToken *nameList, GzPointer *valueList);
short	ctoi(float color);
//--------------------------------------------
#define PI  3.14159265
int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	if( mat == NULL ) return GZ_FAILURE;

	float radians = degree * (PI /180.0);
	
	mat[1][1] = cos(radians);
	mat[1][2] = -sin(radians);
	mat[2][1] = sin(radians);
	mat[2][2] = cos(radians);

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	if( mat == NULL ) return GZ_FAILURE;

	float radians = degree * (PI /180.0);
	
	mat[0][0] = cos(radians);
	mat[0][2] = sin(radians);
	mat[2][0] = -sin(radians);
	mat[2][2] = cos(radians);
		
	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	if( mat == NULL ) return GZ_FAILURE;

	float radians = degree * (PI /180.0);
	
	mat[0][0] = cos(radians);
	mat[0][1] = -sin(radians);
	mat[1][0] = sin(radians);
	mat[1][1] = cos(radians);
	
	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	if( mat == NULL ) return GZ_FAILURE;
	
	mat[0][3] = translate[X];
	mat[1][3] = translate[Y];
	mat[2][3] = translate[Z];

	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	if( mat == NULL ) return GZ_FAILURE;

	mat[0][0] = scale[X];
	mat[1][1] = scale[Y];
	mat[2][2] = scale[Z];
	
	return GZ_SUCCESS;
}
//----------------------------------------------------------
int NormalizeVector(GzCoord vec) 
{
	if (vec == NULL) {
		return GZ_FAILURE;
	}
	float norm = sqrt((vec[X]*vec[X]) + (vec[Y]*vec[Y]) + (vec[Z]*vec[Z]));
	vec[X] = vec[X] / norm;
	vec[Y] = vec[Y] / norm;
	vec[Z] = vec[Z] / norm;
	return GZ_SUCCESS;
}

int PutVector( GzCoord v1, GzCoord v2 )
{
	v1[X] = v2[X];
	v1[Y] = v2[Y];
	v1[Z] = v2[Z];
	return GZ_SUCCESS;
}

int SubtractVector( GzCoord v1, GzCoord v2, GzCoord Result )
{
	Result[X] = v1[X]-v2[X];
	Result[Y] = v1[Y]-v2[Y];
	Result[Z] = v1[Z]-v2[Z];
	return GZ_SUCCESS;
}

int AddVector( GzCoord v1, GzCoord v2, float additional ) 
{
	v1[X] = v1[X] + v2[X]*additional;
	v1[Y] = v1[Y] + v2[Y]*additional;
	v1[Z] = v1[Z] + v2[Z]*additional;
	return GZ_SUCCESS;
}

int MultiConstToVector( GzCoord result, GzCoord vec, float value )
{
	result[X] = vec[X]*value;
	result[Y] = vec[Y]*value;
	result[Z] = vec[Z]*value;
	return GZ_SUCCESS;
}

float DotProduct(GzCoord v1, GzCoord v2) 
{
	return v1[X]*v2[X] + v1[Y]*v2[Y] + v1[Z]*v2[Z];
}

int CrossProduct( GzCoord v1, GzCoord v2, GzCoord v3 )
{
	v3[X] = v1[Y]*v2[Z] - v1[Z]*v2[Y];
	v3[Y] = v1[Z]*v2[X] - v1[X]*v2[Z];
	v3[Z] = v1[X]*v2[Y] - v1[Y]*v2[X];
	return GZ_SUCCESS;
}

//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay *display)
{
/*  
- malloc a renderer struct 
- setup Xsp and anything only done once 
- save the pointer to display 
- init default camera 
*/ 
	//malloc 
	*render = (GzRender*)malloc(sizeof(GzRender));
	(*render)->display = display;
	(*render)->matlevel = -1;  //empty stack
	//setup Xsp matrix
	float xs = (*render)->display->xres;
	float ys = (*render)->display->yres;
	GzMatrix tmpXsp = { 
			xs/2,     0.0,     0.0,    xs/2,
			 0.0, -(ys/2),     0.0,    ys/2, 
			 0.0,     0.0, INT_MAX,     0.0, 
			 0.0,     0.0,     0.0,     1.0
	};
	memcpy( (*render)->Xsp, tmpXsp, sizeof(GzMatrix) );

	//init default camera
	(*render)->camera.FOV = DEFAULT_FOV;//35
	/* world coords for image plane origin */
	(*render)->camera.position[X] = DEFAULT_IM_X;
	(*render)->camera.position[Y] = DEFAULT_IM_Y;
	(*render)->camera.position[Z] = DEFAULT_IM_Z;
	/* default look-at point = 0,0,0 */
	(*render)->camera.lookat[X] = 0.0;
	(*render)->camera.lookat[Y] = 0.0;
	(*render)->camera.lookat[Z] = 0.0;
	/* default world-up point = 0,1,0 */
	(*render)->camera.worldup[X] = 0.0;
	(*render)->camera.worldup[Y] = 1.0;
	(*render)->camera.worldup[Z] = 0.0;

	(*render)->numlights = 0;

	return GZ_SUCCESS;

}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	free( render );
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render)
{
/*  
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 
	GzInitDisplay( render->display );
	/*--------------------Compute Xiw --------------------*/
	//Camera Z-axis Z = cl/|cl| compute cl and normalize cl
	GzCoord cameraZ;
	SubtractVector( render->camera.lookat, render->camera.position, cameraZ );
	NormalizeVector(cameraZ);
	//Camera Y-axis Y = up - (up.Z)*Z [up == world up]
	GzCoord cameraY;
	float dotProduct_UpZ = DotProduct( render->camera.worldup, cameraZ );
	MultiConstToVector( cameraY, cameraZ, dotProduct_UpZ ); //tmp cameraY = dotProduct_UpZ*cameraZ
	SubtractVector( render->camera.worldup, cameraY, cameraY ); //cameraT = worldup - tmp cameraY
	NormalizeVector(cameraY);
	//Camera X-axis X = Y X Z
	GzCoord cameraX;
	CrossProduct( cameraY, cameraZ, cameraX ); //X = Y X Z
	NormalizeVector(cameraX);
	//------------Building Matrix Xiw
	float dotProduct_XC = DotProduct( cameraX, render->camera.position );
	float dotProduct_YC = DotProduct( cameraY, render->camera.position );
	float dotProduct_ZC = DotProduct( cameraZ, render->camera.position );
	GzMatrix tmpXiw = {
		cameraX[X], cameraX[Y], cameraX[Z], -(dotProduct_XC),
		cameraY[X], cameraY[Y], cameraY[Z], -(dotProduct_YC),
		cameraZ[X], cameraZ[Y], cameraZ[Z], -(dotProduct_ZC),
		0.0,		   0.0,		   0.0,				1.0
	};

	memcpy( render->camera.Xiw, tmpXiw, sizeof(GzMatrix) );
	
	/*--------------------Compute Xpi --------------------*/
	// 1/d = tan(FOV/2)
	GzMatrix tmpXpi = {
		1,	0,	0,									0, 
		0,	1,	0,									0, 
		0,	0,	tan((render->camera.FOV/2)*PI/180),	0,
		0,	0,	tan((render->camera.FOV/2)*PI/180), 1
	};

	memcpy( render->camera.Xpi, tmpXpi, sizeof(GzMatrix) );

	/*---init Ximage - put Xsp at base of stack, push on Xpi and Xiw---*/
	//push Xsp into Ximage 
	GzPushMatrix( render, render->Xsp );
	//push Xpi into Ximage 
	GzPushMatrix( render, render->camera.Xpi );
	//push Xiw into Both Stack
	GzPushMatrix( render, render->camera.Xiw );

	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	render->camera.FOV = camera->FOV;
	PutVector( render->camera.lookat, camera->lookat );
	PutVector( render->camera.position, camera->position );
	PutVector( render->camera.worldup, camera->worldup );
	NormalizeVector(render->camera.worldup);

	return GZ_SUCCESS;	
}


int GzPushNormalsMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Xnorm stack
- stick with Ximage stack level
- matlevel =0 and =1 push I matrix into the stack
*/
	GzMatrix tmpTArray;
	GzMatrix Ident = { 
		1.0,	0.0,	0.0,	0.0, 
		0.0,	1.0,	0.0,	0.0, 
		0.0,	0.0,	1.0,	0.0, 
		0.0,	0.0,	0.0,	1.0 
	};
	//Do unitary rotation
	float K = 1/sqrt((matrix[0][0]*matrix[0][0] + matrix[0][1]*matrix[0][1] + matrix[0][2]*matrix[0][2]));
	for( int i=0; i<3; i++ ){// 3X3 R
		for( int j=0; j<3; j++ ){
			tmpTArray[i][j] = matrix[i][j]*K;
		}
	}
	//No translation 
	for( int i=0; i<3; i++ ){//UR 3 
		tmpTArray[i][3] = 0;
	}
	//Do multiplication with previous
	GzMatrix tmpNArray;
	int previous = render->matlevel - 1;
	if( previous >= 0 ){
		//render->Xnorm[render->matlevel]
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				tmpNArray[i][j]=0;
				for(int k=0; k<4; k++){
					tmpNArray[i][j]=tmpNArray[i][j]+render->Xnorm[previous][i][k]*tmpTArray[k][j];
				}
			}
		}
	}
	/*---init Xnorm -  put I & I at base of stack---*/
	if( render->matlevel == 0 || render->matlevel == 1 )//push I into stacks
		memcpy(render->Xnorm[render->matlevel], Ident, sizeof(GzMatrix));
	else if( render->matlevel > 1 )
		memcpy(render->Xnorm[render->matlevel], tmpNArray, sizeof(GzMatrix));
	
	return GZ_SUCCESS;
}


int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	//Do multiplication with previous
	float tmpArray[4][4];
	if(render->matlevel >= 0){
		//render->Ximage[render->matlevel]
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				tmpArray[i][j]=0;
				for(int k=0; k<4; k++){
					tmpArray[i][j]=tmpArray[i][j]+render->Ximage[render->matlevel][i][k]*matrix[k][j];
				}
			}
		}
	}

	render->matlevel++;

	//Overflow
	if(render->matlevel >= MATLEVELS){
		return GZ_FAILURE;
	}
	else{
		if( render->matlevel > 0 )
			memcpy(render->Ximage[render->matlevel], tmpArray, sizeof(GzMatrix));
		else if( render->matlevel == 0 )
			memcpy(render->Ximage[render->matlevel], matrix, sizeof(GzMatrix));
		GzPushNormalsMatrix( render, matrix );
	}
	
	return GZ_SUCCESS;
}


int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (render->matlevel < 0) {
		return GZ_FAILURE;
	}
	else{
		render->matlevel--;
	}

	return GZ_SUCCESS;
}


int GzShader(GzRender *render, GzColor Color, GzCoord N ) 
{
	//L,E,N,R
	//--L = render->lights[i].direction //--N = N that pass in
	//--Construct E ----------(0, 0, -1)----------//
	GzCoord E = { 0, 0, -1 };
	NormalizeVector(E);
	//N
	NormalizeVector(N);
	//--Construct R ----------R=2(N.L)N-L----------// 
	GzCoord* R = new GzCoord[render->numlights];
	int* ShadingCase = new int[render->numlights];
	float dotProduct_NL; 
	float dotProduct_NE = DotProduct( N, E );
	for( int i=0; i<render->numlights; i++ ){
		dotProduct_NL = DotProduct( N, render->lights[i].direction );
		if( dotProduct_NL >= 0 && dotProduct_NE >= 0 ){
			ShadingCase[i] = 1;
			R[i][X] = 2*dotProduct_NL*N[X] - render->lights[i].direction[X];
			R[i][Y] = 2*dotProduct_NL*N[Y] - render->lights[i].direction[Y];
			R[i][Z] = 2*dotProduct_NL*N[Z] - render->lights[i].direction[Z];
			NormalizeVector(R[i]);
		}
		else if( dotProduct_NL < 0 && dotProduct_NE < 0 ){
			ShadingCase[i] = -1;
			R[i][X] = 2*dotProduct_NL*(-N[X]) - render->lights[i].direction[X];
			R[i][Y] = 2*dotProduct_NL*(-N[Y]) - render->lights[i].direction[Y];
			R[i][Z] = 2*dotProduct_NL*(-N[Z]) - render->lights[i].direction[Z];
			NormalizeVector(R[i]);
		}
		else {
			ShadingCase[i] = 0;
		}
	}
	//Shading Equation--C=Specular+Diffuse+Ambient
	GzCoord tmpKs, tmpKd, tmpKa;
	GzCoord oneMatirx = {1.0, 1.0, 1.0};
	if( render->interp_mode == GZ_COLOR && render->tex_fun != NULL ){
		PutVector( tmpKs, oneMatirx );
		PutVector( tmpKd, oneMatirx );
		PutVector( tmpKa, oneMatirx );
	}
	else{
		PutVector( tmpKs, render->Ks );
		PutVector( tmpKd, render->Kd );
		PutVector( tmpKa, render->Ka );
	}
	
	//Compute Specular --Ks*sum(le(R.E)^spec)
	GzColor Specular = {0, 0, 0};
	
	for(int i=0; i<render->numlights; i++){
		if( ShadingCase[i] != 0 ){
			float dotProduct_RE = DotProduct( R[i], E );
			if( dotProduct_RE < 0 ) dotProduct_RE = 0;
			if( dotProduct_RE > 1 ) dotProduct_RE = 1;
			Specular[RED]	+= tmpKs[RED]   * render->lights[i].color[RED]	 * pow( dotProduct_RE, render->spec );
			Specular[GREEN] += tmpKs[GREEN] * render->lights[i].color[GREEN] * pow( dotProduct_RE, render->spec );
			Specular[BLUE]	+= tmpKs[BLUE]  * render->lights[i].color[BLUE]  * pow( dotProduct_RE, render->spec );
		}
	}
	//Compute Diffuse --Kd*sum(le(N.L))
	GzColor Diffuse = {0, 0, 0};
	dotProduct_NL = 0.0;
	for(int i=0; i<render->numlights; i++){
		if(ShadingCase[i] == 1){
			dotProduct_NL = DotProduct(N, render->lights[i].direction);
			Diffuse[0] += tmpKd[0] * render->lights[i].color[0] * dotProduct_NL;
			Diffuse[1] += tmpKd[1] * render->lights[i].color[1] * dotProduct_NL;
			Diffuse[2] += tmpKd[2] * render->lights[i].color[2] * dotProduct_NL;
		}
		else if(ShadingCase[i] == -1){
			GzCoord flipN = {-N[X], -N[Y], -N[Z]};
			dotProduct_NL = DotProduct(flipN, render->lights[i].direction);
			Diffuse[0] += tmpKd[0] * render->lights[i].color[0] * dotProduct_NL;
			Diffuse[1] += tmpKd[1] * render->lights[i].color[1] * dotProduct_NL;
			Diffuse[2] += tmpKd[2] * render->lights[i].color[2] * dotProduct_NL;
		}
	}
	//Compute Ambient --Ka*la
	GzColor Ambient = {0, 0, 0};
	Ambient[RED]   = tmpKa[RED]   * render->ambientlight.color[RED];
	Ambient[GREEN] = tmpKa[GREEN] * render->ambientlight.color[GREEN];
	Ambient[BLUE]  = tmpKa[BLUE]  * render->ambientlight.color[BLUE];
	//Complete the Equation
	
	Color[RED]	 = Specular[RED]   + Diffuse[RED]   + Ambient[RED];
	Color[GREEN] = Specular[GREEN] + Diffuse[GREEN] + Ambient[GREEN];
	Color[BLUE]  = Specular[BLUE]  + Diffuse[BLUE]  + Ambient[BLUE];

	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for(int i=0; i<numAttributes; i++){
		switch( nameList[i] ){
			case GZ_RGB_COLOR:
				PutVector( render->flatcolor, ((float*)valueList[i]) );
				break;
			case GZ_INTERPOLATE:
				render->interp_mode = *(int*)valueList[i];
				break;
			case GZ_DIRECTIONAL_LIGHT:
				if( render->numlights < MAX_LIGHTS ){
					render->lights[render->numlights] = *(GzLight*)valueList[i]; 
					render->numlights++;
				}
				break;
			case GZ_AMBIENT_LIGHT:
				render->ambientlight = *(GzLight*)valueList[i];//?
				break;
			case GZ_AMBIENT_COEFFICIENT:
				PutVector( render->Ka, ((float*)valueList[i]) );
				break;
			case GZ_DIFFUSE_COEFFICIENT:
				PutVector( render->Kd, ((float*)valueList[i]) );
				break;
			case GZ_SPECULAR_COEFFICIENT:
				PutVector( render->Ks, ((float*)valueList[i]) );
				break;
			case GZ_DISTRIBUTION_COEFFICIENT:
				render->spec = *(float*)valueList[i];
				break;
			case GZ_TEXTURE_MAP:
				render->tex_fun = (GzTexture)valueList[i];
				break;
			case GZ_AASHIFTX:
				render->offset[0] = *(float*)valueList[i];
				render->offset[2] = 0;
				break;
			case GZ_AASHIFTY:
				render->offset[1] = *(float*)valueList[i];
				break;
			default:
				break;
		}
	}//end of for 
	
	return GZ_SUCCESS;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts using matrix on top of stack 
- Clip - just discard any triangle with any vert(s) behind view plane 
       - optional: test for triangles with all three verts off-screen (trivial frustum cull)
- invoke triangle rasterizer  
*/ 
	GzCoord* triVertex = (GzCoord*) valueList[0];
	GzCoord* triNormal = (GzCoord*) valueList[1];
	GzTextureIndex* triTextUV = (GzTextureIndex*) valueList[2];
	GzCoord xformVertexList[3];
	GzCoord xformNormalList[3];
	GzCoord xformUVIndxList[3];
	GzPointer xformValueListTriangle[3];
	GzCoord R2;
	GzCoord surfNorm;
	GzCoord	xformSurfNorm;
	GzCoord R1;
	
	bool checkBehindView = false;
	for( int i=0; i<numParts; i++ ){
		if( nameList[i] == GZ_POSITION ) {
			//Compute Xform
			float w_xform = 0.0;
			for(int a=0; a<3; a++){//3 vertice
				//multiplication
				xformVertexList[a][X] = (render->Ximage[render->matlevel])[0][0]*triVertex[a][X]
									   +(render->Ximage[render->matlevel])[0][1]*triVertex[a][Y]
									   +(render->Ximage[render->matlevel])[0][2]*triVertex[a][Z]
									   +(render->Ximage[render->matlevel])[0][3]*1.0;
				xformVertexList[a][Y] = (render->Ximage[render->matlevel])[1][0]*triVertex[a][X]
									   +(render->Ximage[render->matlevel])[1][1]*triVertex[a][Y]
									   +(render->Ximage[render->matlevel])[1][2]*triVertex[a][Z]
									   +(render->Ximage[render->matlevel])[1][3]*1.0;
				xformVertexList[a][Z] = (render->Ximage[render->matlevel])[2][0]*triVertex[a][X]
									   +(render->Ximage[render->matlevel])[2][1]*triVertex[a][Y]
									   +(render->Ximage[render->matlevel])[2][2]*triVertex[a][Z]
									   +(render->Ximage[render->matlevel])[2][3]*1.0;
				w_xform = (render->Ximage[render->matlevel])[3][0]*triVertex[a][X]
						 +(render->Ximage[render->matlevel])[3][1]*triVertex[a][Y]
						 +(render->Ximage[render->matlevel])[3][2]*triVertex[a][Z]
						 +(render->Ximage[render->matlevel])[3][3]*1.0;
				//if( xformVertexList[a][Z] < render->camera.position[Z] ){
				if( xformVertexList[a][Z] < 0 ){
					checkBehindView = true;
					break;
				}
				MultiConstToVector( xformVertexList[a], xformVertexList[a], 1/w_xform  );
			}//end of for vertex	
		}
		if( nameList[i] == GZ_NORMAL ){
			//Compute Xform for normals
			for( int a=0; a<3; a++ ){// 3 Normals
				//multiplication 3x3 3x1->3x1 
				for( int i=0; i<3; i++ ){
					xformNormalList[a][i] = 0;
					for( int j=0; j<3; j++ ){
						xformNormalList[a][i] += (render->Xnorm[render->matlevel])[i][j]*triNormal[a][j];
					}
				}
				NormalizeVector(xformNormalList[a]);
			}
		}
		if( nameList[i] == GZ_TEXTURE_INDEX ){
			for( int a=0; a<3; a++ ){
				if( checkBehindView != true ){
					//not yet transform to the right perspective
					float tmpU = triTextUV[a][0];
					float tmpV = triTextUV[a][1];
					
					//if( tmpU > 1 ) tmpU = 1;//tmpU - floor(tmpU);
					//if( tmpU < 0 ) tmpU = 0;//1 - ceil(tmpU) + tmpU;
					//if( tmpV > 1 ) tmpV = 1;//tmpV - floor(tmpV);
					//if( tmpV < 0 ) tmpV = 0;//1 - ceil(tmpV) + tmpV;
					
					float Vz = xformVertexList[a][Z] / (INT_MAX - xformVertexList[a][Z]);
					xformUVIndxList[a][0] = tmpU / ( Vz+1.0 );
					xformUVIndxList[a][1] = tmpV / ( Vz+1.0 );
					xformUVIndxList[a][2] = 0;
				}
			}
		}
		if( nameList[i] == GZ_NULL_TOKEN ){
			return GZ_SUCCESS;
		}

	}//for i< numparts

	//Clipping
	if( checkBehindView != true ){
		//Compute Surface Normal
		if( render->interp_mode == GZ_FLAT ){
			SubtractVector( triVertex[1], triVertex[0], R1 );
			SubtractVector( triVertex[2], triVertex[1], R2 );
			CrossProduct( R1, R2, surfNorm );
			//Transform to image plane
			for( int i=0; i<3; i++ ){
				xformSurfNorm[i] = 0;
				for( int j=0; j<3; j++ ){
					xformSurfNorm[i] += (render->Xnorm[render->matlevel])[i][j]*surfNorm[j];
				}
			}
			NormalizeVector(xformSurfNorm);
			//GetColor
			GzColor Flat;
			GzShader( render, Flat, xformSurfNorm );
			PutVector( render->flatcolor, Flat );
		}

		xformValueListTriangle[0] = (GzPointer)xformVertexList;
		xformValueListTriangle[1] = (GzPointer)xformNormalList;
		xformValueListTriangle[2] = (GzPointer)xformUVIndxList;
		RunScanLine( render, 1, nameList, xformValueListTriangle );
	}
	return GZ_SUCCESS;
}

/* Run Scan Line */

int RunScanLine( GzRender *render, int numParts, GzToken *nameList, GzPointer *valueList  ){
	///====== Decide Min y and Max y vetex and Decide 1 2 3 ======///
	/////////////////////////////////////////////////////////////////
	float minY = ((GzCoord*)valueList[0])[0][1];
	float maxY = ((GzCoord*)valueList[0])[0][1];
	int topIndex=0, bottomIndex=0, midIndex=0;
	for( int i=0; i<3; i++ ){
		if( minY > ((GzCoord*)valueList[0])[i][1] ){
			minY = ((GzCoord*)valueList[0])[i][1];
			topIndex = i;
		}
		if( maxY < ((GzCoord*)valueList[0])[i][1] ){
			maxY = ((GzCoord*)valueList[0])[i][1];
			bottomIndex = i;
		}
		if( i == 2 ){
			for( int k=0; k<3; k++ ){
				if( k != topIndex && k != bottomIndex )
					midIndex = k;
			}
		}
	}
	///-- May at same Y -- look at their X value --///
	if( ((GzCoord*)valueList[0])[midIndex][1] == minY ){
		if( ((GzCoord*)valueList[0])[midIndex][0] > ((GzCoord*)valueList[0])[topIndex][0] ){
			midIndex--;
			topIndex++;
		}
	}
	else if( ((GzCoord*)valueList[0])[midIndex][1] == maxY ){
		if( ((GzCoord*)valueList[0])[midIndex][0] > ((GzCoord*)valueList[0])[bottomIndex][0] ){
			midIndex++;
			bottomIndex--;
		}
	}
	//////////////////////////////////////////////////////////////////////////////
	//----------------------- assign vertice norms UV --------------------------//
	GzCoord vertice[3]={0};
	GzCoord normals[3]={0};
	GzCoord texture[3]={0};
	float aX=0, aY=0, aZ=0;
	for( int i=0; i<3; i++ ){
		if( i == topIndex ){
			PutVector( vertice[0], ((GzCoord*)valueList[0])[i] );
			PutVector( normals[0], ((GzCoord*)valueList[1])[i] );
			PutVector( texture[0], ((GzCoord*)valueList[2])[i] );
		}
		else if( i == bottomIndex ){
			PutVector( vertice[2], ((GzCoord*)valueList[0])[i] );
			PutVector( normals[2], ((GzCoord*)valueList[1])[i] );
			PutVector( texture[2], ((GzCoord*)valueList[2])[i] );
		}
		else{
			PutVector( vertice[1], ((GzCoord*)valueList[0])[i] );
			PutVector( normals[1], ((GzCoord*)valueList[1])[i] );
			PutVector( texture[1], ((GzCoord*)valueList[2])[i] );
		}
	}

	for( int i=0; i<3; i++ ){
		vertice[i][X] = vertice[i][X] - render->offset[X];
		vertice[i][Y] = vertice[i][Y] - render->offset[Y];
	}
	
	//====== Setup Edges ======//
	GzEdge edges[3];
	// (1-2)
	PutVector( edges[0].Spoint, vertice[0] );
	PutVector( edges[0].Epoint, vertice[1] );
	// (2-3)
	PutVector( edges[1].Spoint, vertice[1] );
	PutVector( edges[1].Epoint, vertice[2] );
	// (1-3)
	PutVector( edges[2].Spoint, vertice[0] );
	PutVector( edges[2].Epoint, vertice[2] );
	// init the slope and currentPoint
	GzCoord tmpResult;
	for( int i=0; i<3; i++ ){
		PutVector( edges[i].currentPoint, edges[i].Spoint );
		//SubtractVector( edges[i].Epoint, edges[i].Spoint, tmpResult );
		//MultiConstToVector( edges[i].slope, tmpResult, 1 / (edges[i].Epoint[Y]-edges[i].Spoint[Y]) );
		edges[i].slope[X] = (edges[i].Epoint[X] - edges[i].Spoint[X])/
							(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(xe-xs)/(ye-ys) 
		edges[i].slope[Y] = (edges[i].Epoint[Y] - edges[i].Spoint[Y])/
							(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(ye-ys)/(ye-ys)
		edges[i].slope[Z] = (edges[i].Epoint[Z] - edges[i].Spoint[Z])/
							(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(ze-zs)/(ye-ys)
	}
	//====== Setup Normal Edges =====//
	GzEdge norms[3];
	// (N0-N1)
	PutVector( norms[0].Spoint, normals[0] );
	PutVector( norms[0].Epoint, normals[1] );
	// (N1-N2)
	PutVector( norms[1].Spoint, normals[1] );
	PutVector( norms[1].Epoint, normals[2] );
	// (N0-N2)
	PutVector( norms[2].Spoint, normals[0] );
	PutVector( norms[2].Epoint, normals[2] );
	// init 
	for( int i=0; i<3; i++ ){
		PutVector( norms[i].currentPoint, norms[i].Spoint );
		//SubtractVector( norms[i].Epoint, norms[i].Spoint, tmpResult );
		//MultiConstToVector( norms[i].slope, tmpResult, 1/(edges[i].Epoint[Y]-edges[i].Spoint[Y]) );
		norms[i].slope[X] = (norms[i].Epoint[X] - norms[i].Spoint[X])/
							(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(xe-xs)/(ye-ys) 
		norms[i].slope[Y] = (norms[i].Epoint[Y] - norms[i].Spoint[Y])/
							(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(ye-ys)/(ye-ys)
		norms[i].slope[Z] = (norms[i].Epoint[Z] - norms[i].Spoint[Z])/
							(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(ze-zs)/(ye-ys)

	}
	//--- Gourand Mode Get the Color and then interpolation ---//
	GzColor gourandColor[3];
	GzShader( render, gourandColor[0], normals[0] );
	GzShader( render, gourandColor[1], normals[1] );
	GzShader( render, gourandColor[2], normals[2] );
	//--- Setup Color Edges
	GzEdge colorEdge[3];
	// (C0-C1)
	PutVector( colorEdge[0].Spoint, gourandColor[0] );
	PutVector( colorEdge[0].Epoint, gourandColor[1] );
	// (C1-C2)
	PutVector( colorEdge[1].Spoint, gourandColor[1] );
	PutVector( colorEdge[1].Epoint, gourandColor[2] );
	// (C0-C2)
	PutVector( colorEdge[2].Spoint, gourandColor[0] );
	PutVector( colorEdge[2].Epoint, gourandColor[2] );
	// init
	for( int i=0; i<3; i++ ){
		PutVector( colorEdge[i].currentPoint, colorEdge[i].Spoint );
		//SubtractVector( colorEdge[i].Epoint, colorEdge[i].Spoint, tmpResult );
		//MultiConstToVector( colorEdge[i].slope, tmpResult, 1/(edges[i].Epoint[Y]-edges[i].Spoint[Y]) );
		colorEdge[i].slope[X] = (colorEdge[i].Epoint[X] - colorEdge[i].Spoint[X])/
												(edges[i].Epoint[Y] - edges[i].Spoint[Y]);
		colorEdge[i].slope[Y] = (colorEdge[i].Epoint[Y] - colorEdge[i].Spoint[Y])/
												(edges[i].Epoint[Y] - edges[i].Spoint[Y]);
		colorEdge[i].slope[Z] = (colorEdge[i].Epoint[Z] - colorEdge[i].Spoint[Z])/
												(edges[i].Epoint[Y] - edges[i].Spoint[Y]);

	}
	//====== Setup Texture Edges =====//
	GzEdge uvTexture[3];
	// (uv0-uv1)
	PutVector( uvTexture[0].Spoint, texture[0] );
	PutVector( uvTexture[0].Epoint, texture[1] );
	// (uv1-uv2)
	PutVector( uvTexture[1].Spoint, texture[1] );
	PutVector( uvTexture[1].Epoint, texture[2] );
	// (uv0-uv2)
	PutVector( uvTexture[2].Spoint, texture[0] );
	PutVector( uvTexture[2].Epoint, texture[2] );
	//init
	for( int i=0; i<3; i++ ){
		PutVector( uvTexture[i].currentPoint, uvTexture[i].Spoint );
		//SubtractVector( uvTexture[i].Epoint, uvTexture[i].Spoint, tmpResult );
		//MultiConstToVector( uvTexture[i].slope, tmpResult, 1/(edges[i].Epoint[Y]-edges[i].Spoint[Y]) );
		uvTexture[i].slope[X] = (uvTexture[i].Epoint[X] - uvTexture[i].Spoint[X])/
										(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(xe-xs)/(ye-ys) 
		uvTexture[i].slope[Y] = (uvTexture[i].Epoint[Y] - uvTexture[i].Spoint[Y])/
										(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(ye-ys)/(ye-ys)
		uvTexture[i].slope[Z] = (uvTexture[i].Epoint[Z] - uvTexture[i].Spoint[Z])/
										(edges[i].Epoint[Y] - edges[i].Spoint[Y]); //(ze-zs)/(ye-ys)
	}
	//====== Decide L or R triangle ======//
	int triCase = 2; // L case 1
	//(1-2) < (1-3) >> R >> case 2
	if( edges[0].slope[X] < edges[2].slope[X] )
		triCase = 0;
	if( vertice[0][1] == vertice[1][1] )
		triCase = 1;
	if( vertice[1][1] == vertice[2][1] )
		triCase = 3;
	//====== Advance (1-2) and (1-3) ======//
	float deltaY = 0; float deltaX=0;
	int edgeIndex=0, secondEdgeIndex=0; int spanIndex =0;
	bool checkOutOfBound;
	//GzSpan span; 
	GzEdge spanNorm; GzEdge spanColor; GzEdge spanEdge; GzEdge spanTexture;
	int pos_z_tmp=0;
	//====== Run ScanLine ======
	while( edges[2].currentPoint[1] < edges[2].Epoint[1] ){ //(1-3)
		// Move down the edges
		checkOutOfBound = false;
		if( triCase == 1 ){
			if( edgeIndex == 0 ){
				deltaY = ceil(vertice[0][1]) - vertice[0][1];
				//--- Edges
				AddVector( edges[2].currentPoint, edges[2].slope, deltaY );
				//--- Norms
				AddVector( norms[2].currentPoint, norms[2].slope, deltaY );
				//--- Color
				AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, deltaY );
				//--- UVTexture
				AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, deltaY );
				//-----------------------------------------------------------------//
				deltaY = ceil(vertice[1][1]) - vertice[1][1];
				//--- Edges
				AddVector( edges[1].currentPoint, edges[1].slope, deltaY );
				//--- Norms
				AddVector( norms[1].currentPoint, norms[1].slope, deltaY );
				//--- Color
				AddVector( colorEdge[1].currentPoint, colorEdge[1].slope, deltaY );
				//--- UVTexture
				AddVector( uvTexture[1].currentPoint, uvTexture[1].slope, deltaY );
			}
			else{
				//--- Edges
				AddVector( edges[2].currentPoint, edges[2].slope, 1 );
				AddVector( edges[1].currentPoint, edges[1].slope, 1 );
				//--- Norms
				AddVector( norms[2].currentPoint, norms[2].slope, 1 );
				AddVector( norms[1].currentPoint, norms[1].slope, 1 );
				//--- Color
				AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, 1 );
				AddVector( colorEdge[1].currentPoint, colorEdge[1].slope, 1 );
				//--- UVTexture
				AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, 1 );
				AddVector( uvTexture[1].currentPoint, uvTexture[1].slope, 1 );
			}
			edgeIndex++;
			//--- Edges Interpolation Set the Span---//
			PutVector( spanEdge.Epoint, edges[2].currentPoint );//End == Right
			PutVector( spanEdge.Spoint, edges[1].currentPoint );//Start == Left
			//--- Norms Interpolation Set the Span---//
			PutVector( spanNorm.Epoint, norms[2].currentPoint );//End == Right
			PutVector( spanNorm.Spoint, norms[1].currentPoint );//Start == Left
			//--- Color Interpolation Set the Span---//
			PutVector( spanColor.Epoint, colorEdge[2].currentPoint );
			PutVector( spanColor.Spoint, colorEdge[1].currentPoint );
			//--- UVtxt Interpolation Set the Span---//
			PutVector( spanTexture.Epoint, uvTexture[2].currentPoint );
			PutVector( spanTexture.Spoint, uvTexture[1].currentPoint );
		}
		else if( triCase == 3 ){
			if( edgeIndex == 0 ){
				deltaY = ceil(vertice[0][1]) - vertice[0][1];
				//--- Edges
				AddVector( edges[2].currentPoint, edges[2].slope, deltaY );
				AddVector( edges[0].currentPoint, edges[0].slope, deltaY );
				//--- Norms
				AddVector( norms[2].currentPoint, norms[2].slope, deltaY );
				AddVector( norms[0].currentPoint, norms[0].slope, deltaY );
				//--- Color
				AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, deltaY );
				AddVector( colorEdge[0].currentPoint, colorEdge[0].slope, deltaY );
				//--- UVTexture
				AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, deltaY );
				AddVector( uvTexture[0].currentPoint, uvTexture[0].slope, deltaY );
			}
			else{
				//--- Edges
				AddVector( edges[2].currentPoint, edges[2].slope, 1 );
				AddVector( edges[0].currentPoint, edges[0].slope, 1 );
				//--- Norms
				AddVector( norms[2].currentPoint, norms[2].slope, 1 );
				AddVector( norms[0].currentPoint, norms[0].slope, 1 );
				//--- Color
				AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, 1 );
				AddVector( colorEdge[0].currentPoint, colorEdge[0].slope, 1 );
				//--- UVTexture
				AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, 1 );
				AddVector( uvTexture[0].currentPoint, uvTexture[0].slope, 1 );
			}
			edgeIndex++;
			//--- Edges Interpolate set the Span ---//spanEdge
			PutVector( spanEdge.Epoint, edges[2].currentPoint );//Right == (1-3)
			PutVector( spanEdge.Spoint, edges[0].currentPoint );//Left == (1-2)
			//--- Norms Interpolate Set the Span ---//
			PutVector( spanNorm.Epoint, norms[2].currentPoint );
			PutVector( spanNorm.Spoint, norms[0].currentPoint );
			//--- Color Interpolate Set the Span ---//
			PutVector( spanColor.Epoint, colorEdge[2].currentPoint );
			PutVector( spanColor.Spoint, colorEdge[0].currentPoint );
			//--- UVtxt Interpolation Set the Span---//
			PutVector( spanTexture.Epoint, uvTexture[2].currentPoint );
			PutVector( spanTexture.Spoint, uvTexture[0].currentPoint );
		}// End of special Case 1 && 3 
		else{
			if( edges[2].currentPoint[Y] < vertice[1][Y] ){
				if( edgeIndex == 0 ){
					//--- First Round---
					deltaY = ceil(vertice[0][Y]) - vertice[0][Y]; //ceil(Y1)-Y1
					//--- Edges
					AddVector( edges[2].currentPoint, edges[2].slope, deltaY );
					AddVector( edges[0].currentPoint, edges[0].slope, deltaY );
					//--- Norms
					AddVector( norms[2].currentPoint, norms[2].slope, deltaY );
					AddVector( norms[0].currentPoint, norms[0].slope, deltaY );
					//--- Color
					AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, deltaY );
					AddVector( colorEdge[0].currentPoint, colorEdge[0].slope, deltaY );
					//--- UVTexture
					AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, deltaY );
					AddVector( uvTexture[0].currentPoint, uvTexture[0].slope, deltaY );
				}
				else {
					//--- Round 2 3 ... ---
					//--- Edges
					AddVector( edges[2].currentPoint, edges[2].slope, 1 );
					AddVector( edges[0].currentPoint, edges[0].slope, 1 );
					//--- Norms
					AddVector( norms[2].currentPoint, norms[2].slope, 1 );
					AddVector( norms[0].currentPoint, norms[0].slope, 1 );
					//--- Color
					AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, 1 );
					AddVector( colorEdge[0].currentPoint, colorEdge[0].slope, 1 );
					//--- UVTexture
					AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, 1 );
					AddVector( uvTexture[0].currentPoint, uvTexture[0].slope, 1 );
				}
				// Set up Span
				if( triCase == 2 ){ // left == (1-3)
					if( floor(edges[0].currentPoint[Y]) > render->display->yres ) break;
					if( ceil(edges[0].currentPoint[Y]) < 0 )	checkOutOfBound = true;
					//--- Edges 
					PutVector( spanEdge.Spoint, edges[2].currentPoint );//Lefe  (1-3)
					PutVector( spanEdge.Epoint, edges[0].currentPoint );//Right (1-2)
					//--- Norms
					PutVector( spanNorm.Spoint, norms[2].currentPoint );//Left = Start
					PutVector( spanNorm.Epoint, norms[0].currentPoint );//Right = End
					//--- Color
					PutVector( spanColor.Spoint, colorEdge[2].currentPoint );
					PutVector( spanColor.Epoint, colorEdge[0].currentPoint );
					//--- UVtxt 
					PutVector( spanTexture.Spoint, uvTexture[2].currentPoint );
					PutVector( spanTexture.Epoint, uvTexture[0].currentPoint );
				}
				else if( triCase == 0 ){ // right == (1-3)
					if( floor(edges[2].currentPoint[Y]) > render->display->yres ) break;
					if( ceil(edges[2].currentPoint[Y]) < 0 )	checkOutOfBound = true;
					//--- Edges
					PutVector( spanEdge.Epoint, edges[2].currentPoint );//Right (1-3)
					PutVector( spanEdge.Spoint, edges[0].currentPoint );//Right (1-3)
					//--- Norms
					PutVector( spanNorm.Epoint, norms[2].currentPoint );//Right = End
					PutVector( spanNorm.Spoint, norms[0].currentPoint );//Left = start
					//--- Color
					PutVector( spanColor.Epoint, colorEdge[2].currentPoint );
					PutVector( spanColor.Spoint, colorEdge[0].currentPoint );
					//--- UVtxt 
					PutVector( spanTexture.Epoint, uvTexture[2].currentPoint );
					PutVector( spanTexture.Spoint, uvTexture[0].currentPoint );
				}
			}
			edgeIndex++;
			//---------------- Switch ----------------//
			////////////////////////////////////////////
			if( edges[2].currentPoint[1] >= vertice[1][1] ) { 
				if( secondEdgeIndex == 0 ){
					deltaY = ceil(vertice[1][Y]) - (vertice[1][Y]);
					//--- Edges
					AddVector( edges[1].currentPoint, edges[1].slope, deltaY );
					//--- Norms
					AddVector( norms[1].currentPoint, norms[1].slope, deltaY );
					//--- Color
					AddVector( colorEdge[1].currentPoint, colorEdge[1].slope, deltaY );
					//--- UVTexture
					AddVector( uvTexture[1].currentPoint, uvTexture[1].slope, deltaY );
				}
				else{
					//--- Edges
					AddVector( edges[1].currentPoint, edges[1].slope, 1 );
					AddVector( edges[2].currentPoint, edges[2].slope, 1 );
					//--- Norm
					AddVector( norms[1].currentPoint, norms[1].slope, 1 );
					AddVector( norms[2].currentPoint, norms[2].slope, 1 );
					//--- Color
					AddVector( colorEdge[1].currentPoint, colorEdge[1].slope, 1 );
					AddVector( colorEdge[2].currentPoint, colorEdge[2].slope, 1 );
					//--- UVTexture
					AddVector( uvTexture[1].currentPoint, uvTexture[1].slope, 1 );
					AddVector( uvTexture[2].currentPoint, uvTexture[2].slope, 1 );
				}
				secondEdgeIndex++;
				//Set up Span
				if( triCase == 2 ){ //Left == (1-3)
					if( floor(edges[1].currentPoint[Y]) > render->display->yres ) break;
					if( ceil(edges[1].currentPoint[Y]) < 0 )	checkOutOfBound = true;
					//--- Edges
					PutVector( spanEdge.Spoint, edges[2].currentPoint );
					PutVector( spanEdge.Epoint, edges[1].currentPoint );
					//--- Norms
					PutVector( spanNorm.Spoint, norms[2].currentPoint );//Left == start
					PutVector( spanNorm.Epoint, norms[1].currentPoint );//Right == end
					//--- Color
					PutVector( spanColor.Spoint, colorEdge[2].currentPoint );
					PutVector( spanColor.Epoint, colorEdge[1].currentPoint );
					//--- UVtxt 
					PutVector( spanTexture.Spoint, uvTexture[2].currentPoint );
					PutVector( spanTexture.Epoint, uvTexture[1].currentPoint );
				}
				else if( triCase == 0 ){ //Right == (1-3)
					if( floor(edges[2].currentPoint[Y]) > render->display->yres ) break;
					if( ceil(edges[2].currentPoint[Y]) < 0 )	checkOutOfBound = true;
					//--- Edges
					PutVector( spanEdge.Epoint, edges[2].currentPoint );
					PutVector( spanEdge.Spoint, edges[1].currentPoint );
					//--- Norms
					PutVector( spanNorm.Epoint, norms[2].currentPoint );
					PutVector( spanNorm.Spoint, norms[1].currentPoint );
					//--- Color
					PutVector( spanColor.Epoint, colorEdge[2].currentPoint );
					PutVector( spanColor.Spoint, colorEdge[1].currentPoint );
					//--- UVtxt 
					PutVector( spanTexture.Epoint, uvTexture[2].currentPoint );
					PutVector( spanTexture.Spoint, uvTexture[1].currentPoint );
				}
			}
			/////////////////////////////////////////////////////////
			//---------------- End of Switch-----------------------//
		}// end else [triCase]
		//===============================================================//
		//SubtractVector( spanEdge.Epoint, spanEdge.Spoint, tmpResult );
		//MultiConstToVector( spanEdge.slope, tmpResult, 1.0/( spanEdge.Epoint[0]-spanEdge.Spoint[0] ) );
		spanEdge.slope[X] = (spanEdge.Epoint[X] -spanEdge.Spoint[X])/(spanEdge.Epoint[X] - spanEdge.Spoint[X]);
		spanEdge.slope[Y] = (spanEdge.Epoint[Y] -spanEdge.Spoint[Y])/(spanEdge.Epoint[X] - spanEdge.Spoint[X]);
		spanEdge.slope[Z] = (spanEdge.Epoint[Z] -spanEdge.Spoint[Z])/(spanEdge.Epoint[X] - spanEdge.Spoint[X]);
		PutVector( spanEdge.currentPoint, spanEdge.Spoint );
		//--------------Normals------------------------
		SubtractVector( spanNorm.Epoint, spanNorm.Spoint, tmpResult );
		MultiConstToVector( spanNorm.slope, tmpResult, 1.0/( spanEdge.Epoint[0]-spanEdge.Spoint[0] ) );
		PutVector( spanNorm.currentPoint, spanNorm.Spoint);
		//--------------Colors------------------------
		SubtractVector( spanColor.Epoint, spanColor.Spoint, tmpResult );
		MultiConstToVector( spanColor.slope, tmpResult, 1.0/( spanEdge.Epoint[0]-spanEdge.Spoint[0] ) );
		PutVector( spanColor.currentPoint, spanColor.Spoint );
		//--------------Texture------------------------
		SubtractVector( spanTexture.Epoint, spanTexture.Spoint, tmpResult );
		MultiConstToVector( spanTexture.slope, tmpResult, 1.0/( spanEdge.Epoint[0]-spanEdge.Spoint[0] ) );
		PutVector( spanTexture.currentPoint, spanTexture.Spoint );
		//--------------------------------------
		// Run Span
		deltaX = ceil(spanEdge.Spoint[0]) - spanEdge.Spoint[0];
		spanIndex=0;
		while( spanEdge.currentPoint[0] < spanEdge.Epoint[0] ){
			if( spanEdge.Epoint[0] < 0 ) break;
			if( checkOutOfBound ) break;
			if( floor(spanEdge.currentPoint[0]) > render->display->xres ) break;
			if( spanIndex == 0 ){
				//--- Edges
				AddVector( spanEdge.currentPoint, spanEdge.slope, deltaX);
				//--- Norms
				AddVector( spanNorm.currentPoint, spanNorm.slope, deltaX);
				//--- Color
				AddVector( spanColor.currentPoint, spanColor.slope, deltaX );
				//--- UVTxt
				AddVector( spanTexture.currentPoint, spanTexture.slope, deltaX );
			}
			else{
				//--- Edges
				AddVector( spanEdge.currentPoint, spanEdge.slope, 1);
				//--- Norms
				AddVector( spanNorm.currentPoint, spanNorm.slope, 1);
				//--- Color
				AddVector( spanColor.currentPoint, spanColor.slope, 1 );
				//--- UVTxt
				AddVector( spanTexture.currentPoint, spanTexture.slope, 1 );
			}
			if( spanEdge.currentPoint[0] < spanEdge.Epoint[0] && spanEdge.currentPoint[0] >= 0){
				// Write color valuw in FB
				// Compare z value
				GzColor result={0};
				GzColor textureColor={0};
				GzTextureIndex UV={0};
				pos_z_tmp = (edges[2].currentPoint[1]*render->display->xres + spanEdge.currentPoint[0]);
				if( spanEdge.currentPoint[2] < (render->display->fbuf[pos_z_tmp].z) ){
					float V_z = spanEdge.currentPoint[2] / (INT_MAX - spanEdge.currentPoint[2]);
					UV[0] = spanTexture.currentPoint[0] * (V_z + 1.0);
					UV[1] = spanTexture.currentPoint[1] * (V_z + 1.0);
					switch( render->interp_mode ){
						case GZ_COLOR:
							if( render->tex_fun != NULL ){
								render->tex_fun(UV[0], UV[1], textureColor);
								result[0] = textureColor[0] * spanColor.currentPoint[0];
								result[1] = textureColor[1] * spanColor.currentPoint[1];
								result[2] = textureColor[2] * spanColor.currentPoint[2];
							}
							else
								PutVector( result, spanColor.currentPoint );
							break;
						case GZ_NORMALS:
							if (render->tex_fun != NULL) {	
								render->tex_fun(UV[0], UV[1], textureColor); 
								PutVector( render->Ka, textureColor );
								PutVector( render->Kd, textureColor );
							}
							GzShader( render, result, spanNorm.currentPoint);
							break;
						case GZ_FLAT:
							PutVector(result, render->flatcolor);
							break;
						default:
							break;
					}

					GzPutDisplay( render->display, spanEdge.currentPoint[0], edges[2].currentPoint[1], 
									(GzIntensity)ctoi(result[RED]), 
									(GzIntensity)ctoi(result[GREEN]), 
									(GzIntensity)ctoi(result[BLUE]), 
									1, 
									(GzDepth)spanEdge.currentPoint[2]);
				}
			}
			spanIndex++;
		}//end while span
		
	}//end while edge

	return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

