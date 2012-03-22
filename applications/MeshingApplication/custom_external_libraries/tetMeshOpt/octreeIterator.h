/*

MATH3D Functions
 */

#include <iostream>
#include <algorithm>
#include <math.h>
#include <windows.h>




typedef float float32[32];

typedef unsigned int uint;

 
 
////////////////////////////////////////////////
//---- Script
/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
 
 /*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 * 
 * Tridiagonal solvers.
 * Device code for parallel cyclic reduction (PCR).
 *
 * Original CUDA kernels: UC Davis, Yao Zhang & John Owens, 2009
 * 
 * NVIDIA, Nikolai Sakharnykh, 2009
 */

uint rgbaFloatToInt(float4 rgba)
{
   /* rgba.x = clamp(rgba.x,0.0f,1.0f);  
    rgba.y = clamp(rgba.y,0.0f,1.0f);  
    rgba.z = clamp(rgba.z,0.0f,1.0f);  
    rgba.w = clamp(rgba.w,0.0f,1.0f);  */
    return ((uint)(255.0f)<<24) | ((uint)(rgba.z*255.0f)<<16) | ((uint)(rgba.y*255.0f)<<8) | (uint)(rgba.x*255.0f);
}


int pop(int stack[],int *indexCell)
    {
      (*indexCell)-- ;
      return stack[(*indexCell)];
    }

    void push(int n,int stack[],int *indexCell)
    {
       if ((*indexCell)>=100) return;
	   stack[*indexCell] = n;
       (*indexCell)++;
    }

	void sort(int nodesIndexes[],int nc, float* nodes,float4 rayPos)
	{
         for (int i=0;i<nc;i++)
			for (int j=i+1;j<nc;j++)
			{
				 GridCell n1 = packToCell(nodes,32*nodesIndexes[i]);
				 GridCell n2 = packToCell(nodes,32*nodesIndexes[j]);
				 float d1 =distance(n1.center,rayPos);
                 float d2 =distance(n2.center,rayPos);
				 if (d1>d2)
				 {
					int swap = nodesIndexes[i];					
					 nodesIndexes[i] = nodesIndexes[j];
					 nodesIndexes[j] = swap;
				 }
			}
	}

	
    
 void octreeSegmentIntersection(  float4* out_data, float* nodes,  float4* vertexes , 
 int* classFaces , uint* imageOut,
 GridCell* outnodes ,
float4 dim,float4 minV,
int number_Of_Vertexes,int number_of_cells, float4 lightPos,float4 _outdef)
{
    int i,j,k ;
    float4 outP,outPos ;
    float minD ,t1,t2,d;
  
    BoundBox rb; 
    int n;
    int stack[100];
	int list[100];

    int indexCell , nout;
	int indexList;
	GridCell node;

	float4 from, _to ,outNormal;

	 // get index into global data array
   uint x = get_global_id(0);
   uint y = get_global_id(1);
   
   uint maxX = get_global_size(0);
   uint maxY = get_global_size(1);
   
   //----    
	//--

   from = Float4( (float)(x*dim.x/maxX)+0.5,
                           (float)(y*dim.y/maxY)+0.5,
                           dim.z+1200.0) ;

   _to = Float4( (float)(x*dim.x/maxX)+0.5,
                       (float)(y*dim.y/maxY)+0.5,
                           minV.z-1200.0) ;
    

  indexCell = 0;
  indexList = 0;
  rb =  calcBound(from,_to-from);
  push(0,stack,&indexCell);

  minD = 50000000;
  int result = -1 ;
  bool found;
  int intersection = 0;
  float4 int0,int1,v0,v1,v2;


  for (int iter = 0;iter<2;iter++)
  {
      found= false;
	  
		 while (indexCell>0)
		 {
			n = pop(stack, &indexCell);
			if (n<0)  break;
			if (n>=number_of_cells) break;
			
			node = packToCell(nodes,32*n); 	

			// is Leaf
			if (node.childIndex<0)
			 {
			   push(n,list,&indexList);
			}
			else
			{

			  for (int i = 0; i<8;i++)
			  {
				 int nch = (int)(node.childIndex+i);
				 if (nch<0) continue;
				 if (nch>=number_of_cells) continue;
				  GridCell c = packToCell(nodes,32*nch); 	
				 
				  if ( boxOverlap( node.center,node.size,rb.center,rb.size) &&
 					   RayBoxIntersection(c.minBB,c.maxBB,from,_to,&t1,&t2) )
				  {
					push(node.childIndex+i,stack, &indexCell);           
				  }

			 }
			}
		 }

		 
		 sort(list,indexList, nodes,from);

		 for ( i = 0;i<indexList;i++)
		 {
		   node =packToCell(nodes,32*list[i]);// nodes[list[i]];
		   if (node.childIndex<0)
		   {
		               
				for (j = 0;j<node.content - 1;j++)
				{
				   int it =classFaces[ (int)( node.triangleIndex+j )];
				   if ((3*it+2)>=number_Of_Vertexes) continue;
								 
					 v0 = vertexes[3*it];
					 v1 = vertexes[3*it+1];
					 v2 = vertexes[3*it+2];

					outP =intersect(from,_to,v0,v1,v2);
					if ( outP.w<1.0 ) continue;
				  //-- tomo la esfera de los vecinos
				  //if not rbOverlap(rb,t.bound) then continue;
		          
				  d = distance(from, outP)  ;
				  if (d<0.5) continue;

				  if (d<minD)
				  {
					result = it;
					outNormal = sNormal(v0,v1,v2);
					outPos = outP;
					//outPos.w = node.index;
					minD =  d;
					found = true;						
					continue;
				  }
				 }
				 if (result>=0)  break;

		   }
		 }
         if (found)
		 {
		       if (intersection==0) 
						 {
						   int0 = outPos;
						   from = outPos;
						   _to = lightPos;
						   outNormal = sNormal(v0,v1,v2);
						   push(0,stack,&indexCell);
						   intersection = 1;
						   
						 }
						 else
						 if (intersection==1) 
						 {
						    int1 = outP;
							intersection = 2;
						    
						 }	
		 }
  }
 
 //barrier(CLK_GLOBAL_MEM_FENCE);
	// Save Data
	if (intersection==0)
	{
					 imageOut[(int)(y*_outdef.x+x)] = 0;
	}	
	else
	{
						
			float dd = dot( sNormalize(int0-lightPos) , outNormal*-1.0) *0.5+0.5; 
            float4 fColor; 
			if (intersection>1)
				fColor = Float4(dd/4.0,dd/4.0,dd/4.0);
			else
				fColor = Float4(dd,dd,dd);			
				
			out_data[(int)(y*_outdef.x*3+x*3)] = int0;
			out_data[(int)(y*_outdef.x*3+x*3)].w = intersection;
            out_data[(int)(y*_outdef.x*3+x*3)+1] = outNormal;
            imageOut[(int)(y*_outdef.x+x)] = rgbaFloatToInt(fColor) ;
     
	}
		
 


}