#include "stdafx.h"
#include "u_colorScale.h"

float4 TColorToRealColor(int ct)
{
 float4 result;
 result.z = (float)((abs(ct) / (256*256))/255);
 result.y = (float)(((abs(ct) % (256*256)) / 256)/255)  ;
 result.x = (float)((abs(ct) % 256)/255);
 result.w =  1.0f;
 return result;
}


int TVector4fToTColor(float4 aV )
{
 return (int)((aV.z*255) * 256*256 + (aV.y*255) *256 + (aV.x*255));
}
