Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData



Begin Properties  0
End Properties

Begin Nodes
*tcl(DEM::WriteNonSphereAndCircleNodes *FileId)*\
End Nodes
*tcl(DEM::ReleaseSphereAndCircleNodes)

*set elems(All)
*if( GenData(Domain_Dimension,int) == 3 )
Begin Conditions RigidFace3D3N
*loop elems *all
*if(strcmp(ElemsTypeName,"Triangle") == 0)
*Set var i=0
*set var j= ElemsNnode
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*operation(ElemsConec(*i))*\
*end for

*endif
*end elems
End Conditions
*endif

*set elems(All)
*if( GenData(Domain_Dimension,int) == 3 )
Begin Conditions RigidFace3D4N
*loop elems *all
*if(strcmp(ElemsTypeName,"Quadrilateral") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Conditions
*endif

*set elems(All)
*if( GenData(Domain_Dimension,int) == 2 )
Begin Conditions RigidEdge3D2N
*loop elems *all
*if(strcmp(ElemsTypeName,"Linear") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Conditions
*endif


Begin Table 1 TIME TEMPERATURE
0.0  0.0
1.0  0.0
End Table

Begin Table 2 TIME TEMPERATURE
0.0  0.0
1.0  0.0
End Table

Begin Table 3 TIME TEMPERATURE
0.0  0.0
1.0  0.0
End Table


Begin Mesh 1
Begin MeshNodes

End MeshNodes
End Mesh

Begin ConditionalData GROUP_ID
End ConditionalData