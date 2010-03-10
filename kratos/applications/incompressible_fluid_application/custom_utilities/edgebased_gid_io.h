/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER. 

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2009-01-13 10:35:33 $
//   Revision:            $Revision: 1.4 $
//
//

#if !defined(KRATOS_EDGEBASED_GID_IO_BASE_H_INCLUDED)
#define  KRATOS_EDGEBASED_GID_IO_BASE_H_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <cstddef>


// Project includes
#include "includes/define.h"
//#include "includes/datafile_io.h"
#include "includes/gid_io.h"

namespace Kratos
{
    /**
     * Type definitions
     */
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef GeometryData::IntegrationMethod IntegrationMethodType;
    typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;
    
    /**
     * Auxiliary class to store gauss point containers and perform result printing
     * on gauss points
     */
    class EdgebasedGidGaussPointsContainer
    {
        public:
            ///Constructor
            EdgebasedGidGaussPointsContainer( const char* gp_title, KratosGeometryFamily geometryFamily,
                                     GiD_ElementType gid_element_type,
                                     int number_of_integration_points,
                                     std::vector<int> index_container )
            :mGPTitle(gp_title),mKratosElementFamily(geometryFamily),
                      mGidElementFamily(gid_element_type), mSize(number_of_integration_points),
                                        mIndexContainer(index_container){}
            
            bool AddElement( const ModelPart::ElementsContainerType::iterator pElemIt )
            {
                KRATOS_TRY
                if( pElemIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily
                    && pElemIt->GetGeometry().IntegrationPoints( 
                                            pElemIt->GetIntegrationMethod() ).size() == mSize )
                {
                    mMeshElements.push_back( *(pElemIt.base() ) );
                    return true;
                }
                else return false;
                KRATOS_CATCH("")
            }
            
            bool AddCondition( const ModelPart::ConditionsContainerType::iterator pCondIt )
            {
                KRATOS_TRY
                if( pCondIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily
                    && pCondIt->GetGeometry().IntegrationPoints().size() == mSize )
                {
                    mMeshConditions.push_back( *(pCondIt.base() ) );
                    return true;
                }
                else return false;
                KRATOS_CATCH("")
            }
            
            void PrintResults( Variable<double> rVariable, ModelPart& r_model_part, 
                               double SolutionTag, int value_index = 0 )
            {
                if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
                {
                    WriteGaussPoints();
                    GiD_BeginResult( (char *)(rVariable.Name()).c_str(), (const char *)("Kratos"), SolutionTag,
                                      GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
                    std::vector<double> ValuesOnIntPoint(mSize);
                    if( mMeshElements.size() != 0 )
                    {
                        for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); 
                             it != mMeshElements.end(); it++ )
                        {
                            it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                    r_model_part.GetProcessInfo() );
                            for(unsigned int i=0; i<mIndexContainer.size(); i++)
                            {
                                int index = mIndexContainer[i];
                                GiD_WriteScalar( it->Id(), ValuesOnIntPoint[index] );
                            }
                        }
                    }
                    if( mMeshConditions.size() != 0 )
                    {
                        for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                             it != mMeshConditions.end(); it++ )
                        {
                            it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                    r_model_part.GetProcessInfo() );
                            for(unsigned int i=0; i<mIndexContainer.size(); i++)
                            {
                                int index = mIndexContainer[i];
                                GiD_WriteScalar( it->Id(), ValuesOnIntPoint[index] );
                            }
                        }
                    }
                    GiD_EndResult();
                }
            }
            
            void PrintResults( Variable<Vector> rVariable, ModelPart& r_model_part, 
                               double SolutionTag, int value_index = 0 )
            {
                if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
                {
                    WriteGaussPoints();
                    GiD_BeginResult( (char *)(rVariable.Name()).c_str(), (const char*)("Kratos"), SolutionTag,
                                      GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
                    std::vector<Vector> ValuesOnIntPoint(mSize);
                    if( mMeshElements.size() != 0 )
                    {
                        for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); 
                         it != mMeshElements.end(); ++it )
                        {
                            it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                    r_model_part.GetProcessInfo() );
                            for(unsigned int i=0; i<mIndexContainer.size(); i++)
                            {
                                int index = mIndexContainer[i];
                                if( ValuesOnIntPoint[0].size() == 3 )
                                    GiD_WriteVector( it->Id(), ValuesOnIntPoint[index][0],
                                        ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                            }
                        }
                    }
                    if( mMeshConditions.size() != 0 )
                    {
                        for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                             it != mMeshConditions.end(); it++ )
                        {
                            it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                    r_model_part.GetProcessInfo() );
                            for(unsigned int i=0; i<mIndexContainer.size(); i++)
                            {
                                int index = mIndexContainer[i];
                                if( ValuesOnIntPoint[0].size() == 3 )
                                    GiD_WriteVector( it->Id(), ValuesOnIntPoint[index][0],
                                        ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                            }
                        }
                    }
                    GiD_EndResult();
                }
            }
            
            void PrintResults( Variable<Matrix> rVariable, ModelPart& r_model_part, 
                               double SolutionTag, int value_index = 0 )
            {
                if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
                {
                    WriteGaussPoints();
                    GiD_BeginResult( (char *)(rVariable.Name()).c_str(), (const char*)("Kratos"), SolutionTag,
                                      GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
                    std::vector<Matrix> ValuesOnIntPoint(mSize);
                    if( mMeshElements.size() != 0 )
                    {
                        for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); 
                             it != mMeshElements.end(); ++it )
                        {
                            it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                    r_model_part.GetProcessInfo() );
                            for(unsigned int i=0; i<mIndexContainer.size(); i++)
                            {
                                int index = mIndexContainer[i];
                                if(ValuesOnIntPoint[index].size1() ==3 
                                   && ValuesOnIntPoint[index].size2() ==3)
                                    GiD_Write3DMatrix( it->Id(), ValuesOnIntPoint[index](0,0),
                                     ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                     ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                     ValuesOnIntPoint[index](0,2) );
                                if(ValuesOnIntPoint[index].size1() ==1 
                                   && ValuesOnIntPoint[index].size2() ==6)
                                    GiD_Write3DMatrix( it->Id(), ValuesOnIntPoint[index](0,0),
                                     ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                     ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                     ValuesOnIntPoint[index](0,5) );
                            }
                        }
                    }
                    if( mMeshConditions.size() != 0 )
                    {
                        for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                             it != mMeshConditions.end(); it++ )
                        {
                            it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                    r_model_part.GetProcessInfo() );
                            for(unsigned int i=0; i<mIndexContainer.size(); i++)
                            {
                                int index = mIndexContainer[i];
                                if(ValuesOnIntPoint[index].size1() ==3 
                                   && ValuesOnIntPoint[index].size2() ==3)
                                    GiD_Write3DMatrix( it->Id(), ValuesOnIntPoint[index](0,0),
                                        ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                        ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                        ValuesOnIntPoint[index](0,2) );
                                if(ValuesOnIntPoint[index].size1() ==1 
                                   && ValuesOnIntPoint[index].size2() ==6)
                                    GiD_Write3DMatrix( it->Id(), ValuesOnIntPoint[index](0,0),
                                        ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                        ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                        ValuesOnIntPoint[index](0,5) );
                            }
                        }
                    }
                    GiD_EndResult();
                }
            }
            
            void Reset()
            {
                mMeshElements.clear();
                mMeshConditions.clear();
            }
        
        protected:
            void WriteGaussPoints()
            {
                //setting up gauss points
                if( mGidElementFamily == GiD_Tetrahedra && mSize == 5 )
                {
                    GiD_BeginGaussPoint( mGPTitle, GiD_Tetrahedra, NULL, 4, 0, 0 );
                    GiD_WriteGaussPoint3D( 1.0/6.0, 1.0/6.0, 1.0/6.0 );
                    GiD_WriteGaussPoint3D( 1.0/2.0, 1.0/6.0, 1.0/6.0 );
                    GiD_WriteGaussPoint3D( 1.0/6.0, 1.0/2.0, 1.0/6.0 );
                    GiD_WriteGaussPoint3D( 1.0/6.0, 1.0/6.0, 1.0/2.0 );
                    GiD_EndGaussPoint();
                }
                else if( mGidElementFamily == GiD_Tetrahedra && mSize == 10 )
                {
                    GiD_BeginGaussPoint("tet10_element_gp", GiD_Tetrahedra, NULL, 4, 0, 0);
                    GiD_WriteGaussPoint3D( 1.0/14.0, 1.0/14.0, 1.0/14.0 );
                    GiD_WriteGaussPoint3D( 11.0/14.0, 1.0/14.0, 1.0/14.0 );
                    GiD_WriteGaussPoint3D( 1.0/14.0, 11.0/14.0, 1.0/14.0 );
                    GiD_WriteGaussPoint3D( 1.0/14.0, 1.0/14.0, 11.0/14.0 );
                    GiD_EndGaussPoint();
                }
                else
                {
                    GiD_BeginGaussPoint(mGPTitle, mGidElementFamily, NULL, 
                                        mSize, 0, 1);
                    GiD_EndGaussPoint();
                }
            }
            
            ///member variables
            const char* mGPTitle;
            KratosGeometryFamily mKratosElementFamily;
            GiD_ElementType mGidElementFamily;
            unsigned int mSize;
            std::vector<int> mIndexContainer;
            ModelPart::ElementsContainerType mMeshElements;
            ModelPart::ConditionsContainerType mMeshConditions;
    };//class EdgebasedGidGaussPointsContainer
    
    /**
     * Auxiliary class to store meshes of different element types and to
     * write these meshes to an output file
     */
    class EdgebasedGidMeshContainer
    {
        public:
            ///Constructor
            EdgebasedGidMeshContainer( GeometryData::KratosGeometryType geometryType, 
                              GiD_ElementType elementType, const char* mesh_title )
            :mGeometryType(geometryType), mGidElementType(elementType), mMeshTitle(mesh_title){}
            
            bool AddElement( const ModelPart::ElementsContainerType::iterator pElemIt )
            {
                KRATOS_TRY
                if( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
                {
                    mMeshElements.push_back( *(pElemIt.base() ) );
                    Geometry<Node<3> >&geom = pElemIt->GetGeometry();
                    for( Element::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
                    {
                        mMeshNodes.push_back( *(it.base() ) );
                    }
                    return true;
                }
                else
                    return false;
                KRATOS_CATCH("")
            }
            
            bool AddCondition(const ModelPart::ConditionsContainerType::iterator pCondIt)
            {
                KRATOS_TRY
                        if( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
                {
                    mMeshConditions.push_back( *(pCondIt.base() ) );
                    Geometry<Node<3> >&geom = pCondIt->GetGeometry();
                    for( Condition::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
                    {
                        mMeshNodes.push_back( *(it.base() ) );
                    }
                    return true;
                }
                else
                    return false;
                KRATOS_CATCH("")
            }
            
            void FinalizeMeshCreation()
            {
                if( mMeshElements.size() != 0 )
                {
                    mMeshNodes.Unique();
                }
            }
            
            void WriteMesh(bool deformed)
            {
                KRATOS_TRY
                if( mMeshElements.size() != 0 )
                {
                    if( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        std::cout << "writing a 2D mesh" << std::endl;
                        GiD_BeginMesh( "Volume mesh", GiD_2D, mGidElementType,
                                       mMeshElements.begin()->GetGeometry().size() );
                    }
                    else if( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        std::cout << "writing a 3D mesh" << std::endl;
                        GiD_BeginMesh( "Volume mesh", GiD_3D, mGidElementType,
                                       mMeshElements.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_ERROR(std::logic_error,"check working space dimension of model","");
                    //printing nodes
                    GiD_BeginCoordinates();
                    for( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin(); 
                         it != mMeshNodes.end(); ++it )
                    {
                        if( deformed )
                            GiD_WriteCoordinates( (it)->Id(), (it)->X(),
                                (it)->Y(), (it)->Z());
                        else
                            GiD_WriteCoordinates( (it)->Id(), (it)->X0(),
                                (it)->Y0(), (it)->Z0());
                    }
                    GiD_EndCoordinates();
                    //printing elements
                    GiD_BeginElements();
                    int* nodes_id = new int[mMeshElements.begin()->GetGeometry().size()+1];
                    for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); 
                         it != mMeshElements.end(); ++it )
                    {
                        for( unsigned int i=0; i<(it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();
                        
                        //setting the color for either fluid or porous medium:
			double n_dist=0;
			unsigned int n_por=0;
			unsigned int color=13;
			
			for ( unsigned int i=0; i<(it)->GetGeometry().size(); i++ )
			{
				n_dist += it->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				if(     it->GetGeometry()[i].FastGetSolutionStepValue(POROSITY) == 1.0)
				  {    n_por = (it)->GetGeometry().size();}
// 				n_por += it->GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
			}
//  			n_dist /= (it)->GetGeometry().size();
// 			n_por /= (it)->GetGeometry().size();



			if (n_dist> 0.0 && n_por== (it)->GetGeometry().size())//air
				{
				    color=7;			
				}
			else if (n_dist<= 0.0 && n_por== (it)->GetGeometry().size())//free fluid 
				{
				    color=4;				
				}
			else if (n_dist> 0.0 && n_por< (it)->GetGeometry().size())//porous medium 
				{
				    color=5;				
				}
			else //fluid into the porous medium 
				{
				    color=20;		
				}

		
			
			nodes_id[(it)->GetGeometry().size()]= color;
			
                        GiD_WriteElementMat((it)->Id(), nodes_id);
                    }
                    delete [] nodes_id;
                    GiD_EndElements();
                    GiD_EndMesh(); 
                }
		if( mMeshConditions.size() != 0 ) 
                {
                KRATOS_WATCH( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() )   
			 
		    if( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        std::cout << "writing a 3D mesh of the faces" << std::endl;
                        GiD_BeginMesh( "Surface Structure Mesh", GiD_3D, GiD_Triangle,  3);
                    }
		    else                  
			   KRATOS_ERROR(std::logic_error,"Check your space dimensions","");
                    //printing nodes
                    GiD_BeginCoordinates();
			
                    GiD_EndCoordinates();
                    //printing elements
                    GiD_BeginElements();
			//for every face of tetrahedron we create a list of its nodes	
                    int* nodes_id = new int[4];
                   
		   
                    for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); 
                         it != mMeshConditions.end(); ++it )
                    {
                        for( unsigned int i=0; i<(it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();

			unsigned int n_str=0;		
			
			unsigned int n_free_surf=0;	
			
			for (int i=0;i<3;i++)
			{
			n_free_surf+=(it)->GetGeometry()[i].FastGetSolutionStepValue(IS_FREE_SURFACE);
			n_str+=(it)->GetGeometry()[i].FastGetSolutionStepValue(IS_STRUCTURE);
			}       
			if (n_str==it->GetGeometry().size() && n_free_surf!=it->GetGeometry().size())
				{
				nodes_id[3]=3;
				 GiD_WriteElementMat((it)->Id(), nodes_id);			
				}
		
                    }
                    delete [] nodes_id;
                    GiD_EndElements();
                    GiD_EndMesh();
		}
		if( mMeshConditions.size() != 0 ) 
	        {
	                KRATOS_WATCH( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() )    
		    if( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        std::cout << "writing a 3D mesh of the faces" << std::endl;
                        GiD_BeginMesh( "Surface Fluid Mesh", GiD_3D, GiD_Triangle,  3);
                    }
		    else                  
			   KRATOS_ERROR(std::logic_error,"Check your space dimensions","");
			
			//now writing the fluid surface mesh
			//printing nodes
                    GiD_BeginCoordinates();
			
                    GiD_EndCoordinates();
                    //printing elements
                    GiD_BeginElements();
			//for every face of tetrahedron we create a list of its nodes	
                    //int* nodes_id = new int[4];
                   
		   
			 int* nodes_id = new int[4];
			unsigned int n_fl = 0;
			unsigned int n_str = 0;
	

                    for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); 
                         it != mMeshConditions.end(); ++it )
                    {
                        for( unsigned int i=0; i<(it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();

			
			for (int i=0;i<3;i++)
			{
			n_str+=(it)->GetGeometry()[i].FastGetSolutionStepValue(IS_STRUCTURE);
			n_fl+=(it)->GetGeometry()[i].FastGetSolutionStepValue(IS_FLUID);
			}       
			//if (n_free_surf==it->GetGeometry().size())
			//	nodes_id[3]=1;
			if (n_str!=it->GetGeometry().size() && n_fl>0)
				{	
				//the color of water
				nodes_id[3]=4;
				 GiD_WriteElementMat((it)->Id(), nodes_id);			
				}
		
                    }
                    delete [] nodes_id;
                    GiD_EndElements();
                    GiD_EndMesh();
                }
                
                KRATOS_CATCH("")
            }

            void Reset()
            {
                mMeshNodes.clear();
                mMeshElements.clear();
                mMeshConditions.clear();
            }
            
        private:
            ///member variables
            GeometryData::KratosGeometryType mGeometryType;
            GiD_ElementType mGidElementType;
            ModelPart::NodesContainerType mMeshNodes;
            ModelPart::ElementsContainerType mMeshElements;
            ModelPart::ConditionsContainerType mMeshConditions;
            const char* mMeshTitle;
    };//class EdgebasedGidMeshContainer

}// namespace Kratos.


#endif // KRATOS_EDGEBASED_GID_IO_BASE_H_INCLUDED  defined 
