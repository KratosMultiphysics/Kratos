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
//   Date:                $Date: 2009-01-13 16:51:36 $
//   Revision:            $Revision: 1.39 $
//
//

#if !defined(KRATOS_GID_IO_BASE_H_INCLUDED)
#define  KRATOS_GID_IO_BASE_H_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <cstddef>
#include <iomanip>

// External includes 
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/define.h"
#include "includes/datafile_io.h"
#include "geometries/geometry_data.h"

#include "includes/gid_gauss_point_container.h"
#include "includes/gid_mesh_container.h"

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
    
    ///Flags for mesh writing
    enum WriteDeformedMeshFlag{WriteDeformed, WriteUndeformed};
    enum WriteConditionsFlag{WriteConditions, WriteElementsOnly};
    enum MultiFileFlag{SingleFile, MultipleFiles};
            
    
    /**
     * This class defines an interface to the GiDPost library
     * in order to provide GiD compliant I/O functionality
     */
    template<class TGaussPointContainer = GidGaussPointsContainer, class TMeshContainer = GidMeshContainer>
    class GidIO : public DatafileIO
    {
        public:
            ///pointer definition of GidIO
            KRATOS_CLASS_POINTER_DEFINITION(GidIO);
            
            ///typedefs
            typedef DatafileIO BaseType;
            
            ///Flags for mesh writing
//             enum WriteDeformedMeshFlag{WriteDeformed, WriteUndeformed};
//             enum WriteConditionsFlag{WriteConditions, WriteElementsOnly};
//             enum MultiFileFlag{SingleFile, MultipleFiles};
            
            ///Constructor
            GidIO( const std::string& rNodeDatafilename, 
                   const std::string& rPropertiesDatafilename,
                   const std::string& rElementDatafilename, 
                   const std::string& rConditionDatafilename, 
                   const std::string& rInitialValueDatafilename,
                   const  std::string& rResultFilename, 
                   GiD_PostMode Mode,
                   MultiFileFlag use_multiple_files_flag,
                   WriteDeformedMeshFlag write_deformed_flag,
                   WriteConditionsFlag write_conditions_flag
                 )
            : BaseType( rNodeDatafilename, rPropertiesDatafilename, 
                        rElementDatafilename, rConditionDatafilename, 
                        rInitialValueDatafilename
                      )
            {
                mMode = Mode;
                mResultFileOpen = false;
                mMeshFileOpen = false;
                mWriteDeformed = write_deformed_flag;
                mWriteConditions = write_conditions_flag;
                mUseMultiFile = use_multiple_files_flag;
                mResultFileName = rResultFilename;
                InitializeResultFile(mResultFileName);
                mMeshFileName = mResultFileName;
                mMeshFileName += ".post.msh";
                SetUpMeshContainers();
                SetUpGaussPointContainers();
            }
            
            ///single stream IO constructor
            GidIO( const std::string& rDatafilename, 
                   GiD_PostMode Mode,
                   MultiFileFlag use_multiple_files_flag,
                   WriteDeformedMeshFlag write_deformed_flag,
                   WriteConditionsFlag write_conditions_flag
                 )
            : BaseType(rDatafilename)
            {
                mMode = Mode;
                mResultFileOpen = false;
                mMeshFileOpen = false;
                mWriteDeformed = write_deformed_flag;
                mWriteConditions = write_conditions_flag;
                mUseMultiFile = use_multiple_files_flag;
                mResultFileName = rDatafilename;
                InitializeResultFile(mResultFileName);
                mMeshFileName = rDatafilename;
//                 mMeshFileName += ".post.msh";
                SetUpMeshContainers();
                SetUpGaussPointContainers();
            }
            
            ///Destructor.
            virtual ~GidIO()
            {
                if( mResultFileOpen )
                {
                    GiD_ClosePostResultFile();
                    mResultFileOpen = false;
                }
            }
            
            ///initialization functions
            /**
             * creates the mesh containers for all different element types.
             * Note that the containers are not filled yet in here!
             */
            void SetUpMeshContainers()
            {
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Hexahedra3D20, 
                                              GiD_Hexahedra, "Kratos_Hexahedra3D20_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Hexahedra3D27, 
                                              GiD_Hexahedra, "Kratos_Hexahedra3D27_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Hexahedra3D8, 
                                              GiD_Hexahedra, "Kratos_Hexahedra3D8_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Prism3D15, 
                                              GiD_Prism, "Kratos_Prism3D15_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Prism3D6, 
                                              GiD_Prism, "Kratos_Prism3D6_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Quadrilateral2D4, 
                                              GiD_Quadrilateral, "Kratos_Quadrilateral2D4_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Quadrilateral2D8, 
                                              GiD_Quadrilateral, "Kratos_Quadrilateral2D8_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Quadrilateral2D9, 
                                              GiD_Quadrilateral, "Kratos_Quadrilateral2D9_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Quadrilateral3D4, 
                                              GiD_Quadrilateral, "Kratos_Quadrilateral3D4_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Quadrilateral3D8, 
                                              GiD_Quadrilateral, "Kratos_Quadrilateral3D8_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Quadrilateral3D9, 
                                              GiD_Quadrilateral, "Kratos_Quadrilateral3D9_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Tetrahedra3D10, 
                                              GiD_Tetrahedra, "Kratos_Tetrahedra3D10_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Tetrahedra3D4, 
                                              GiD_Tetrahedra, "Kratos_Tetrahedra3D4_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Triangle2D3, 
                                              GiD_Triangle, "Kratos_Triangle2D3_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Triangle2D6, 
                                              GiD_Triangle, "Kratos_Triangle2D6_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Triangle3D3, 
                                              GiD_Triangle, "Kratos_Triangle3D3_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Line2D2, 
                                              GiD_Linear, "Kratos_Line2D2_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Line3D2, 
                                              GiD_Linear, "Kratos_Line3D2_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Line2D3, 
                                              GiD_Linear, "Kratos_Line2D3_Mesh" ) );
                mGidMeshContainers.push_back( TMeshContainer(
                                              GeometryData::Kratos_Line3D3, 
                                              GiD_Linear, "Kratos_Line3D3_Mesh" ) );
            }//SetUpMeshContainers
            
            /**
             * creates the gauss point containers for all different element types.
             * Note that the containers are not filled yet in here!
             */
            virtual void SetUpGaussPointContainers()
            {
                //elements with 1 gauss point
                std::vector<int> gp_indices(1);
                gp_indices[0] = 0;
                //case Triangle with 1 gauss point
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tri1_element_gp",
                        GeometryData::Kratos_Triangle, GiD_Triangle, 1, gp_indices ) );
                //case Quadrilateral with 1 gauss point
                mGidGaussPointContainers.push_back( TGaussPointContainer( "quad1_element_gp",
                        GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, 1, gp_indices ) );
                //case Tetrahedra with 1 gauss point
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tet1_element_gp",
                        GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, 1, gp_indices ) );
                //case Hexahedra with 1 gauss point
                mGidGaussPointContainers.push_back( TGaussPointContainer( "hex1_element_gp",
                        GeometryData::Kratos_Hexahedra, GiD_Hexahedra, 1, gp_indices ) );
                //case Prism with 1 gauss point
                mGidGaussPointContainers.push_back( TGaussPointContainer( "prism1_element_gp",
                        GeometryData::Kratos_Prism, GiD_Prism, 1, gp_indices ) );
                
                //elements with 3 gauss points
                gp_indices.resize(3);
                //case Triangle with 3 gauss points
                gp_indices[0] = 0;      gp_indices[1] = 1;      gp_indices[2] = 2;
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tri3_element_gp",
                        GeometryData::Kratos_Triangle, GiD_Triangle, 3, gp_indices ) );
                
                //elements with 4 gauss points
                gp_indices.resize(4);
                gp_indices[0] = 0;      gp_indices[1] = 1;
                gp_indices[2] = 2;      gp_indices[3] = 3;
                //case Quadrilateral with 4 gauss points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "quad4_element_gp",
                        GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, 4, gp_indices ) );
                //case Tetrahedra with 4 gauss points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tet4_element_gp",
                        GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, 4, gp_indices ) );
                //case Triangle with 4 gauss points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tri4_element_gp",
                        GeometryData::Kratos_Triangle, GiD_Triangle, 4, gp_indices ) );
                //case Tetrahedra with 5 gauss points (4 gauss points will be created for GiD)
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tet5_element_gp",
                        GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, 5, gp_indices ) );
                //case Tetrahedra with 11 gauss points (4 gauss points will be created for GiD)
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tet11_element_gp",
                        GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, 11, gp_indices ) );
                //case Tetrahedra with 10 gauss points (4 gauss points will be created for GiD)
                gp_indices.resize(10);
                gp_indices[0] = 0;      gp_indices[1] = 1;      gp_indices[2] = 2;
                gp_indices[3] = 3;      gp_indices[4] = 4;      gp_indices[5] = 5;
                gp_indices[6] = 6;      gp_indices[7] = 7;      gp_indices[8] = 8;
                gp_indices[9] = 9;
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tet10_element_gp",
                        GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, 10, gp_indices ) );
                
                //elements with 6 gauss points
                gp_indices.resize(6);
                gp_indices[0] = 0;      gp_indices[1] = 1;      gp_indices[2] = 2;
                gp_indices[3] = 3;      gp_indices[4] = 4;      gp_indices[5] = 5;
                //case Triangle with 6 gauss points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "tri6_element_gp",
                        GeometryData::Kratos_Triangle, GiD_Triangle, 6, gp_indices ) );
                //case Prism with 6 Gauss Points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "prism6_element_gp",
                        GeometryData::Kratos_Prism, GiD_Prism, 6, gp_indices ) );
                
                //elements with 8 gauss points
                gp_indices.resize(8);
                gp_indices[0] = 0;      gp_indices[1] = 1;      gp_indices[2] = 2;
                gp_indices[3] = 3;      gp_indices[4] = 4;      gp_indices[5] = 5;
                gp_indices[6] = 6;      gp_indices[7] = 7;
                //case Hexahedra with 8 gauss points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "hex8_element_gp",
                        GeometryData::Kratos_Hexahedra, GiD_Hexahedra, 8, gp_indices ) );
                
                //elements with 9 gauss points
                gp_indices.resize(9);
                gp_indices[0] = 0;      gp_indices[1] = 1;      gp_indices[2] = 2;
                gp_indices[3] = 3;      gp_indices[4] = 4;      gp_indices[5] = 5;
                gp_indices[6] = 6;      gp_indices[7] = 7;      gp_indices[8] = 8;
                //case Prism with 9 Gauss Points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "prism9_element_gp",
                        GeometryData::Kratos_Prism, GiD_Prism, 9, gp_indices ) );
                // case quadrilateral with 9 Gauss Points 
                mGidGaussPointContainers.push_back( TGaussPointContainer( "quad9_element_gp",
                        GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, 9, gp_indices ) );
                
                //elements with 27 gauss points
                gp_indices.resize(27);
                gp_indices[0] = 0;      gp_indices[8] = 1;      gp_indices[1] = 2;
                gp_indices[11] = 3;     gp_indices[20] = 4;     gp_indices[9] = 5;
                gp_indices[3] = 6;      gp_indices[10] = 7;     gp_indices[2] = 8;
                gp_indices[12] = 9;     gp_indices[21] = 10;    gp_indices[13] = 11;
                gp_indices[24] = 12;    gp_indices[26] = 13;    gp_indices[22] = 14;
                gp_indices[15] = 15;    gp_indices[23] = 16;    gp_indices[14] = 17;
                gp_indices[4] = 18;     gp_indices[16] = 19;    gp_indices[5] = 20;
                gp_indices[19] = 21;    gp_indices[25] = 22;    gp_indices[17] = 23;
                gp_indices[7] = 24;     gp_indices[18] = 25;    gp_indices[6] = 26;
                //case Hexahedra with 27 Gauss Points
                mGidGaussPointContainers.push_back( TGaussPointContainer( "hex27_element_gp",
                        GeometryData::Kratos_Hexahedra, GiD_Hexahedra, 27, gp_indices ) );
            }//SetUpGaussPointContainers
            
            
            ///general GidIO related functions
            /**
             * TODO: to be removed
             */
            void ChangeOutputName(const std::string& rDatafilename )
            {
                KRATOS_TRY
                mMeshFileName = rDatafilename;
                mResultFileName = rDatafilename;
                KRATOS_CATCH("")
            }
            
            /**
             * sets up the file names and opens the result file in case there
             * is ASCII mode and only one file written
             */
            void InitializeResultFile( std::string const& rResultFileName )
            {
                std::cout << "initializing result files" << std::endl;
                mResultFileName = rResultFileName;
            }
            
            /**
             * TODO: check whether this is still necessary!
             */
            void  CloseResultFile()
            {
                if( mResultFileOpen )
                    GiD_ClosePostResultFile();
            }
            
            /**
             * TODO: check whether this is still necessary!
             */
            void Flush()
            {
                GiD_FlushPostFile();
            }
            
            /**
             * Turn back information as a string.
             */
            virtual std::string Info() const
            {
                return "gid io";
            }
            
            /**
             * Print information about this object.
             */
            virtual void PrintInfo(std::ostream& rOStream) const
            {
                rOStream << Info();
            }
            
            /**
             * Print object's data.
             */
            virtual void PrintData(std::ostream& rOStream) const
            {
                BaseType::PrintData(rOStream);
            }

            ///result functions
            /**
             * This has to be called for each solution step BEFORE any results
             * (on nodes and on gauss points) is written
             * @param SolutionTag the current solution step (i.e. time)
             * @param conditions_flag states whether results should also be written on conditions
             */
            virtual void InitializeResults( double name, MeshType rThisMesh )
            {
                if( mMode == GiD_PostAscii && ! mResultFileOpen )
                {
                    std::stringstream file_name;
                    file_name << mResultFileName << "_" << name << ".post.res";
                    GiD_OpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
                    mResultFileOpen = true;
                }
                if( mUseMultiFile == SingleFile && mMode == GiD_PostBinary )
                {
                    std::stringstream step_index;
                    step_index << "step" << name;
                    GiD_BeginOnMeshGroup( (char *)(step_index.str()).c_str() );
                }
                //initializing gauss points containers
                for( MeshType::ElementIterator element_iterator = rThisMesh.ElementsBegin();
                     element_iterator != rThisMesh.ElementsEnd(); ++element_iterator )
                    for( typename std::vector<TGaussPointContainer>::iterator it =
                         mGidGaussPointContainers.begin();
                         it != mGidGaussPointContainers.end(); it++ )
                        if( it->AddElement( element_iterator ) )
                            break;
                
                if( mWriteConditions == WriteConditions )
                    for( MeshType::ConditionsContainerType::iterator conditions_iterator =
                         rThisMesh.ConditionsBegin(); conditions_iterator 
                         != rThisMesh.ConditionsEnd(); conditions_iterator++ )
                        for( typename std::vector<TGaussPointContainer>::iterator it =
                             mGidGaussPointContainers.begin();
                             it != mGidGaussPointContainers.end(); it++ )
                            if( it->AddCondition( conditions_iterator ) )
                                break;
            }
            
            /**
             * This has to be called for each solution step AFTER all the results
             * have been written
             */
            void FinalizeResults()
            {
                if( mUseMultiFile == MultipleFiles || mMode == GiD_PostAscii )
                {
                    GiD_ClosePostResultFile();
                    mResultFileOpen = false;
                }
                if( mUseMultiFile == SingleFile && mMode == GiD_PostBinary )
                {
                    GiD_EndOnMeshGroup();
                }
                //resetting gauss point containers
                for( typename std::vector<TGaussPointContainer>::iterator it =
                     mGidGaussPointContainers.begin();
                     it != mGidGaussPointContainers.end(); it++ )
                    it->Reset();
            }
            
            ///functions for writing nodal results
            /**
             * writes nodal results for variables of type double
             */
            void WriteNodalResults( Variable<double> const& rVariable, 
                                    NodesContainerType& rNodes, double SolutionTag, 
                                    std::size_t SolutionStepNumber)
            {
                GiD_BeginResult( (char*)(rVariable.Name().c_str()), "Kratos", 
                                  SolutionTag, GiD_Scalar, 
                                  GiD_OnNodes, NULL, NULL, 0, NULL );
                for( NodesContainerType::iterator i_node = rNodes.begin();
                     i_node != rNodes.end() ; ++i_node)
                    GiD_WriteScalar( i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                                     SolutionStepNumber) );
                GiD_EndResult();
            }
            
            /**
             * writes nodal results for variables of type array_1d<double, 3>
             * (e.g. DISPLACEMENT)
             */
            void WriteNodalResults( Variable<array_1d<double, 3> > const& rVariable,
                                    NodesContainerType& rNodes, 
                                    double SolutionTag, std::size_t SolutionStepNumber)
            {
                GiD_BeginResult( (char*)(rVariable.Name().c_str()), "Kratos", 
                                  SolutionTag, GiD_Vector, 
                                  GiD_OnNodes, NULL, NULL, 0, NULL );
                for (NodesContainerType::iterator i_node = rNodes.begin(); 
                     i_node != rNodes.end() ; ++i_node)
                {
                    array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                            SolutionStepNumber );
                    GiD_WriteVector( i_node->Id(), temp[0], temp[1], temp[2] );
                }
                GiD_EndResult();
            }
            
            
            /**
             * writes nodal results for variables of type Vector
             * (note that only vectors with 3 components can be printed)
             */
            void WriteNodalResults( Variable<Vector> const& rVariable, 
                                    NodesContainerType& rNodes, 
                                    double SolutionTag, std::size_t SolutionStepNumber)
            {
                GiD_BeginResult( (char*)(rVariable.Name().c_str()), "Kratos", 
                                  SolutionTag, GiD_Matrix, 
                                  GiD_OnNodes, NULL, NULL, 0, NULL );
                for (NodesContainerType::iterator i_node = rNodes.begin(); 
                     i_node != rNodes.end() ; ++i_node)
                {
                    Vector& tempVector = i_node->GetSolutionStepValue(rVariable,
                            SolutionStepNumber);
                    if(tempVector.size() ==3 )
                        GiD_WriteVector(  i_node->Id(), tempVector(0), tempVector(1), tempVector(2) );
                    else if(tempVector.size() == 6 )
                        GiD_Write3DMatrix(  i_node->Id(), tempVector(0), tempVector(1), tempVector(2),
                                            tempVector(3), tempVector(4), tempVector(5) );
                }
                GiD_EndResult();
            }
            
            /**
             * writes nodal results for variables of type Matrix
             */
            void WriteNodalResults( Variable<Matrix> const& rVariable, 
                                    NodesContainerType& rNodes, 
                                    double SolutionTag, std::size_t SolutionStepNumber)
            {
                GiD_BeginResult( (char*)(rVariable.Name().c_str()), "Kratos", 
                                  SolutionTag, GiD_Matrix, 
                                  GiD_OnNodes, NULL, NULL, 0, NULL );
                for (NodesContainerType::iterator i_node = rNodes.begin(); 
                     i_node != rNodes.end() ; ++i_node)
                {
                    Matrix& tempMatrix = i_node->GetSolutionStepValue(rVariable,
                            SolutionStepNumber);
                    if(tempMatrix.size1() ==3 && tempMatrix.size2() ==3)
                        GiD_Write3DMatrix(  i_node->Id(), tempMatrix(0,0), tempMatrix(1,1),
                                            tempMatrix(2,2), tempMatrix(0,1), tempMatrix(1,2),
                                                    tempMatrix(0,2) );
                }
                GiD_EndResult();
            }
            
            ///mesh writing functions
            /**
             * opens a new mesh group
             */
            void InitializeMesh( double name )
            {
                if( mUseMultiFile == MultipleFiles )
                {
                    if( mMode == GiD_PostAscii && ! mMeshFileOpen )
                    {
                        std::stringstream file_name;
                        file_name << std::setprecision(12) << mMeshFileName << "_" << name << ".post.msh";
                        GiD_OpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                        mMeshFileOpen = true;
                    }
                    if( mMode == GiD_PostBinary && ! mResultFileOpen )
                    {
                        std::stringstream file_name;
                        file_name << mResultFileName << "_" << name << ".post.bin";
                        if( ! mResultFileOpen )
                        {
                            GiD_OpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
                            mResultFileOpen = true;
                        }
                    }
                }
                if( mUseMultiFile == SingleFile )
                {
                    if( mMode == GiD_PostBinary && ! mResultFileOpen )
                    {
                        std::stringstream file_name;
                        file_name << mResultFileName << ".post.bin";
                        if ( GiD_OpenPostResultFile((char*)(file_name.str()).c_str(), mMode) )
                        {
                            std::stringstream buffer;
                            buffer << "error opening results file:" << "/" <<  mResultFileName   << "/";
                            KRATOS_ERROR(std::runtime_error, buffer.str(), "");	
                        }
                        mResultFileOpen = true;
                    }
                    if( mMode == GiD_PostBinary )
                    {
                        std::stringstream step_index;
                        step_index << "step" << name;
                        GiD_BeginMeshGroup( (char *)(step_index.str()).c_str() );
                    }
                    if( mMode == GiD_PostAscii && ! mMeshFileOpen )
                    {
                        std::stringstream file_name;
                        file_name << mMeshFileName << "_" << name << ".post.msh";
                        GiD_OpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                        mMeshFileOpen = true;
                    }
                    
                }
            }
            
            /**
             * closes a mesh group
             */
            void FinalizeMesh()
            {
                if( mUseMultiFile == MultipleFiles && mMode == GiD_PostAscii )
                {
                    GiD_ClosePostMeshFile();
                    mMeshFileOpen = false;
                }
                if( mUseMultiFile == SingleFile && mMode == GiD_PostBinary )
                    GiD_EndMeshGroup();
                if( mUseMultiFile == SingleFile && mMode == GiD_PostAscii )
                {
                    GiD_ClosePostMeshFile();
                    mMeshFileOpen = false;
                }
            }
            
            /**
             * Writes a node mesh.
             * @param rThisMesh the given mesh to be written to the output file
             * @param solution_step the current solution step
             * @param deformed_flag indicates whether the mesh shall be written in deformed
             * or undeformed state
             * @param Mode either GiD_PostAscii (default) or GiD_PostBinary
             */
            void WriteNodeMesh( MeshType& rThisMesh )
            {
                KRATOS_TRY
                GiD_BeginMesh("Kratos Mesh",GiD_3D,GiD_Point,1);
                GiD_BeginCoordinates();
                for( MeshType::NodeIterator node_iterator = rThisMesh.NodesBegin();
                     node_iterator != rThisMesh.NodesEnd();
                     ++node_iterator)
                {
                    if( mWriteDeformed == WriteUndeformed )
                        GiD_WriteCoordinates( node_iterator->Id(), node_iterator->X0(),
                                              node_iterator->Y0(), node_iterator->Z0() );
                    else if( mWriteDeformed == WriteDeformed )
                        GiD_WriteCoordinates( node_iterator->Id(), node_iterator->X(),
                                              node_iterator->Y(), node_iterator->Z() );
                    else 
                        KRATOS_ERROR( std::logic_error,"undefined WriteDeformedMeshFlag","" );
                }
                GiD_EndCoordinates();
                int nodes_id[1];
                GiD_BeginElements();
                for( MeshType::NodeIterator node_iterator = rThisMesh.NodesBegin();
                     node_iterator != rThisMesh.NodesEnd();
                     ++node_iterator)
                {
                    nodes_id[0] = node_iterator->Id();
                    GiD_WriteElement(node_iterator->Id(), nodes_id);
                }
                GiD_EndElements();
                GiD_EndMesh();
                KRATOS_CATCH("")
            }//WriteNodeMesh
            
            
            /**
             * This is a multi-purpose function that writes arbitrary meshes of elements
             * and conditions in either deformed or undeformed state
             * @param rThisMesh the current mesh to be written
             * @param deformed_flag states whether the mesh should be written in deformed configuration
             * @param conditions_flag states whether conditions should also be written
             */
            void WriteMesh( MeshType& rThisMesh )
            {
                KRATOS_TRY
                for( MeshType::ElementIterator element_iterator = rThisMesh.ElementsBegin();
                element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
                        for( typename std::vector<TMeshContainer>::iterator it = mGidMeshContainers.begin();
                             it != mGidMeshContainers.end(); it++ )
                        if( it->AddElement( element_iterator ) )
                            break;
                if( mWriteConditions == WriteConditions )
                    for( MeshType::ConditionsContainerType::iterator conditions_iterator =
                         rThisMesh.ConditionsBegin(); 
                         conditions_iterator != rThisMesh.ConditionsEnd(); conditions_iterator++ )
                        for( typename std::vector<TMeshContainer>::iterator it = mGidMeshContainers.begin();
                             it != mGidMeshContainers.end(); it++ )
                            if( it->AddCondition( conditions_iterator ) )
                                break;
                for( typename std::vector<TMeshContainer>::iterator it = mGidMeshContainers.begin();
                     it != mGidMeshContainers.end(); it++ )
                {
                    it->FinalizeMeshCreation();
                    if( mWriteDeformed == WriteDeformed )
                        it->WriteMesh(true);
                    else if( mWriteDeformed == WriteUndeformed )
                        it->WriteMesh(false);
                    else
                        KRATOS_ERROR( std::logic_error, "undefined WriteDeformedMeshFlag" , "" );
                    it->Reset();
                }
                KRATOS_CATCH("")
            }//WriteMesh
            
            
            ///functions for printing results on gauss points
            /**
             * Prints variables of type double on gauss points of the complete mesh
             * @param rVariable the given variable name
             * @param r_model_part the current model part
             */
            virtual void PrintOnGaussPoints( const Variable<double>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index = 0 )
            {
                KRATOS_TRY;
                for( typename std::vector<TGaussPointContainer>::iterator it =
                     mGidGaussPointContainers.begin();
                     it != mGidGaussPointContainers.end(); it++ )
                {
                    it->PrintResults( rVariable, r_model_part, SolutionTag, value_index );
                }
                KRATOS_CATCH("");
            }
            
            /**
             * Prints variables of type double on gauss points of the complete mesh
             * @param rVariable the given variable name
             * @param r_model_part the current model part
             */
            void PrintOnGaussPoints( const Variable<array_1d<double,3> >& rVariable, ModelPart& r_model_part, double SolutionTag, int value_index = 0 )
            {
                KRATOS_TRY;
                for( typename std::vector<TGaussPointContainer>::iterator it =
                     mGidGaussPointContainers.begin();
                     it != mGidGaussPointContainers.end(); it++ )
                {
                    it->PrintResults( rVariable, r_model_part, SolutionTag, value_index );
                }
                KRATOS_CATCH("");
            }
            
            /**
             * Prints variables of type double on gauss points of the complete mesh
             * @param rVariable the given variable name
             * @param r_model_part the current model part
             */
            virtual void PrintOnGaussPoints( const Variable<Vector>& rVariable, ModelPart& r_model_part,
                                             double SolutionTag, int value_index = 0 )
            {
                KRATOS_TRY;
                for( typename std::vector<TGaussPointContainer>::iterator it =
                     mGidGaussPointContainers.begin();
                     it != mGidGaussPointContainers.end(); it++ )
                {
                    it->PrintResults( rVariable, r_model_part, SolutionTag, value_index );
		    //KRATOS_WATCH(rVariable)
                }
                KRATOS_CATCH("");
            }
            
            /**
             * Prints variables of type double on gauss points of the complete mesh
             * @param rVariable the given variable name
             * @param r_model_part the current model part
             */
            virtual void PrintOnGaussPoints( const Variable<Matrix>& rVariable, ModelPart& r_model_part,
                                             double SolutionTag, int value_index = 0 )
            {
                KRATOS_TRY;
                for( typename std::vector<TGaussPointContainer>::iterator it =
                     mGidGaussPointContainers.begin();
                     it != mGidGaussPointContainers.end(); it++ )
                {
                    it->PrintResults( rVariable, r_model_part, SolutionTag, value_index );
		    //KRATOS_WATCH(r_model_part)
                }
                KRATOS_CATCH("");
            }
        
        protected:
            /**
             * File names
             */
            std::string mResultFileName;
            std::string mMeshFileName;
            
            /**
             * Flags
             */
            WriteDeformedMeshFlag mWriteDeformed;
            WriteConditionsFlag mWriteConditions;
            MultiFileFlag mUseMultiFile;
            GiD_PostMode mMode;
            
            /**
             * member variables
             */
            std::vector<TMeshContainer> mGidMeshContainers;
            std::vector<TGaussPointContainer> mGidGaussPointContainers;
            bool mMeshFileOpen;
            bool mResultFileOpen;
        
        private:    
            /**
             * assignment operator
             */
            GidIO& operator=(GidIO const& rOther);
            
            /**
             * Copy constructor
             */
            GidIO(GidIO const& rOther);
    }; // Class GidIO 


    /**
     * Input and output 
     */
/*    GidIO& operator >> (GidIO& rInput, IO::NodeType& rNode)
    {
        rInput.ReadNode(rNode);
        return rInput;
    }
    
    GidIO& operator >> (GidIO& rInput, IO::NodesContainerType& rNodes)
    {
        rInput.ReadNodes(rNodes);
        return rInput;
    }
    
    GidIO& operator >> (GidIO& rInput, IO::PropertiesContainerType& rProperties)
    {
        rInput.ReadProperties(rProperties);
        return rInput;
    }
    
    GidIO& operator >> (GidIO& rInput, IO::MeshType& rMesh)
    {
        rInput.ReadMesh(rMesh);
        return rInput;
    }
    
    GidIO& operator << (GidIO& rOutput, IO::NodesContainerType& rNodes)
    {
        rOutput.WriteNodes(rNodes);
        return rOutput;
    }
    
    GidIO& operator << (GidIO& rOutput, IO::ElementsContainerType& rElements)
    {
        rOutput.WriteElements(rElements);
        return rOutput;
    }*/
    
    /**
     * output stream function
     */
    inline std::ostream& operator << (std::ostream& rOStream, const GidIO<>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

}// namespace Kratos.

#undef KRATOS_INDEX_PARSER
#undef KRATOS_COORDINATES_PARSER
#undef KRATOS_NODE_PARSER
#undef KRATOS_NODE_DATA_PARSER
#undef KRATOS_PROPERTIES_LHS_PARSER
#undef KRATOS_PROPERTIES_TEMPORARY_VARIABLES
#undef KRATOS_ARRAY_1D_3_PARSER
#undef KRATOS_VECTOR_PARSER
#undef KRATOS_MATRIX_PARSER
#undef KRATOS_CONDITIONS_TEMPORARY_VARIABLES
#undef KRATOS_CONDITIONS_FIX_PARSER

#endif // KRATOS_GID_IO_BASE_H_INCLUDED  defined 
