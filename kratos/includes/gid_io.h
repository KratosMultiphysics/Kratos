//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Pooyan Dadvand
//

#pragma once

// System includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstddef>
#include <iomanip>

// External includes
#define USE_CONST
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/gid_gauss_point_container.h"
#include "includes/gid_mesh_container.h"
#include "utilities/timer.h"
#include "containers/flags.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///Flags for mesh writing
enum WriteDeformedMeshFlag {WriteDeformed, WriteUndeformed};
enum WriteConditionsFlag {WriteConditions, WriteElementsOnly, WriteConditionsOnly};
enum MultiFileFlag {SingleFile, MultipleFiles};

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GidIOBase
 * @ingroup KratosCore
 * @brief Base class for GidIO
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) GidIOBase
    : public IO
{
protected:
    /**
     * Counter of live GidIO instances
     * (to ensure GiD_PostInit and GiD_PostDone are properly called)
     */
    int data;

    // Private constructor so that no objects can be created.
    GidIOBase() {
        data = 0;
    }

public:
    static GidIOBase& GetInstance();

    int GetData();
    void SetData(int data);

private:
    static void Create();

    static GidIOBase* mpInstance;
};

/**
 * @class GidIO
 * @ingroup KratosCore
 * @brief This class defines an interface to the GiDPost library in order to provide GiD compliant I/O functionality
 * @tparam TGaussPointContainer The gauss point container considered
 * @tparam TMeshContainer The mesh container considered
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Pooyan Dadvand
 */
template<class TGaussPointContainer = GidGaussPointsContainer, class TMeshContainer = GidMeshContainer>
class KRATOS_API(KRATOS_CORE) GidIO
    : public GidIOBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GidIO
    KRATOS_CLASS_POINTER_DEFINITION(GidIO);

    /// Base class definition
    typedef IO BaseType;

    /// Containers definition
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// Integration method definition
    typedef GeometryData::IntegrationMethod IntegrationMethodType;

    /// Geometry family definition
    typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. Single stream IO constructor
     */
    GidIO(
        const std::string& rDatafilename,
        const GiD_PostMode Mode,
        const MultiFileFlag UseMultipleFilesFlag,
        const WriteDeformedMeshFlag WriteDeformedFlag,
        const WriteConditionsFlag WriteConditions,
        const bool InitializeGaussPointContainers=true
         ) : mResultFileName(rDatafilename),
        mMeshFileName(rDatafilename),
        mWriteDeformed(WriteDeformedFlag),
        mWriteConditions(WriteConditions),
        mUseMultiFile(UseMultipleFilesFlag),
        mMode(Mode)
    {
        mResultFileOpen = false;
        mMeshFileOpen = false;

        InitializeResultFile(mResultFileName);
        SetUpMeshContainers();
        if (InitializeGaussPointContainers) {
            SetUpGaussPointContainers();
        }

        GidIOBase& r_gid_io_base = GidIOBase::GetInstance();

        if (r_gid_io_base.GetData() == 0){
            GiD_PostInit();
        }
        GiD_PostSetFormatReal("%g");
        r_gid_io_base.SetData(r_gid_io_base.GetData() + 1);
    }

    ///Destructor.
    ~GidIO() override
    {
        if ( mResultFileOpen ) {
            GiD_fClosePostResultFile( mResultFile );
            mResultFileOpen = false;
        }

        GidIOBase& r_gid_io_base = GidIOBase::GetInstance();

        r_gid_io_base.SetData(r_gid_io_base.GetData() - 1);

        if (r_gid_io_base.GetData() == 0) {
            GiD_PostDone();
        }
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates the mesh containers for all different element types.
     * @note The containers are not filled yet in here!
     */
    void SetUpMeshContainers()
    {
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,
                                          GiD_Hexahedra, "Kratos_Hexahedra3D20_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D27,
                                          GiD_Hexahedra, "Kratos_Hexahedra3D27_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,
                                          GiD_Hexahedra, "Kratos_Hexahedra3D8_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Prism3D15,
                                          GiD_Prism, "Kratos_Prism3D15_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Prism3D6,
                                          GiD_Prism, "Kratos_Prism3D6_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral2D4_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral2D8_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral2D9_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral3D4_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral3D8_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral3D9_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,
                                          GiD_Tetrahedra, "Kratos_Tetrahedra3D10_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,
                                          GiD_Tetrahedra, "Kratos_Tetrahedra3D4_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle2D3,
                                          GiD_Triangle, "Kratos_Triangle2D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle2D6,
                                          GiD_Triangle, "Kratos_Triangle2D6_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle3D3,
                                          GiD_Triangle, "Kratos_Triangle3D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle3D6,
                                          GiD_Triangle, "Kratos_Triangle3D6_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line2D2,
                                          GiD_Linear, "Kratos_Line2D2_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line3D2,
                                          GiD_Linear, "Kratos_Line3D2_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line2D3,
                                          GiD_Linear, "Kratos_Line2D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line3D3,
                                          GiD_Linear, "Kratos_Line3D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Point2D,
                                          GiD_Point, "Kratos_Point2D_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Point3D,
                                          GiD_Point, "Kratos_Point3D_Mesh" ) );


    }//SetUpMeshContainers

    /**
     * @brief Creates the gauss point containers for all different element types.
     * @note The containers are not filled yet in here!
     */
    virtual void SetUpGaussPointContainers()
    {
        //elements with 1 gauss point
        std::vector<int> gp_indices(1, 0);

        //case Triangle with 1 gauss point
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri1_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Triangle, GiD_Triangle, 1, gp_indices ) );
        //case Quadrilateral with 1 gauss point
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad1_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Quadrilateral, GiD_Quadrilateral, 1, gp_indices ) );
        //case Tetrahedra with 1 gauss point
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet1_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Tetrahedra, GiD_Tetrahedra, 1, gp_indices ) );
        //case Hexahedra with 1 gauss point
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex1_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Hexahedra, GiD_Hexahedra, 1, gp_indices ) );
        //case Prism with 1 gauss point
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism1_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 1, gp_indices ) );
        //case Linear with 1 gauss point
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin1_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 1, gp_indices ) );

        //elements with 2 gauss points
        gp_indices.resize(2);
        gp_indices[1] = 1;
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin2_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 2, gp_indices ) );


        //elements with 3 gauss points
        gp_indices.resize(3);
        //case Triangle with 3 gauss points
        gp_indices[2] = 2;
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri3_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Triangle, GiD_Triangle, 3, gp_indices ) );
        //case Linear with 3 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin3_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 3, gp_indices ) );

        //elements with 4 gauss points
        gp_indices.resize(4);
        gp_indices[3] = 3;

        //case Linear with 4 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin3_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 4, gp_indices ) );
        //case Quadrilateral with 4 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad4_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Quadrilateral, GiD_Quadrilateral, 4, gp_indices ) );
        //case Tetrahedra with 4 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet4_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Tetrahedra, GiD_Tetrahedra, 4, gp_indices ) );
        //case Triangle with 4 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri4_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Triangle, GiD_Triangle, 4, gp_indices ) );
        gp_indices[0] = 1;
        gp_indices[1] = 2;
        gp_indices[2] = 3;
        gp_indices[3] = 4;
        //case Tetrahedra with 5 gauss points (4 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet5_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Tetrahedra, GiD_Tetrahedra, 5, gp_indices ) );
        //case Tetrahedra with 11 gauss points (4 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet11_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Tetrahedra, GiD_Tetrahedra, 11, gp_indices ) );


        //elements with 5 gauss points
        gp_indices.resize(5);
        gp_indices[0] = 0;
        gp_indices[1] = 1;
        gp_indices[2] = 2;
        gp_indices[3] = 3;
        gp_indices[4] = 4;
        //case Linear with 5 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin5_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 5, gp_indices ) );

        //case Tetrahedra with 10 gauss points (4 gauss points will be created for GiD)
        gp_indices.resize(10);
        gp_indices[5] = 5;
        gp_indices[6] = 6;
        gp_indices[7] = 7;
        gp_indices[8] = 8;
        gp_indices[9] = 9;
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet10_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Tetrahedra, GiD_Tetrahedra, 10, gp_indices ) );

        //elements with 6 gauss points
        gp_indices.resize(6);
        gp_indices[0] = 0;
        gp_indices[1] = 1;
        gp_indices[2] = 2;
        gp_indices[3] = 3;
        gp_indices[4] = 4;
        gp_indices[5] = 5;
        //case Triangle with 6 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri6_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Triangle, GiD_Triangle, 6, gp_indices ) );
        //case Prism with 6 Gauss Points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism6_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 6, gp_indices ) );

        /* START: Adding manually the custom prism */
        //case Prism with 2 Gauss Points (6 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism2_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 2, gp_indices ) );
        //case Prism with 3 Gauss Points (6 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism3_element_gp",
                                                                  GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 3, gp_indices ) );
        //case Prism with 5 Gauss Points (6 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism5_element_gp",
                                                                  GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 5, gp_indices ) );
        //case Prism with 7 Gauss Points (6 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism7_element_gp",
                                                                  GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 7, gp_indices ) );
        //case Prism with 11 Gauss Points (6 gauss points will be created for GiD)
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism11_element_gp",
                                                                  GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 11, gp_indices ) );
        /* END: Adding manually the custom prism */

        //elements with 7 gauss points
        gp_indices.resize(7);
        gp_indices[6] = 6;
        //case Linear with 7 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin7_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 7, gp_indices ) );
        //elements with 8 gauss points
        gp_indices.resize(8);
        gp_indices[7] = 7;
        //case Hexahedra with 8 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex8_element_gp",
                                                                  GeometryData::KratosGeometryFamily::Kratos_Hexahedra, GiD_Hexahedra, 8, gp_indices ) );

        //elements with 9 gauss points
        gp_indices.resize(9);
        gp_indices[8] = 8;
        //case Linear with 9 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin9_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 9, gp_indices ) );
        //case Prism with 9 Gauss Points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "prism9_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Prism, GiD_Prism, 9, gp_indices ) );
        // case quadrilateral with 9 Gauss Points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad9_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Quadrilateral, GiD_Quadrilateral, 9, gp_indices ) );

        //elements with 11 gauss points
        gp_indices.resize(11);
        gp_indices[10] = 10;
        //case Linear with 11 gauss points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "lin11_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Linear, GiD_Linear, 11, gp_indices ) );

        //elements with 27 gauss points
        gp_indices.resize(27);
        gp_indices[0] = 0;
        gp_indices[8] = 1;
        gp_indices[1] = 2;
        gp_indices[11] = 3;
        gp_indices[20] = 4;
        gp_indices[9] = 5;
        gp_indices[3] = 6;
        gp_indices[10] = 7;
        gp_indices[2] = 8;
        gp_indices[12] = 9;
        gp_indices[21] = 10;
        gp_indices[13] = 11;
        gp_indices[24] = 12;
        gp_indices[26] = 13;
        gp_indices[22] = 14;
        gp_indices[15] = 15;
        gp_indices[23] = 16;
        gp_indices[14] = 17;
        gp_indices[4] = 18;
        gp_indices[16] = 19;
        gp_indices[5] = 20;
        gp_indices[19] = 21;
        gp_indices[25] = 22;
        gp_indices[17] = 23;
        gp_indices[7] = 24;
        gp_indices[18] = 25;
        gp_indices[6] = 26;
        //case Hexahedra with 27 Gauss Points
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex27_element_gp",
                                            GeometryData::KratosGeometryFamily::Kratos_Hexahedra, GiD_Hexahedra, 27, gp_indices ) );

    }//SetUpGaussPointContainers


    /// General GidIO related functions ///

    /**
     * @brief
     * @todo To be removed
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
        mResultFileName = rResultFileName;
    }

    /**
     * TODO: check whether this is still necessary!
     */
    void  CloseResultFile()
    {
        if ( mResultFileOpen )
        {
            GiD_fClosePostResultFile( mResultFile );
            mResultFileOpen = false;
        }
    }

    /**
     * TODO: check whether this is still necessary!
     */
    void Flush()
    {
        GiD_fFlushPostFile( mResultFile );
    }


    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        return "gid io";
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///result functions
    /**
     * This has to be called for each solution step BEFORE any results
     * (on nodes and on gauss points) is written
     * @param SolutionTag the current solution step (i.e. time)
     * @param conditions_flag states whether results should also be written on conditions
     */
    virtual void InitializeResults(
        const double name,
        const MeshType& rThisMesh
        )
    {
        if ( mMode == GiD_PostAscii && ! mResultFileOpen )
        {
            std::stringstream file_name;
            if ( mUseMultiFile == SingleFile ) {
                file_name << mResultFileName << ".post.res";
            } else {
                file_name << mResultFileName << std::setprecision(12) << "_" << name << ".post.res";
            }
            mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
            mResultFileOpen = true;
        }
        //initializing gauss points containers
        if ( mWriteConditions != WriteConditionsOnly )
        {
            for ( auto element_iterator = rThisMesh.ElementsBegin(); element_iterator != rThisMesh.ElementsEnd(); ++element_iterator )
            {
                for ( auto it = mGidGaussPointContainers.begin();  it != mGidGaussPointContainers.end(); it++ )
                {
                    if ( it->AddElement( element_iterator ) )
                        break;
                }

            }
        }

        if ( mWriteConditions == WriteConditionsFlag::WriteConditions || mWriteConditions == WriteConditionsOnly )
            for ( auto conditions_iterator = rThisMesh.ConditionsBegin(); conditions_iterator != rThisMesh.ConditionsEnd(); conditions_iterator++ )
            {
                for ( auto it =  mGidGaussPointContainers.begin();  it != mGidGaussPointContainers.end(); it++ )
                {

                    if ( it->AddCondition( conditions_iterator ) )
                        break;
                }
            }

        // Writing gauss points definitions
        for ( auto it = mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {
            it->WriteGaussPoints(mResultFile);
        }
    }
    /**
     * This has to be called for each solution step AFTER all the results
     * have been written
     */
    void FinalizeResults()
    {
        if ( mUseMultiFile == MultipleFiles || mMode == GiD_PostAscii )
        {
            GiD_fClosePostResultFile( mResultFile );
            mResultFileOpen = false;
        }
        //resetting gauss point containers
        for ( auto it =
                    mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {
            it->Reset();
        }
    }

    ///functions for writing nodal results

    ///////////////////////////////////////////////////////////////////////
    //////                  HISTORICAL DATABASE BLOCK                 /////
    ///////////////////////////////////////////////////////////////////////
     /**
     * writes nodal results for variables of type bool
     */
    void WriteNodalResults(
        Variable<bool> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin();  it_node != rNodes.end() ; ++it_node)
            GiD_fWriteScalar( mResultFile, it_node->Id(), static_cast<double>(it_node->GetSolutionStepValue(rVariable,
                                SolutionStepNumber)) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

    ///functions for writing nodal results
    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResults(
        Variable<double> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin();  it_node != rNodes.end() ; ++it_node)
            GiD_fWriteScalar( mResultFile, it_node->Id(), it_node->GetSolutionStepValue(rVariable,
                                SolutionStepNumber) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type int
     */
    void WriteNodalResults(
        Variable<int> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
            GiD_fWriteScalar( mResultFile, it_node->Id(), it_node->GetSolutionStepValue(rVariable,
                                SolutionStepNumber) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }


    /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResults(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult(mResultFile,(char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Vector,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for (auto it_node = rNodes.begin();
                it_node != rNodes.end() ; ++it_node)
        {
            const array_1d<double, 3>& temp = it_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteVector( mResultFile, it_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 components can be printed)
     */
    void WriteNodalResults(
        Variable<Vector> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Matrix,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const Vector& temp_vector = it_node->FastGetSolutionStepValue(rVariable,
                                    SolutionStepNumber);
            if (temp_vector.size() ==3 )
                GiD_fWrite2DMatrix(mResultFile, it_node->Id(), temp_vector[0], temp_vector[1], temp_vector[2]);
            else if (temp_vector.size() == 6 )
                GiD_fWrite3DMatrix( mResultFile, it_node->Id(), temp_vector[0], temp_vector[1], temp_vector[2],
                                    temp_vector[3], temp_vector[4], temp_vector[5] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Matrix
     */
    void WriteNodalResults(
        Variable<Matrix> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Matrix,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const Matrix& temp_matrix = it_node->GetSolutionStepValue(rVariable,
                    SolutionStepNumber);
            if (temp_matrix.size1() ==3 && temp_matrix.size2() ==3)
            {
                GiD_fWrite3DMatrix( mResultFile,  it_node->Id(), temp_matrix(0,0), temp_matrix(1,1),
                                    temp_matrix(2,2), temp_matrix(0,1), temp_matrix(1,2),
                                    temp_matrix(0,2) );
            }
            else if (temp_matrix.size1() ==2 && temp_matrix.size2() ==2)
            {
                GiD_fWrite2DMatrix( mResultFile, it_node->Id(), temp_matrix(0,0), temp_matrix(1,1), temp_matrix(0,1));
            }

            else if (temp_matrix.size1() ==1 && temp_matrix.size2() ==3)
            {

                GiD_fWrite3DMatrix( mResultFile, it_node->Id(), temp_matrix(0,0), temp_matrix(0,1), 0.00,
                                    temp_matrix(0,2), 0.00, 0.00);
            }
            else if (temp_matrix.size1() ==1 && temp_matrix.size2() ==6)
            {
                GiD_fWrite3DMatrix( mResultFile, it_node->Id(), temp_matrix(0,0), temp_matrix(0,1), temp_matrix(0,2),
                                    temp_matrix(0,3), temp_matrix(0,4), temp_matrix(0,5) );
            }
            //it_node->GetValue(rVariable) = temp_matrix;

        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    void WriteLocalAxesOnNodes(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_LocalAxes,
                            GiD_OnNodes, NULL, NULL, 0, NULL );

        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const array_1d<double, 3>& temp = it_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteLocalAxes( mResultFile, it_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    ///////////////////////////////////////////////////////////////////////
    //////                 NON- HISTORICAL DATABASE BLOCK             /////
    ///////////////////////////////////////////////////////////////////////

   /**
    * Writes nodal flags
    */
    void WriteNodalFlags(
        const Kratos::Flags& rFlag,
        const std::string& rFlagName,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rFlagName.c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            GiD_fWriteScalar( mResultFile, it_node->Id(),  static_cast<double>(it_node->Is(rFlag)));
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type bool
     */
    void WriteNodalResultsNonHistorical(
        Variable<bool> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
            GiD_fWriteScalar( mResultFile, it_node->Id(), static_cast<double>(it_node->GetValue(rVariable)) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    ///functions for writing nodal results
    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResultsNonHistorical(
        Variable<double> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
            GiD_fWriteScalar( mResultFile, it_node->Id(), it_node->GetValue(rVariable) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type int
     */
    void WriteNodalResultsNonHistorical(
        Variable<int> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Scalar,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
            GiD_fWriteScalar( mResultFile, it_node->Id(), it_node->GetValue(rVariable) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResultsNonHistorical(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult(mResultFile,(char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Vector,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const array_1d<double, 3>& temp = it_node->GetValue( rVariable);
            GiD_fWriteVector( mResultFile, it_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 components can be printed)
     */
    void WriteNodalResultsNonHistorical(
        Variable<Vector> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Matrix,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const Vector& temp_vector = it_node->GetValue(rVariable);
            if (temp_vector.size() ==3 )
                GiD_fWrite2DMatrix(mResultFile, it_node->Id(), temp_vector[0], temp_vector[1], temp_vector[2]);
            else if (temp_vector.size() == 6 )
                GiD_fWrite3DMatrix( mResultFile, it_node->Id(), temp_vector[0], temp_vector[1], temp_vector[2],
                                    temp_vector[3], temp_vector[4], temp_vector[5] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Matrix
     */
    void WriteNodalResultsNonHistorical(
        Variable<Matrix> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_Matrix,
                            GiD_OnNodes, NULL, NULL, 0, NULL );
        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const Matrix& temp_matrix = it_node->GetValue(rVariable);
            if (temp_matrix.size1() ==3 && temp_matrix.size2() ==3)
            {
                GiD_fWrite3DMatrix( mResultFile,  it_node->Id(), temp_matrix(0,0), temp_matrix(1,1),
                                    temp_matrix(2,2), temp_matrix(0,1), temp_matrix(1,2),
                                    temp_matrix(0,2) );
            }
            else if (temp_matrix.size1() ==2 && temp_matrix.size2() ==2)
            {
                GiD_fWrite2DMatrix( mResultFile, it_node->Id(), temp_matrix(0,0), temp_matrix(1,1), temp_matrix(0,1));
            }

            else if (temp_matrix.size1() ==1 && temp_matrix.size2() ==3)
            {

                GiD_fWrite3DMatrix( mResultFile, it_node->Id(), temp_matrix(0,0), temp_matrix(0,1), 0.00,
                                    temp_matrix(0,2), 0.00, 0.00);
            }
            else if (temp_matrix.size1() ==1 && temp_matrix.size2() ==6)
            {
                GiD_fWrite3DMatrix( mResultFile, it_node->Id(), temp_matrix(0,0), temp_matrix(0,1), temp_matrix(0,2),
                                    temp_matrix(0,3), temp_matrix(0,4), temp_matrix(0,5) );
            }

        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    void WriteLocalAxesOnNodesNonHistorical(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                            SolutionTag, GiD_LocalAxes,
                            GiD_OnNodes, NULL, NULL, 0, NULL );

        for (auto it_node = rNodes.begin(); it_node != rNodes.end() ; ++it_node)
        {
            const array_1d<double, 3>& temp = it_node->GetSolutionStepValue( rVariable);
            GiD_fWriteLocalAxes( mResultFile, it_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    ///mesh writing functions
    /**
     * opens a new mesh group
     */
    void InitializeMesh( const double name )
    {
        if ( mUseMultiFile == MultipleFiles )
        {
            if ( mMode == GiD_PostAscii && ! mMeshFileOpen )
            {
                std::stringstream file_name;
                file_name << std::setprecision(12) << mMeshFileName << "_" << name << ".post.msh";
                mMeshFile = GiD_fOpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                mMeshFileOpen = true;
            }
            if ( (mMode == GiD_PostBinary || mMode == GiD_PostHDF5) && ! mResultFileOpen )
            {
                std::stringstream file_name;
                file_name << std::setprecision(12) << mResultFileName << "_" << name << ".post.bin";
                if ( ! mResultFileOpen )
                {
                    mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
                    mResultFileOpen = true;
                }
                mMeshFile = mResultFile;
            }
        }
        if ( mUseMultiFile == SingleFile )
        {
            if ( (mMode == GiD_PostBinary || mMode == GiD_PostHDF5) && ! mResultFileOpen )
            {
                std::stringstream file_name;
                file_name << mResultFileName << ".post.bin";
                mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
                if ( mResultFile == 0) //error handler can not be zero
                {
                    std::stringstream buffer;
                    buffer << "error opening results file:" << "/" <<  file_name.str()   << "/";
                    KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
                }
                mResultFileOpen = true;
                mMeshFile = mResultFile;
            }
            if ( mMode == GiD_PostAscii && ! mMeshFileOpen )
            {
                std::stringstream file_name;
                file_name << mMeshFileName << ".post.msh";
                mMeshFile = GiD_fOpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                mMeshFileOpen = true;
            }
        }
    }

    /**
     * closes a mesh group
     */
    void FinalizeMesh()
    {
        if ( mUseMultiFile == MultipleFiles && mMode == GiD_PostAscii )
        {
            GiD_fClosePostMeshFile(mMeshFile);
            mMeshFileOpen = false;
        }
        if ( mUseMultiFile == SingleFile && mMode == GiD_PostAscii )
        {
            GiD_fClosePostMeshFile(mMeshFile);
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
    void WriteNodeMesh(MeshType& rThisMesh) override
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        GiD_fBeginMesh(mMeshFile,  "Kratos Mesh",GiD_3D,GiD_Point,1);
        GiD_fBeginCoordinates(mMeshFile);
        for ( auto node_iterator = rThisMesh.NodesBegin();
                node_iterator != rThisMesh.NodesEnd();
                ++node_iterator)
        {
            if ( mWriteDeformed == WriteUndeformed )
                GiD_fWriteCoordinates(mMeshFile, node_iterator->Id(), node_iterator->X0(),
                                        node_iterator->Y0(), node_iterator->Z0() );
            else if ( mWriteDeformed == WriteDeformed )
                GiD_fWriteCoordinates(mMeshFile, node_iterator->Id(), node_iterator->X(),
                                        node_iterator->Y(), node_iterator->Z() );
            else
                KRATOS_ERROR << "Undefined WriteDeformedMeshFlag" << std::endl;
        }
        GiD_fEndCoordinates(mMeshFile);
        int nodes_id[1];
        GiD_fBeginElements(mMeshFile);

        for ( auto node_iterator = rThisMesh.NodesBegin(); node_iterator != rThisMesh.NodesEnd(); ++node_iterator)
        {
            nodes_id[0] = node_iterator->Id();
            GiD_fWriteElement(mMeshFile,node_iterator->Id(), nodes_id);
        }
        GiD_fEndElements(mMeshFile);
        GiD_fEndMesh(mMeshFile);

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }

    void WriteSphereMesh(const MeshType& rThisMesh)
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        GiD_fBeginMesh(mMeshFile, "Kratos Mesh",GiD_3D,GiD_Sphere,1);
        GiD_fBeginCoordinates(mMeshFile);
        for ( auto node_iterator = rThisMesh.NodesBegin(); node_iterator != rThisMesh.NodesEnd(); ++node_iterator)
        {
            if ( mWriteDeformed == WriteUndeformed )
                GiD_fWriteCoordinates( mMeshFile, node_iterator->Id(), node_iterator->X0(),
                                        node_iterator->Y0(), node_iterator->Z0() );
            else if ( mWriteDeformed == WriteDeformed )
                GiD_fWriteCoordinates( mMeshFile, node_iterator->Id(), node_iterator->X(),
                                        node_iterator->Y(), node_iterator->Z() );
            else
                KRATOS_ERROR << "Undefined WriteDeformedMeshFlag" << std::endl;
        }
        GiD_fEndCoordinates( mMeshFile );

        GiD_fBeginElements( mMeshFile );

        // DEM variables
        const Variable<int>& particle_material = KratosComponents<Variable<int>>::Get("PARTICLE_MATERIAL");
        const Variable<double>& radius = KratosComponents<Variable<double>>::Get("RADIUS");

        for ( auto element_iterator = rThisMesh.ElementsBegin();  element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
        {
            const unsigned int node_id = element_iterator->GetGeometry()[0].Id();
            GiD_fWriteSphereMat(mMeshFile, node_id, node_id, element_iterator->GetGeometry()[0].FastGetSolutionStepValue(radius), element_iterator->GetGeometry()[0].FastGetSolutionStepValue(particle_material));
        }
        GiD_fEndElements( mMeshFile );
        GiD_fEndMesh( mMeshFile);

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }

    void WriteCircleMesh(const MeshType& rThisMesh)
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");
        GiD_fBeginMesh(mMeshFile, "Kratos Mesh",GiD_2D,GiD_Circle,1);
        GiD_fBeginCoordinates(mMeshFile);
        for ( auto node_iterator = rThisMesh.NodesBegin(); node_iterator != rThisMesh.NodesEnd(); ++node_iterator)
        {
            if ( mWriteDeformed == WriteUndeformed )
                GiD_fWriteCoordinates( mMeshFile, node_iterator->Id(), node_iterator->X0(),
                                        node_iterator->Y0(), node_iterator->Z0() );
            else if ( mWriteDeformed == WriteDeformed )
                GiD_fWriteCoordinates( mMeshFile, node_iterator->Id(), node_iterator->X(),
                                        node_iterator->Y(), node_iterator->Z() );
            else
                KRATOS_ERROR << "Undefined WriteDeformedMeshFlag" << std::endl;
        }
        GiD_fEndCoordinates( mMeshFile );
        int nodes_id[1];
        GiD_fBeginElements( mMeshFile );

        double nx = 0.0;
        double ny = 0.0;
        double nz = 1.0;

        // DEM variables
        const Variable<int>& particle_material = KratosComponents<Variable<int>>::Get("PARTICLE_MATERIAL");
        const Variable<double>& radius = KratosComponents<Variable<double>>::Get("RADIUS");

        for ( auto node_iterator = rThisMesh.NodesBegin();  node_iterator != rThisMesh.NodesEnd(); ++node_iterator)
        {
            nodes_id[0] = node_iterator->Id();
            GiD_fWriteCircleMat(mMeshFile, node_iterator->Id(), nodes_id[0], node_iterator->FastGetSolutionStepValue(radius), nx, ny, nz, node_iterator->FastGetSolutionStepValue(particle_material));
        }
        GiD_fEndElements( mMeshFile );
        GiD_fEndMesh( mMeshFile);
        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }

    void WriteClusterMesh(const MeshType& rThisMesh)
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        GiD_fBeginMesh(mMeshFile, "Kratos Mesh",GiD_3D,GiD_Cluster,1);
        GiD_fBeginCoordinates(mMeshFile);
        for ( auto node_iterator = rThisMesh.NodesBegin(); node_iterator != rThisMesh.NodesEnd(); ++node_iterator)
        {
            if ( mWriteDeformed == WriteUndeformed )
                GiD_fWriteCoordinates( mMeshFile, node_iterator->Id(), node_iterator->X0(),
                                        node_iterator->Y0(), node_iterator->Z0() );
            else if ( mWriteDeformed == WriteDeformed )
                GiD_fWriteCoordinates( mMeshFile, node_iterator->Id(), node_iterator->X(),
                                        node_iterator->Y(), node_iterator->Z() );
            else
                KRATOS_ERROR << "Undefined WriteDeformedMeshFlag" << std::endl;
        }
        GiD_fEndCoordinates( mMeshFile );

        GiD_fBeginElements( mMeshFile );

        // DEM variables
        const Variable<int>& particle_material = KratosComponents<Variable<int>>::Get("PARTICLE_MATERIAL");

        for ( auto element_iterator = rThisMesh.ElementsBegin();  element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
        {
            const unsigned int node_id = element_iterator->GetGeometry()[0].Id();
            GiD_fWriteClusterMat(mMeshFile, node_id, node_id, element_iterator->GetGeometry()[0].FastGetSolutionStepValue(particle_material));
        }
        GiD_fEndElements( mMeshFile );
        GiD_fEndMesh( mMeshFile);

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }//WriteClusterMesh


    /**
     * This is a multi-purpose function that writes arbitrary meshes of elements
     * and conditions in either deformed or undeformed state
     * @param rThisMesh the current mesh to be written
     * @param deformed_flag states whether the mesh should be written in deformed configuration
     * @param conditions_flag states whether conditions should also be written
     */
    void WriteMesh(MeshType& rThisMesh) override
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        if ( mWriteConditions != WriteConditionsOnly )
        {
            for ( auto element_iterator = rThisMesh.ElementsBegin(); element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
                for ( auto it = mGidMeshContainers.begin(); it != mGidMeshContainers.end(); it++ )
                    if ( it->AddElement( element_iterator ) )
                        break;
        }
        if ( mWriteConditions == WriteConditionsFlag::WriteConditions || mWriteConditions == WriteConditionsOnly )
        {
            for ( auto conditions_iterator = rThisMesh.ConditionsBegin(); conditions_iterator != rThisMesh.ConditionsEnd(); conditions_iterator++ )
                for ( auto it = mGidMeshContainers.begin(); it != mGidMeshContainers.end(); it++ )
                    if ( it->AddCondition( conditions_iterator ) )
                        break;
        }

        for ( auto it = mGidMeshContainers.begin(); it != mGidMeshContainers.end(); it++ )
        {
            it->FinalizeMeshCreation();
            if ( mWriteDeformed == WriteDeformed )
                it->WriteMesh(mMeshFile,true);
            else if ( mWriteDeformed == WriteUndeformed )
                it->WriteMesh(mMeshFile,false);
            else
                KRATOS_ERROR << "Undefined WriteDeformedMeshFlag" << std::endl;

            it->Reset();
        }

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }

    ///functions for printing results on gauss points

    /**
    * @brief Writes elemental and conditional flags
    * @param rFlag the flag
    * @param rFlagName the given flag name
    * @param rModelPart the current model part
    */
    void PrintFlagsOnGaussPoints(
        const Kratos::Flags& rFlag,
        const std::string& rFlagName,
        const ModelPart& rModelPart,
        const double SolutionTag
        )
    {
        Timer::Start("Writing Results");

        for ( auto it =  mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintFlagsResults( mResultFile, rFlag, rFlagName, rModelPart, SolutionTag );
        }

        Timer::Stop("Writing Results");
    }

    /**
     * Prints variables of type int on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<bool>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag, const int ValueIndex = 0 )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( auto it = mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintResults( mResultFile, rVariable, rModelPart, SolutionTag, ValueIndex );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type int on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<int>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( auto it = mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintResults( mResultFile, rVariable, rModelPart, SolutionTag, ValueIndex );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<double>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( auto it = mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintResults( mResultFile, rVariable, rModelPart, SolutionTag, ValueIndex );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<array_1d<double,3> >& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( auto it = mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintResults(  mResultFile, rVariable, rModelPart, SolutionTag, ValueIndex );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<Vector>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        )
    {
        KRATOS_TRY;
        Timer::Start("Writing Results");

        for ( auto it = mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintResults(  mResultFile, rVariable, rModelPart, SolutionTag, ValueIndex );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<Matrix>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        )
    {
        KRATOS_TRY;
        Timer::Start("Writing Results");

        for ( auto it = mGidGaussPointContainers.begin(); it != mGidGaussPointContainers.end(); it++ ) {
            it->PrintResults(  mResultFile, rVariable, rModelPart, SolutionTag, ValueIndex );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

protected:
    /**
     * File names
     */
    std::string mResultFileName;
    std::string mMeshFileName;

    GiD_FILE mMeshFile;
    GiD_FILE mResultFile;

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

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) GidIO<GidGaussPointsContainer,GidMeshContainer>;

///@}
///@name Input and output
///@{

inline std::ostream& operator << (std::ostream& rOStream, const GidIO<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

///@}
