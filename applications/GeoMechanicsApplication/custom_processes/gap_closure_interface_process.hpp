// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_PERIODIC_INTERFACE_PROCESS )
#define  KRATOS_GEO_PERIODIC_INTERFACE_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/math_utilities.hpp"
#include "utilities/math_utils.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class GapClosureInterfaceProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(GapClosureInterfaceProcess);

    ///definition of the geometry type with given NodeType
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    GapClosureInterfaceProcess( ModelPart& model_part,
                              Parameters rParameters ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "consider_gap_closure": true,
                "gap_width_threshold": 0.01
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mGapWidthThreshold = rParameters["gap_width_threshold"].GetDouble();
        mConsiderGapClosure = rParameters["consider_gap_closure"].GetBool();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~GapClosureInterfaceProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the GapClosureInterfaceProcess algorithms.
    void Execute() override
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mConsiderGapClosure)
        {
            const int nelements = mr_model_part.GetMesh(0).Elements().size();

            if (nelements > 0)
            {
                ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();

                std::vector<bool> SetToDeactive(nelements);
                #pragma omp parallel for
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    SetToDeactive[k] = IsGapCreated(it);
                }


                // Activation/deactivation of the existing parts:
                // ( User must specify each part through the interface)
                #pragma omp parallel for
                for (int k = 0; k < nelements; ++k)
                {
                    if (SetToDeactive[k])
                    {
                        ModelPart::ElementsContainerType::iterator it = el_begin + k;
                        it->Set(ACTIVE, false);
                    } 
                    else
                    {
                        ModelPart::ElementsContainerType::iterator it = el_begin + k;
                        it->Set(ACTIVE, true);
                    }
                }
            }
        }
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GapClosureInterfaceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GapClosureInterfaceProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    double mGapWidthThreshold;
    bool mConsiderGapClosure;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    GapClosureInterfaceProcess& operator=(GapClosureInterfaceProcess const& rOther);

    bool IsGapCreated(const ModelPart::ElementsContainerType::iterator &it)
    {
        const GeometryType& rGeom = it->GetGeometry();
        const SizeType Dim = rGeom.WorkingSpaceDimension();

        if (Dim == N_DIM_2D)
        {
            return IsGapCreated2D4N(it);
        } 
        else if (Dim == N_DIM_3D)
        {
            const SizeType NumNodes = rGeom.PointsNumber();
            if (NumNodes == 6)
            {
                return IsGapCreated3D6N(it);
            }
            else if (NumNodes == 8)
            {
                return IsGapCreated3D8N(it);
            }
            else
            {
                KRATOS_ERROR << "undefined number of nodes in gap_closure_interface_process"
                             << std::endl;
            }
        }
        else
        {
            KRATOS_ERROR << "undefined dimension in gap_closure_interface_process"
                         << std::endl;
         }
    }

//----------------------------------------------------------------------------------------------------
    bool IsGapCreated2D4N(const ModelPart::ElementsContainerType::iterator &it)
    {

        KRATOS_TRY;
        const std::vector<std::vector<int>> IndexSide{{0,1},{2,3}};
        const unsigned int nSides = 2;
        const GeometryType& rGeom = it->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const unsigned int nPointSide = NumNodes/nSides;


        // find normal vectors
        const array_1d<double,3> OutOfPlane{0.0, 0.0, 1.0};

        double GapSize = 0.0;
        for (unsigned int iSide=0; iSide < nSides; ++iSide)
        {
            // find normal vectors
            array_1d<double,3> TangentVector;
            noalias(TangentVector) = rGeom.GetPoint(IndexSide[iSide][1]) - rGeom.GetPoint(IndexSide[iSide][0]);
            array_1d<double,3> normalVector;
            MathUtils<double>::UnitCrossProduct(normalVector, OutOfPlane, TangentVector);

            double GapSide = 0.0;
            for (unsigned int iPoint=0; iPoint < nPointSide; ++iPoint)
            {
                const array_1d<double,3> &displacement = rGeom[IndexSide[iSide][iPoint]].FastGetSolutionStepValue(DISPLACEMENT, 0);
                double GapPoint = inner_prod(normalVector, displacement);
                GapSide -= GapPoint;
            }
            GapSide /= double(nPointSide);
            GapSize += GapSide;
        }

        return !(GapSize < mGapWidthThreshold);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------
    bool IsGapCreated3D6N(const ModelPart::ElementsContainerType::iterator &it)
    {
        KRATOS_TRY;
        const std::vector<std::vector<int>> IndexSide{{0,1,2},{5,4,3}};
        const unsigned int nSides = 2;
        const GeometryType& rGeom = it->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const unsigned int nPointSide = NumNodes/nSides;

        double GapSize = 0.0;
        for (unsigned int iSide=0; iSide < nSides; ++iSide)
        {
            // find normal vectors
            array_1d<double,3> Vector0;
            array_1d<double,3> Vector1;

            noalias(Vector0) = rGeom.GetPoint(IndexSide[iSide][1]) - rGeom.GetPoint(IndexSide[iSide][0]);
            noalias(Vector1) = rGeom.GetPoint(IndexSide[iSide][2]) - rGeom.GetPoint(IndexSide[iSide][0]);

            array_1d<double,3> normalVector;
            MathUtils<double>::UnitCrossProduct(normalVector, Vector0, Vector1);

            double GapSide = 0.0;
            for (unsigned int iPoint=0; iPoint < nPointSide; ++iPoint)
            {
                const array_1d<double,3> &displacement = rGeom[IndexSide[iSide][iPoint]].FastGetSolutionStepValue(DISPLACEMENT, 0);
                double GapPoint = inner_prod(normalVector, displacement);
                GapSide -= GapPoint;
            }
            GapSide /= double(nPointSide);
            GapSize += GapSide;
        }

        return !(GapSize < mGapWidthThreshold);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------
    bool IsGapCreated3D8N(const ModelPart::ElementsContainerType::iterator &it)
    {
        KRATOS_TRY;

        const std::vector<std::vector<int>> IndexSide{{0,1,2,3},{7,6,5,4}};
        const unsigned int nSides = 2;
        const GeometryType& rGeom = it->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const unsigned int nPointSide = NumNodes/nSides;

        double GapSize = 0.0;
        for (unsigned int iSide=0; iSide < nSides; ++iSide)
        {
            // find normal vectors
            array_1d<double,3> Vector0;
            array_1d<double,3> Vector1;

            noalias(Vector0) = rGeom.GetPoint(IndexSide[iSide][1]) - rGeom.GetPoint(IndexSide[iSide][0]);
            noalias(Vector1) = rGeom.GetPoint(IndexSide[iSide][2]) - rGeom.GetPoint(IndexSide[iSide][0]);

            array_1d<double,3> normalVector;
            MathUtils<double>::UnitCrossProduct(normalVector, Vector0, Vector1);

            double GapSide = 0.0;
            for (unsigned int iPoint=0; iPoint < nPointSide; ++iPoint)
            {
                const array_1d<double,3> &displacement = rGeom[IndexSide[iSide][iPoint]].FastGetSolutionStepValue(DISPLACEMENT, 0);
                double GapPoint = inner_prod(normalVector, displacement);
                GapSide -= GapPoint;
            }
            GapSide /= double(nPointSide);
            GapSize += GapSide;
        }

        return !(GapSize < mGapWidthThreshold);

        KRATOS_CATCH( "" )
    }


    /// Copy constructor.
    //GapClosureInterfaceProcess(GapClosureInterfaceProcess const& rOther);

}; // Class GapClosureInterfaceProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GapClosureInterfaceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GapClosureInterfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_PERIODIC_INTERFACE_PROCESS defined */
