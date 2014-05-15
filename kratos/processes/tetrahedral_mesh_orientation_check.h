#ifndef KRATOS_TETRAHEDRAL_MESH_ORIENTATION_CHECK_H
#define KRATOS_TETRAHEDRAL_MESH_ORIENTATION_CHECK_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/math_utils.h"

namespace Kratos
{


/// Check a triangular or tetrahedral mesh to ensure that local connectivities follow the expected convention.
/** This process checks all elements to verify that their Jacobian has positive determinant and face conditions
 *  to ensure that all face normals point outwards.
 *
 *  Note that, as a side result of the procedure used, nodal normals (not normalized) are computed and stored on
 *  solution step data NORMAL.
 */
class TetrahedralMeshOrientationCheck: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedralMeshOrientationCheck);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for TetrahedralMeshOrientationCheck Process
    /**
     * @param rModelPart The model part to check.
     * @param ThrowErrors If true, an error will be thrown if the input model part contains malformed elements or conditions.
     */
    TetrahedralMeshOrientationCheck(ModelPart& rModelPart,
                                    bool ThrowErrors):
        mrModelPart(rModelPart),
        mThrowErrors(ThrowErrors)
    {
    }

    /// Destructor.
    virtual ~TetrahedralMeshOrientationCheck() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    virtual void Execute()
    {
        KRATOS_TRY;

        // Check that required nodal variables are available
        if(NORMAL.Key() == 0)
            KRATOS_ERROR(std::invalid_argument,"NORMAL Key is 0. There is some issue with the registering of variables.","");

        if ( (mrModelPart.NodesBegin() )->SolutionStepsDataHas(NORMAL) == false )
            KRATOS_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data for node ",mrModelPart.NodesBegin()->Id());


        // Initialize normals as zero
        array_1d<double,3> Zero(3,0.0);
        for(ModelPart::NodeIterator itNode =  mrModelPart.NodesBegin(); itNode != mrModelPart.NodesEnd(); itNode++)
        {
            noalias(itNode->FastGetSolutionStepValue(NORMAL)) = Zero;
        }

        // Main loop for elements
        unsigned int ElemSwitchCount = 0;

        for (ModelPart::ElementIterator itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); itElem++)
        {
            ElementType::GeometryType& rGeom = itElem->GetGeometry();
            GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

            if (GeoType == GeometryData::Kratos_Tetrahedra3D4  || GeoType == GeometryData::Kratos_Triangle2D3)
            {
                bool Switched = this->Orient(rGeom);
                if (Switched)
                    ElemSwitchCount++;

                this->NormalContribution(rGeom);
            }
        }

        // Generate output message, throw error if necessary
        std::stringstream OutMsg;
        if (ElemSwitchCount > 0)
        {
            OutMsg << "Mesh orientation check found " << ElemSwitchCount << " inverted elements." << std::endl;
        }
        else
        {
            OutMsg << "No inverted elements found" << std::endl;
        }


        mrModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

        // Main loop for faces
        unsigned int CondSwitchCount = 0;

        for (ModelPart::ConditionIterator itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); itCond++)
        {
            ConditionType::GeometryType& rGeom = itCond->GetGeometry();
            GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();
            array_1d<double,3> FaceNormal(3,0.0);

            if ( GeoType == GeometryData::Kratos_Triangle3D3 )
                FaceNormal3D(FaceNormal,rGeom);
            else if ( GeoType == GeometryData::Kratos_Line2D2 )
                FaceNormal2D(FaceNormal,rGeom);

            const unsigned int NumNodes = rGeom.PointsNumber();
            unsigned int MismatchCount = 0;
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                const array_1d<double,3>& rNormal = rGeom[i].FastGetSolutionStepValue(NORMAL);
                double Dot = FaceNormal[0]*rNormal[0] + FaceNormal[1]*rNormal[1] + FaceNormal[2]*rNormal[2];

                if (Dot < 0.0)
                    MismatchCount++;
            }

            // Re-orient if the face normal is not aligned with node normals
            if (MismatchCount > 0)
            {
                rGeom(0).swap(rGeom(1));
                CondSwitchCount++;
            }
        }

        if (CondSwitchCount > 0)
        {
            OutMsg << "Mesh orientation check found " << CondSwitchCount << " inverted conditions." << std::endl;
        }
        else
        {
            OutMsg << "No inverted conditions found" << std::endl;
        }


        if (mThrowErrors && (ElemSwitchCount+CondSwitchCount) > 0)
        {
            KRATOS_ERROR(std::runtime_error, OutMsg.str(), "");
        }
        else
        {
            std::cout << OutMsg.str();
        }


        KRATOS_CATCH("");
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "TetrahedralMeshOrientationCheck";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TetrahedralMeshOrientationCheck";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    const bool mThrowErrors;

    ///@}
    ///@name Private Operations
    ///@{

    bool Orient(Geometry< Node<3> >& rGeom)
    {
        const unsigned int PointIndex = 0;
        const GeometryData::IntegrationMethod Method = GeometryData::GI_GAUSS_1;

        // Re-orient the element if needed
        double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);
        if (DetJ < 0.0)
        {
            // swap two nodes to change orientation
            rGeom(0).swap(rGeom(1));
            return true;
        }
        else
            return false;
    }

    void NormalContribution(Geometry< Node<3> >&rGeom)
    {
        const unsigned int NumNodes = rGeom.PointsNumber();
        const unsigned int Dim = rGeom.WorkingSpaceDimension();

        const unsigned int PointIndex = 0;
        const GeometryData::IntegrationMethod Method = GeometryData::GI_GAUSS_1;
        double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);

        Geometry< Node<3> >::ShapeFunctionsGradientsType DN_DX;
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,Method);
        Matrix& rDN_DX = DN_DX[0];

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            array_1d<double,3>& rNormal = rGeom[i].FastGetSolutionStepValue(NORMAL);
            for (unsigned int d = 0; d < Dim; d++)
                rNormal[d] += DetJ*rDN_DX(i,d);
        }
    }


    void FaceNormal2D(array_1d<double,3>& An,
                      Geometry<Node<3> >& rGeometry)
    {
        An[0] =   rGeometry[1].Y() - rGeometry[0].Y();
        An[1] = - (rGeometry[1].X() - rGeometry[0].X());
        An[2] =    0.00;

    }


    void FaceNormal3D(array_1d<double,3>& An,
                      Geometry<Node<3> >& rGeometry)
    {

        array_1d<double,3> v1,v2;
        v1[0] = rGeometry[1].X() - rGeometry[0].X();
        v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
        v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

        v2[0] = rGeometry[2].X() - rGeometry[0].X();
        v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
        v2[2] = rGeometry[2].Z() - rGeometry[0].Z();

        MathUtils<double>::CrossProduct(An,v1,v2);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TetrahedralMeshOrientationCheck& operator=(TetrahedralMeshOrientationCheck const& rOther);

    /// Copy constructor.
    TetrahedralMeshOrientationCheck(TetrahedralMeshOrientationCheck const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TetrahedralMeshOrientationCheck& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TetrahedralMeshOrientationCheck& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_TETRAHEDRAL_MESH_ORIENTATION_PROCESS_H
