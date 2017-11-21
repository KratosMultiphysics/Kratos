//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#ifndef KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_WALL_CONDITION_H
#define KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Implements a wall condition for the Navier-Stokes monolithic formulation.
/**
  It is intended to be used in combination with ASGS Navier-Stokes symbolic elements or their
  derived classes and the ResidualBasedIncrementalUpdateStaticSchemeSlip time scheme, which supports
  slip conditions.
  @see NavierStokes,EmbeddedNavierStokes,ResidualBasedIncrementalUpdateStaticSchemeSlip
 */
template <unsigned int TDim, unsigned int TNumNodes = TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) EmbeddedAusasNavierStokesWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedAusasNavierStokesWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedAusasNavierStokesWallCondition);

    struct ConditionDataStruct
    {
        double wGauss;                 // Gauss point weight
        // double charVel;                // Problem characteristic velocity (used in the outlet inflow prevention)
        // double delta;                  // Non-dimensional positive sufficiently small constant (used in the outlet inflow prevention)
        array_1d<double, 3> Normal;    // Condition normal
        array_1d<double, TNumNodes> N; // Gauss point shape functions values

        bounded_matrix<double, TNumNodes, TDim> v; // Current step velocity
    };

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new         // Struct to pass around the data
        ConditionDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);condition
      */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId = 0):Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    EmbeddedAusasNavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    EmbeddedAusasNavierStokesWallCondition(EmbeddedAusasNavierStokesWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~EmbeddedAusasNavierStokesWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    EmbeddedAusasNavierStokesWallCondition & operator=(EmbeddedAusasNavierStokesWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new EmbeddedAusasNavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new EmbeddedAusasNavierStokesWallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Create a new EmbeddedAusasNavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return boost::make_shared< EmbeddedAusasNavierStokesWallCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    /// Calculates the LHS and RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix: reference to the LHS matrix
     * @param rRightHandSideVector: reference to the RHS matrix
     * @param rCurrentProcessInfo: reference to the ProcessInfo (unused)
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;
        this->FillElementData(data);

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_gauss;
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_gauss;

        // LHS and RHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Store the outlet inflow prevention constants in the data structure
        // data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        // const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
        // data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussPtsJDet = ZeroVector(NumGauss);
        rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        // Loop on gauss points
        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointRHSContribution(rhs_gauss, data);
            ComputeGaussPointLHSContribution(lhs_gauss, data);

            noalias(rLeftHandSideMatrix) += lhs_gauss;
            noalias(rRightHandSideVector) += rhs_gauss;
        }

        KRATOS_CATCH("")
    }

    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix: reference to the LHS matrix
     * @param rCurrentProcessInfo: reference to the ProcessInfo (unused)
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;
        this->FillElementData(data);

        // Allocate memory needed
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_gauss;

        // LHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize, MatrixSize);

        // Gauss point information
        GeometryType &rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType &IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussPtsJDet = ZeroVector(NumGauss);
        rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        // Loop on gauss points
        for (unsigned int igauss = 0; igauss < NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointLHSContribution(lhs_gauss, data);

            noalias(rLeftHandSideMatrix) += lhs_gauss;
        }

        KRATOS_CATCH("")
    }

    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector: reference to the RHS matrix
     * @param rCurrentProcessInfo: reference to the ProcessInfo (unused)
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;
        this->FillElementData(data);

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_gauss;

        // RHS contributions initialization
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Store the outlet inflow prevention constants in the data structure
        // data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        // const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
        // data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussPtsJDet = ZeroVector(NumGauss);
        rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointRHSContribution(rhs_gauss, data);

            noalias(rRightHandSideVector) += rhs_gauss;
        }

        KRATOS_CATCH("")
    }


    /// Condition check
    /**
     * @param rCurrentProcessInfo: reference to the ProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const GeometryType& rGeom = this->GetGeometry();

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered
            if(VELOCITY.Key() == 0)
                KRATOS_ERROR << "VELOCITY Key is 0. Check if the application was correctly registered.";
            if(MESH_VELOCITY.Key() == 0)
                KRATOS_ERROR << "MESH_VELOCITY Key is 0. Check if the application was correctly registered.";
            if(ACCELERATION.Key() == 0)
                KRATOS_ERROR << "ACCELERATION Key is 0. Check if the application was correctly registered.";
            if(PRESSURE.Key() == 0)
                KRATOS_ERROR << "PRESSURE Key is 0. Check if the application was correctly registered.";
            if(DENSITY.Key() == 0)
                KRATOS_ERROR << "DENSITY Key is 0. Check if the application was correctly registered.";
            if(DYNAMIC_VISCOSITY.Key() == 0)
                KRATOS_ERROR << "DYNAMIC_VISCOSITY Key is 0. Check if the application was correctly registered.";
            if(EXTERNAL_PRESSURE.Key() == 0)
                KRATOS_ERROR << "EXTERNAL_PRESSURE Key is 0. Check if the application was correctly registered.";

            // Checks on nodes
            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<rGeom.size(); ++i)
            {
                if(rGeom[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_ERROR << "Missing VELOCITY variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(PRESSURE) == false)
                    KRATOS_ERROR << "Missing PRESSURE variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_ERROR << "Missing MESH_VELOCITY variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(ACCELERATION) == false)
                    KRATOS_ERROR << "Missing ACCELERATION variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                    KRATOS_ERROR << "Missing EXTERNAL_PRESSURE variable on solution step data for node " << rGeom[i].Id();
                if(rGeom[i].HasDofFor(VELOCITY_X) == false || rGeom[i].HasDofFor(VELOCITY_Y) == false || rGeom[i].HasDofFor(VELOCITY_Z) == false)
                    KRATOS_ERROR << "Missing VELOCITY component degree of freedom on node " << rGeom[i].Id();
                if(rGeom[i].HasDofFor(PRESSURE) == false)
                    KRATOS_ERROR << "Missing PRESSURE component degree of freedom on node " << rGeom[i].Id();
            }

            return Check;
        }

        KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo) override;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EmbeddedAusasNavierStokesWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmbeddedAusasNavierStokesWallCondition";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateNormal(array_1d<double,3>& An);

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ConditionDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);
    
    void ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);
    // void ComputeRHSOutletInflowContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

    // Auxiliar function to fill the element data structure
    void FillElementData(ConditionDataStruct &rData)
    {
        const GeometryType& rGeom = this->GetGeometry();

        // Compute condition normal
        this->CalculateNormal(rData.Normal);    // This already contains the area
        const double A = norm_2(rData.Normal);  // Compute the element area
        rData.Normal /= A;                      // Divide by the area to get the unit normal vector

        // Fill the nodal velocity array
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double, 3> &vel = rGeom[i].FastGetSolutionStepValue(VELOCITY);

            for (unsigned int k = 0; k < TDim; k++)
            {
                rData.v(i, k) = vel[k];
            }
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class EmbeddedAusasNavierStokesWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_EMBEDDED_AUSAS_NAVIER_STOKES_WALL_CONDITION_H
