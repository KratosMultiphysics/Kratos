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

#ifndef KRATOS_NAVIER_STOKES_WALL_CONDITION_H
#define KRATOS_NAVIER_STOKES_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
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
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NavierStokesWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NavierStokesWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(NavierStokesWallCondition);

    struct ConditionDataStruct
    {
        array_1d<double, 3> Normal;     // Condition normal
        double wGauss;                  // Gauss point weight
        array_1d<double, TNumNodes> N;  // Gauss point shape functions values
    };

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new         // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);condition
      */
    NavierStokesWallCondition(IndexType NewId = 0):Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    NavierStokesWallCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    NavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    NavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    NavierStokesWallCondition(NavierStokesWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    virtual ~NavierStokesWallCondition() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    NavierStokesWallCondition & operator=(NavierStokesWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new NavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer(new NavierStokesWallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Create a new NavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    virtual Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared< NavierStokesWallCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @return a Pointer to the new element
     */
    virtual Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
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
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;

        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_gauss;
        array_1d<double,MatrixSize> rhs_gauss;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Compute condition unit normal vector
        this->CalculateNormal(data.Normal); //this already contains the area
        const double A = std::sqrt(data.Normal[0]*data.Normal[0]+data.Normal[1]*data.Normal[1]+data.Normal[2]*data.Normal[2]);
        data.Normal /= A;

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = rGeom.DeterminantOfJacobian(igauss, GeometryData::GI_GAUSS_2);
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
    virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;

        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

        // Compute condition normal
        this->CalculateNormal(data.Normal); //this already contains the area
        const double A = std::sqrt(data.Normal[0]*data.Normal[0]+data.Normal[1]*data.Normal[1]+data.Normal[2]*data.Normal[2]);
        data.Normal /= A;

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = rGeom.DeterminantOfJacobian(igauss, GeometryData::GI_GAUSS_2);
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointLHSContribution(lhs_local, data);

            noalias(rLeftHandSideMatrix) += lhs_local;
        }

        KRATOS_CATCH("")
    }

    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector: reference to the RHS matrix
     * @param rCurrentProcessInfo: reference to the ProcessInfo (unused)
     */
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        // Loop on gauss points
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Compute condition normal
        this->CalculateNormal(data.Normal); //this already contains the area
        const double A = std::sqrt(data.Normal[0]*data.Normal[0]+data.Normal[1]*data.Normal[1]+data.Normal[2]*data.Normal[2]);
        data.Normal /= A;

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = rGeom.DeterminantOfJacobian(igauss, GeometryData::GI_GAUSS_2);
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointRHSContribution(rhs_local, data);

            noalias(rRightHandSideVector) += rhs_local;
        }

        KRATOS_CATCH("")
    }



    // virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix,
    //                         ProcessInfo& rCurrentProcessInfo)
    // {
    //     VectorType RHS;
    //     this->CalculateLocalVelocityContribution(rDampingMatrix,RHS,rCurrentProcessInfo);
    // }



    /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.0
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    // virtual void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
    //         VectorType& rRightHandSideVector,
    //         ProcessInfo& rCurrentProcessInfo);


    /// Check that all data required by this condition is available and reasonable
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered
            if(VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
            if(MESH_VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","");
            if(ACCELERATION.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check if the application was correctly registered.","");
            if(PRESSURE.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
            if(DENSITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
            if(DYNAMIC_VISCOSITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_VISCOSITY Key is 0. Check if the application was correctly registered.","");
            // if(IS_STRUCTURE.Key() == 0)
            //     KRATOS_THROW_ERROR(std::invalid_argument,"IS_STRUCTURE Key is 0. Check if the application was correctly registered.","");
            // if(Y_WALL.Key() == 0)
            //     KRATOS_THROW_ERROR(std::invalid_argument,"Y_WALL Key is 0. Check if the application was correctly registered.","");
            if(EXTERNAL_PRESSURE.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"EXTERNAL_PRESSURE Key is 0. Check if the application was correctly registered.","");

            // Checks on nodes
            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {
                if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing ACCELERATION variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing EXTERNAL_PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                   this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                   this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
                if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                    KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
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
    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo);

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
        std::stringstream buffer;
        buffer << "NavierStokesWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NavierStokesWallCondition";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    // 1D line shape functions values at Gauss points
    // void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,2,2>& Ncontainer)
    // {
    //     const double val1 = 0.5*(1+std::sqrt(1/3))
    //     const double val2 = 0.5*(1-std::sqrt(1/3))
    //     Ncontainer(0,0) = val1; Ncontainer(0,1) = val2;
    //     Ncontainer(1,0) = val2; Ncontainer(1,1) = val1;
    // }
    //
    // // 2D triangle shape functions values at Gauss points
    // void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,3,3>& Ncontainer)
    // {
    //     const double one_sixt = 1.0/6.0;
    //     const double two_third = 2.0/3.0;
    //     Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
    //     Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
    //     Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    // }

    /// Commpute the wall stress and add corresponding terms to the system contributions.
    /**
      @param rLocalMatrix Local system matrix
      @param rLocalVector Local right hand side
      */
    // virtual void ApplyWallLaw(MatrixType& rLocalMatrix,
    //                   VectorType& rLocalVector,
	// 	      ProcessInfo& rCurrentProcessInfo)
    // {
    //     GeometryType& rGeometry = this->GetGeometry();
    //     const size_t BlockSize = TDim + 1;
    //     const double NodalFactor = 1.0 / double(TDim);
    //
    //     double area = NodalFactor * rGeometry.DomainSize();
    //     // DomainSize() is the way to ask the geometry's length/area/volume (whatever is relevant for its dimension) without asking for the number of spatial dimensions first
    //
    //     for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
    //     {
    //         const NodeType& rConstNode = rGeometry[itNode];
    //         const double y = rConstNode.GetValue(Y_WALL); // wall distance to use in stress calculation
    //         if( y > 0.0 && rConstNode.GetValue(IS_STRUCTURE) != 0.0 )
    //         {
    //             array_1d<double,3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
    //             const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
    //             Vel -= VelMesh;
    //             const double Ikappa = 1.0/0.41; // inverse of Von Karman's kappa
    //             const double B = 5.2;
    //             const double limit_yplus = 10.9931899; // limit between linear and log regions
    //
    //             const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
    //             const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);
    //
    //             double wall_vel = 0.0;
    //             for (size_t d = 0; d < TDim; d++)
    //             {
    //                 wall_vel += Vel[d]*Vel[d];
    //             }
    //             wall_vel = sqrt(wall_vel);
    //
    //             if (wall_vel > 1e-12) // do not bother if velocity is zero
    //             {
    //
    //                 // linear region
    //                 double utau = sqrt(wall_vel * nu / y);
    //                 double yplus = y * utau / nu;
    //
    //                 // log region
    //                 if (yplus > limit_yplus)
    //                 {
    //
    //                     // wall_vel / utau = 1/kappa * log(yplus) + B
    //                     // this requires solving a nonlinear problem:
    //                     // f(utau) = utau*(1/kappa * log(y*utau/nu) + B) - wall_vel = 0
    //                     // note that f'(utau) = 1/kappa * log(y*utau/nu) + B + 1/kappa
    //
    //                     unsigned int iter = 0;
    //                     double dx = 1e10;
    //                     const double tol = 1e-6;
    //                     double uplus = Ikappa * log(yplus) + B;
    //
    //                     while(iter < 100 && fabs(dx) > tol * utau)
    //                     {
    //                         // Newton-Raphson iteration
    //                         double f = utau * uplus - wall_vel;
    //                         double df = uplus + Ikappa;
    //                         dx = f/df;
    //
    //                         // Update variables
    //                         utau -= dx;
    //                         yplus = y * utau / nu;
    //                         uplus = Ikappa * log(yplus) + B;
    //                         ++iter;
    //                     }
    //                     if (iter == 100)
    //                     {
    //                         std::cout << "Warning: wall condition Newton-Raphson did not converge. Residual is " << dx << std::endl;
    //                     }
    //                 }
    //                 const double Tmp = area * utau * utau * rho / wall_vel;
    //                 for (size_t d = 0; d < TDim; d++)
    //                 {
    //                     size_t k = itNode*BlockSize+d;
    //                     rLocalVector[k] -= Vel[d] * Tmp;
    //                     rLocalMatrix(k,k) += Tmp;
    //                 }
    //             }
    //         }
    //     }
    // }

    void CalculateNormal(array_1d<double,3>& An);

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ConditionDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
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

}; // Class NavierStokesWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream,
                                  NavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_NAVIER_STOKES_WALL_CONDITION_H
