//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_MONOLITHIC_WALL_CONDITION_H
#define KRATOS_MONOLITHIC_WALL_CONDITION_H

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

/// Implements a wall condition for the monolithic formulation.
/**
  It is intended to be used in combination with ASGS and VMS elements or their derived classes
  and the ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent time scheme, which supports
  slip conditions.
  This condition will add a wall stress term to all nodes identified with IS_STRUCTURE!=0.0 (in the
  non-historic database, that is, assigned using Node.SetValue()). This stress term is determined
  according to the wall distance provided as Y_WALL.
  @see ASGS2D,ASGS3D,VMS,ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent
 */
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) MonolithicWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonolithicWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(MonolithicWallCondition);

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    MonolithicWallCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    MonolithicWallCondition(IndexType NewId,
                            const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    MonolithicWallCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    MonolithicWallCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    MonolithicWallCondition(MonolithicWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~MonolithicWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    MonolithicWallCondition & operator=(MonolithicWallCondition const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new MonolithicWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<MonolithicWallCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< MonolithicWallCondition >(NewId, pGeom, pProperties);
    }
    
    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );
        
        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());
        
        return pNewCondition;
    }

    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
    }

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
    }

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
    }



    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                            ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType RHS;
        this->CalculateLocalVelocityContribution(rDampingMatrix,RHS,rCurrentProcessInfo);
    }



    /// Calculate wall stress term for all nodes with IS_STRUCTURE != 0.0
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override;


    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) override
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
            if(VISCOSITY.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","");
            if(IS_STRUCTURE.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"IS_STRUCTURE Key is 0. Check if the application was correctly registered.","");
            if(Y_WALL.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument,"Y_WALL Key is 0. Check if the application was correctly registered.","");
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
    void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo) override;


    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList,
                            ProcessInfo& CurrentProcessInfo) override;


    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& Values,
                                           int Step = 0) override
    {
        const SizeType LocalSize = (TDim + 1) * TNumNodes;
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
            for (unsigned int d = 0; d < TDim; ++d)
                Values[LocalIndex++] = rVelocity[d];
            Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
        }
    }


    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& Values,
                                            int Step = 0) override
    {
        const SizeType LocalSize = (TDim + 1) * TNumNodes;
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
            for (unsigned int d = 0; d < TDim; ++d)
                Values[LocalIndex++] = rVelocity[d];
            Values[LocalIndex++] = 0.0; // No value on pressure positions
        }
    }


    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                             std::vector<array_1d<double, 3 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;



    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                             std::vector<double>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;


    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                             std::vector<array_1d<double, 6 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                             std::vector<Vector>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;


    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                             std::vector<Matrix>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;


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
        buffer << "MonolithicWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MonolithicWallCondition";
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

    /// Commpute the wall stress and add corresponding terms to the system contributions.
    /**
      @param rLocalMatrix Local system matrix
      @param rLocalVector Local right hand side
      */
    virtual void ApplyWallLaw(MatrixType& rLocalMatrix,
                      VectorType& rLocalVector,
		      ProcessInfo& rCurrentProcessInfo)
    {
        GeometryType& rGeometry = this->GetGeometry();
        const size_t BlockSize = TDim + 1;
        const double NodalFactor = 1.0 / double(TDim);

        double area = NodalFactor * rGeometry.DomainSize();
        // DomainSize() is the way to ask the geometry's length/area/volume (whatever is relevant for its dimension) without asking for the number of spatial dimensions first

        for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
        {
            const NodeType& rConstNode = rGeometry[itNode];
            const double y = rConstNode.GetValue(Y_WALL); // wall distance to use in stress calculation
            if( y > 0.0 && rConstNode.GetValue(IS_STRUCTURE) != 0.0 )
            {
                array_1d<double,3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
                const array_1d<double,3>& VelMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                Vel -= VelMesh;
                const double Ikappa = 1.0/0.41; // inverse of Von Karman's kappa
                const double B = 5.2;
                const double limit_yplus = 10.9931899; // limit between linear and log regions

                const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
                const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);

                double wall_vel = 0.0;
                for (size_t d = 0; d < TDim; d++)
                {
                    wall_vel += Vel[d]*Vel[d];
                }
                wall_vel = sqrt(wall_vel);

                if (wall_vel > 1e-12) // do not bother if velocity is zero
                {

                    // linear region
                    double utau = sqrt(wall_vel * nu / y);
                    double yplus = y * utau / nu;

                    // log region
                    if (yplus > limit_yplus)
                    {

                        // wall_vel / utau = 1/kappa * log(yplus) + B
                        // this requires solving a nonlinear problem:
                        // f(utau) = utau*(1/kappa * log(y*utau/nu) + B) - wall_vel = 0
                        // note that f'(utau) = 1/kappa * log(y*utau/nu) + B + 1/kappa

                        unsigned int iter = 0;
                        double dx = 1e10;
                        const double tol = 1e-6;
                        double uplus = Ikappa * log(yplus) + B;

                        while(iter < 100 && fabs(dx) > tol * utau)
                        {
                            // Newton-Raphson iteration
                            double f = utau * uplus - wall_vel;
                            double df = uplus + Ikappa;
                            dx = f/df;

                            // Update variables
                            utau -= dx;
                            yplus = y * utau / nu;
                            uplus = Ikappa * log(yplus) + B;
                            ++iter;
                        }
                        if (iter == 100)
                        {
                            std::cout << "Warning: wall condition Newton-Raphson did not converge. Residual is " << dx << std::endl;
                        }
                    }
                    const double Tmp = area * utau * utau * rho / wall_vel;
                    for (size_t d = 0; d < TDim; d++)
                    {
                        size_t k = itNode*BlockSize+d;
                        rLocalVector[k] -= Vel[d] * Tmp;
                        rLocalMatrix(k,k) += Tmp;
                    }
                }
            }
        }
    }


    void CalculateNormal(array_1d<double,3>& An );


    void ApplyNeumannCondition(MatrixType &rLocalMatrix, VectorType &rLocalVector);


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

}; // Class MonolithicWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream,
                                  MonolithicWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MonolithicWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_MONOLITHIC_WALL_CONDITION_H
