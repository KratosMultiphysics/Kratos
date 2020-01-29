//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_MOVE_SHALLOW_PARTICLES_PROCESS_H_INCLUDED
#define KRATOS_MOVE_SHALLOW_PARTICLES_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_utilities/convection_operator.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
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

/// Short class definition.
/** Detail class definition.
*/
template<size_t TDim>
class KRATOS_API(SHALLOW_WATER_APPLICATION) MoveShallowParticlesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    /// Pointer definition of MoveShallowParticlesProcess
    KRATOS_CLASS_POINTER_DEFINITION(MoveShallowParticlesProcess<TDim>);

    ///@}
    ///@name  Enum's
    ///@{

    enum ProjectionType {Incompressible, Compressible};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MoveShallowParticlesProcess(
        ModelPart& rModelPart,
        ModelPart& rParticles,
        Variable<array_1d<double,3>>& rVectorVariable,
        Variable<double>& rScalarVariable,
        Parameters Settings);

    /// Destructor.
    virtual ~MoveShallowParticlesProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function will seed the initial particles
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will update the mesh info
     */
    void ExecuteBeforeSolutionLoop() override;

    /**
     * @brief This function will move the particles and transfer the info to the eulerian mesh
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This function will update the particles with the eulerian mesh solution
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function is called to verify that the input is correct
     */
    int Check() override;

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
    virtual std::string Info() const override {return "MoveShallowParticlesProcess";}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << Info();}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    ModelPart& mrModelPart;
    ModelPart& mrParticles;
    Variable<array_1d<double,3>>& mrMomentumVariable;
    Variable<double>& mrMassVariable;
    double mDryTreshold;
    size_t mLastParticleId;
    ConvectionOperator<TDim> mConvectionOperator;
    GeometryData::IntegrationMethod mQuadratureOrder;
    ProjectionType mProjectionType;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeVariables();

    void InitializeParticle(NodeType& rParticle, Element& pElement);

    void ComputeMeanSize();

    void ComputeMeanVelocity();

    void SeedInitialParticles();

    void ComputeNumberOfParticlesPerElement();

    void MoveParticles();

    void TransferLagrangianToEulerian();

    void TransferEulerianToLagrangian();

    void StandardProjection();

    void StandardParticlesUpdate();

    void ProjectionWithMassConservation();

    void ParticlesUpdateWithMassConservation();

    void GetShapeFunctionsValues(Vector& rN, const GeometryType& rGeom, const NodeType& rParticle);

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

    /// Assignment operator.
    // MoveShallowParticlesProcess& operator=(MoveShallowParticlesProcess const& rOther);

    /// Copy constructor.
    // MoveShallowParticlesProcess(MoveShallowParticlesProcess const& rOther);


    ///@}

}; // Class MoveShallowParticlesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<size_t TDim>
inline std::istream& operator >> (std::istream& rIStream,
                MoveShallowParticlesProcess<TDim>& rThis);

/// output stream function
template<size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                const MoveShallowParticlesProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOVE_SHALLOW_PARTICLES_PROCESS_H_INCLUDED  defined
