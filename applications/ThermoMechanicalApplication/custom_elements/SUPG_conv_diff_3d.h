//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//                   Riccardo Rossi
//


#if !defined(KRATOS_SUPG_CONV_DIFF_3D_INCLUDED )
#define  KRATOS_SUPG_CONV_DIFF_3D_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"


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

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// ASGS, Incompressible fluid, Variational multi scale method, Quasi-static subscales, Implicit method.
/**
ASGS is an abriviation for Algebraic Sub-Grid Scale element. It is implemented to solve
Implicitly the NS equations in a variotionally consistant sub-grid scale methid. It also has the OSS swith
to use Orthogonal Sub Scales to use impose explicity the orthogonality condition on subscales' estimation.
The "Dynamic_Tau" swith allows the use of "Dt", time step, in calculation of Tau.
This element just work with Monolithic schemes like "monolithic_solver_eulerian" or "monolithic_solver_lagranigan".
The detailed description of the formulation could be fined in
   "Stabilized finite element approximation of transient incompressible flows using orthogonal subscales, Comput. Methods Appl. Mech. Engrg. 191 (2002) 4295?4321"

*/
class SUPGConvDiff3D
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Fluid2DASGS
    KRATOS_CLASS_POINTER_DEFINITION(SUPGConvDiff3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    SUPGConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry);
    SUPGConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~SUPGConvDiff3D();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// The DOF's are VELOCITY_X, VELOCITY_Y and PRESSURE
    /**
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
    */
    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;



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
        return "SUPGConvDiff3D #" ;
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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
    SUPGConvDiff3D() : Element() {}

    ///@}
    ///@name Protected Operations
    ///@{


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

    // private:

    /// To Calculate tau
    /**
     * @return tau: multiplied by Residual of the momentum equation
    */
    virtual void CalculateTau(array_1d<double, 3 >& ms_adv_vel, double& tau,const double& K, const double time, const double area, const ProcessInfo& rCurrentProcessInfo);
    virtual void CalculateArtifitialViscosity(double& art_visc,
						      BoundedMatrix<double, 4, 3 > DN_DX,
						      array_1d<double, 3 > ms_vel_gauss,
						      const Variable<double>& temperature,
						      const double volume,
						      const double scaled_K);

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
//         ASGS2D() : Element()
//         {
//         }

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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
    //Fluid2DASGS& operator=(const Fluid2DASGS& rOther);

    /// Copy constructor.
    //Fluid2DASGS(const Fluid2DASGS& rOther);


    ///@}

}; // Class Fluid2DASGS

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    Fluid2DASGS& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_ASGS_2D_H_INCLUDED  defined


