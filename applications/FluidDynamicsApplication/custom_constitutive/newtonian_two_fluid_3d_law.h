//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

#if !defined (KRATOS_NEWTONIAN_TWO_FLUID_3D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_TWO_FLUID_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "newtonian_3d_law.h"
#include "custom_utilities/element_size_calculator.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

/**
 * Defines a Newtonian constitutive law in 3D.
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY (read from the nodes!!)
 * 2) C_SMAGORINSKY
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NewtonianTwoFluid3DLaw : public Newtonian3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef std::size_t             SizeType;
    
    /**
     * Counted pointer of NewtonianTwoFluid3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NewtonianTwoFluid3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewtonianTwoFluid3DLaw()
        : Newtonian3DLaw()  {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override {
        return Kratos::make_shared<NewtonianTwoFluid3DLaw>(*this);
    }

    /**
     * Copy constructor.
     */
    NewtonianTwoFluid3DLaw (const NewtonianTwoFluid3DLaw& rOther)
        : Newtonian3DLaw()  {}


    /**
     * Destructor.
     */
    ~NewtonianTwoFluid3DLaw() override {}


    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override {
        return "NewtonianTwoFluid3DLaw";
    }

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
    
    
    double ComputeEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override {

        double viscosity;
        EvaluateInPoint(viscosity, DYNAMIC_VISCOSITY, rParameters);
        
        const Properties& prop = rParameters.GetMaterialProperties();

        if (prop.Has(C_SMAGORINSKY)) {
            const double csmag = prop[C_SMAGORINSKY];
            if (csmag > 0.0) {
                double density;
                EvaluateInPoint(density, DENSITY, rParameters);
                const double strain_rate = EquivalentStrainRate(rParameters);
                const BoundedMatrix<double, 4, 3>& DN_DX = rParameters.GetShapeFunctionsDerivatives();
                const double elem_size = ElementSizeCalculator<3,4>::GradientsElementSize(DN_DX);
                double length_scale = csmag * elem_size;
                length_scale *= length_scale;
                viscosity += 2.0*length_scale * strain_rate * density;
            }
        }
        return viscosity;

    }

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
    
    void EvaluateInPoint(double& rResult,
        const Variable<double>& rVariable,
        ConstitutiveLaw::Parameters& rParameters) const {
        
        const SizeType nnodes = 4;

        const GeometryType& geom = rParameters.GetElementGeometry();

        const array_1d<double,nnodes>& N = rParameters.GetShapeFunctionsValues();

        //compute sign of distance on gauss point
        double dist = 0.0;
        for (unsigned int i = 0; i < nnodes; i++)
            dist += N[i] * geom[i].FastGetSolutionStepValue(DISTANCE);
        
        SizeType navg = 0; //number of nodes on the same side as the gauss point
        double value = 0.0;
        for (unsigned int i = 0; i < nnodes; i++) {
            if ( dist * geom[i].FastGetSolutionStepValue(DISTANCE) > 0.0) {
                navg += 1;
                value += geom[i].FastGetSolutionStepValue(rVariable);
            }
        }

        rResult = value/navg;

    }

    double EquivalentStrainRate(ConstitutiveLaw::Parameters& rParameters) const {
        // Calculate Symetric gradient (Voigt notation)
        const SizeType nnodes = 4;
        const GeometryType& geom = rParameters.GetElementGeometry();
        const Matrix& DN_DX = rParameters.GetShapeFunctionsDerivatives();

        array_1d<double,6> S(6,0.0);
        for (unsigned int n = 0; n < nnodes; ++n)
        {
            const array_1d<double,3>& rVel = geom[n].FastGetSolutionStepValue(VELOCITY);
            S[0] += DN_DX(n,0)*rVel[0];
            S[1] += DN_DX(n,1)*rVel[1];
            S[2] += DN_DX(n,2)*rVel[2];
            S[3] += DN_DX(n,2)*rVel[1] + DN_DX(n,1)*rVel[2];
            S[4] += DN_DX(n,2)*rVel[0] + DN_DX(n,0)*rVel[2];
            S[5] += DN_DX(n,1)*rVel[0] + DN_DX(n,0)*rVel[1];
        }

        // Norm of symetric gradient (cross terms don't get the 2)
        return std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + 2.*S[2]*S[2] + S[3]*S[3] + S[4]*S[4] + S[5]*S[5]);
    }


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

    void load(Serializer& rSerializer) override {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

}; // Class NewtonianTwoFluid3DLaw
}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_LAW_NODAL_PROPERTIES_3D_H_INCLUDED  defined 
