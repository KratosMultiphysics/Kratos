/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Anna Bauer
//                  Thomas Oberbichler
//                  Tobias Teschemacher
*/

#if !defined(KRATOS_IGA_MEMBRANE_ELEMENT_H_INCLUDED)
#define KRATOS_IGA_MEMBRANE_ELEMENT_H_INCLUDED

// System includes
#include "iga_application_variables.h"
//#include "iga_application.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"
#include "utilities/math_utils.h"

// Project includes
#include "iga_base_element.h"

#include "custom_utilities/geometry_utilities/iga_curve_on_surface_utilities.h"

namespace Kratos
{

class IgaMembraneElement
    : public IgaBaseElement
{
protected: 

    struct MetricVariables 
    {
        Vector gab; // covariant metric
        Vector gab_con; // contravariant metric
        Vector curvature; //
        Matrix J; //Jacobian
        double  detJ;
        Vector g1; //base vector 1
        Vector g2; //base vector 2
        Vector g3; //base vector 3
        double dA; //differential area
        Matrix H; //Hessian
        Matrix Q; //Transformation matrix Q from contravariant to cartesian basis
        Matrix T; //Transformation matrix T from contravariant to local cartesian basis
       
        MetricVariables(const unsigned int& Dimension)
        {
            gab = ZeroVector(Dimension);
            gab_con = ZeroVector(Dimension);

            curvature = ZeroVector(Dimension);

            J = ZeroMatrix(Dimension, Dimension);
            detJ = 1.0;

            g1 = ZeroVector(Dimension);
            g2 = ZeroVector(Dimension);
            g3 = ZeroVector(Dimension);

            dA = 1.0;

            Matrix H = ZeroMatrix(3, 3);
            Matrix Q = ZeroMatrix(3, 3);
            Matrix T = ZeroMatrix(3, 3);
        }

    }; 

    MetricVariables mInitialMetric = MetricVariables(3); 

    struct ConstitutiveVariables //aus surface_base_element.h übernommen
    {
        Vector E; //strain
        Vector S; //stress
        Matrix D; //constitutive matrix

        /**
        * The default constructor
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        ConstitutiveVariables(const unsigned int& StrainSize)
        {
            E = ZeroVector(StrainSize);
            S = ZeroVector(StrainSize);
            D = ZeroMatrix(StrainSize, StrainSize);
        }
    };

struct SecondVariations
    {
        Matrix B11;
        Matrix B22;
        Matrix B12;

        /**
        * The default constructor
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        SecondVariations(const int& mat_size)
        {
            B11 = ZeroMatrix(mat_size, mat_size);
            B22 = ZeroMatrix(mat_size, mat_size);
            B12 = ZeroMatrix(mat_size, mat_size);
        }
    };

struct PrestresstransVariables
    {
        Matrix Tpre;
        
        PrestresstransVariables(const unsigned int& Dimension)
        {
           Matrix Tpre = ZeroMatrix(3, 3);
        }
    };

struct transPreVariables
    {
        Vector S_prestress_result;

         transPreVariables(const unsigned int& Dimension)
        {
            Vector S_prestress_result = ZeroVector(Dimension);
        }

    };
    

public:
    KRATOS_CLASS_POINTER_DEFINITION( IgaMembraneElement );

    IgaMembraneElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : IgaBaseElement(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    IgaMembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : IgaBaseElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    IgaMembraneElement() : IgaBaseElement() {};

    ~IgaMembraneElement() override
    {
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    
	/**
	* @brief Creates a new element
	* @param NewId The Id of the new created element
	* @param pGeom The pointer to the geometry of the element
	* @param pProperties The pointer to property
	* @return The pointer to the created element
	*/

	Element::Pointer Create(
		IndexType NewId,
		GeometryType::Pointer pGeom,
		PropertiesType::Pointer pProperties
	) const override
	{
		return Kratos::make_shared<IgaMembraneElement>(
			NewId, pGeom, pProperties);
	};    

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) override;

    void PrintInfo(std::ostream& rOStream) const override;

/**
    * Calculate a double Variable on the Element Constitutive Law
    * @param rVariable: The variable we want to get
    * @param rOutput: The values obtained int the integration points
    * @param rCurrentProcessInfo: the current process info instance
    */
    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    void CalculateStresses(
        Vector& rStresses,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculatePresstressTensor(
        Vector& rPrestressTensor,
        MetricVariables& rMetric); 

private:

ConstitutiveLaw::Pointer mConstitutiveLaw; 

void CalculateMetric( MetricVariables& metric ); 

void IgaMembraneElement::CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& metric);

void IgaMembraneElement::CalculateSecondVariationStrain(
        SecondVariations& rSecondVariationsStrain,
        const MetricVariables& rMetric);

void IgaMembraneElement::CalculateStrain(
        Vector& StrainVector,
        Vector& gab,
        Vector& gab0);

void CalculateConstitutiveVariables( 
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure);

void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight );

 void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight);
    
void CalculateSecondVariationStrainMembrane(
        SecondVariations& rSecondVariationsStrain,
        const MetricVariables& rMetric);

void CalculateTransformationmatrixPrestress(
        const MetricVariables& metric,
        PrestresstransVariables& Prestresstrans
        );

};

} // namespace Kratos

#endif // !defined(KRATOS_IGA_MEMBRANE_ELEMENT_H_INCLUDED)
