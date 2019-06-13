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

// System includes

// External includes

// Project includes
#include "iga_membrane_element.h"


namespace Kratos 
{

Element::Pointer IgaMembraneElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaMembraneElement>(NewId, geometry,
        pProperties);
}
void IgaMembraneElement::Initialize() 
    {
       KRATOS_TRY

       //Constitutive Law initialisation
        if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
            mConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLaw->InitializeMaterial(
                GetProperties(),
                GetGeometry(),
                GetValue(SHAPE_FUNCTION_VALUES)
                );  
        }
        else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        MetricVariables reference_metric(3);
        CalculateMetric(reference_metric);
        mInitialMetric = reference_metric; 
        
        KRATOS_CATCH("")
    }

void IgaMembraneElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    const int number_of_control_points = NumberOfNodes();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(NumberOfDofs());

    for (unsigned int i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    KRATOS_CATCH("")
}

void IgaMembraneElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    const int number_of_control_points = NumberOfNodes();

    if (rResult.size() != 3 * number_of_control_points)
        rResult.resize(3 * number_of_control_points, false);

    for (unsigned int i = 0; i < number_of_control_points; ++i) {
        const unsigned int index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }

    KRATOS_CATCH("")
}

    void IgaMembraneElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY 

        // definition of problem size
        const int number_of_nodes = GetGeometry().size();
        const int mat_size = number_of_nodes * 3;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
        }


        //reading in of shape function values and shape function derivatives
        Vector   N     = this->GetValue(SHAPE_FUNCTION_VALUES);
        Matrix  DN_De  = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        
        //Metric and Constitutive Variables
        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables_membrane(3);
        CalculateConstitutiveVariables(actual_metric,
            constitutive_variables_membrane,
            Values, 
            ConstitutiveLaw::StressMeasure_PK2);

        //Calculate B Matrix 
        Matrix BMembrane = ZeroMatrix (3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);

        // Nonlinear Deformation
        SecondVariations second_variations_strain(mat_size);
        CalculateSecondVariationStrain(
            second_variations_strain,
            actual_metric);

        double integration_weight = this->GetValue(INTEGRATION_WEIGHT) * mInitialMetric.dA;

     //*************************
        // Define Prestress
        double thickness = this->GetProperties().GetValue(THICKNESS);
               
        Vector S_prestress_nichttransformiert = this->GetProperties().GetValue(PRESTRESS)*thickness;
        
        PrestresstransVariables prestresstrans_variables(3);
        CalculateTransformationmatrixPrestress(
            actual_metric,
            prestresstrans_variables 
        );

        Vector S_prestress = prod(prestresstrans_variables.Tpre, S_prestress_nichttransformiert); 
    
        Vector S_total = constitutive_variables_membrane.S + S_prestress;
         
        if (CalculateStiffnessMatrixFlag == true) 
        {
            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(
                rLeftHandSideMatrix,
                BMembrane, 
                constitutive_variables_membrane.D, 
                integration_weight);

            // adding  non-linear-contribution to Stiffness-Matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                second_variations_strain,
                S_total, 
                integration_weight);

        }

        //RIGHT HAND SIDE 
          if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), S_total);
        }

    if(this->Id() == 0){KRATOS_WATCH(NumberOfNodes())};
    if(this->Id() == 0){KRATOS_WATCH(NumberOfDofs())};
    if(this->Id() == 0){KRATOS_WATCH(rLeftHandSideMatrix)};
    if(this->Id() == 0){KRATOS_WATCH(rRightHandSideVector)};    

        KRATOS_CATCH("");
    }
        
      
void IgaMembraneElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaMembraneElement\" #" << Id();
}


//********************************************************************
//********************************************************************
//********************************************************************

void IgaMembraneElement::CalculateMetric( 
        MetricVariables& metric
    )
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        
        //Calculate Jacobi
        IgaGeometryUtilities::CalculateJacobian(
            GetGeometry(),
            DN_De,
            3, //WorkingSpaceDimension,
            2, //LocalSpaceDimension,
            metric.J
        );
       
        metric.g1[0] = metric.J(0, 0);
        metric.g2[0] = metric.J(0, 1);
        metric.g1[1] = metric.J(1, 0);
        metric.g2[1] = metric.J(1, 1);
        metric.g1[2] = metric.J(2, 0);
        metric.g2[2] = metric.J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(metric.g3, metric.g1, metric.g2);
        //differential area dA
        metric.dA = norm_2(metric.g3);
        //normal vector _n
        Vector n = metric.g3 / metric.dA;

        //GetcovariantMetric
        metric.gab[0] = pow(metric.g1[0], 2) + pow(metric.g1[1], 2) + pow(metric.g1[2], 2);
        metric.gab[1] = pow(metric.g2[0], 2) + pow(metric.g2[1], 2) + pow(metric.g2[2], 2);
        metric.gab[2] = metric.g1[0] * metric.g2[0] + metric.g1[1] * metric.g2[1] + metric.g1[2] * metric.g2[2];

       
        //contravariant metric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (metric.gab[0] * metric.gab[1] - metric.gab[2] * metric.gab[2]);
        metric.gab_con[0] = invdetGab*metric.gab[1];
        metric.gab_con[2] = -invdetGab*metric.gab[2];
        metric.gab_con[1] = invdetGab*metric.gab[0];

        array_1d<double, 3> g_con_1 = metric.g1*metric.gab_con[0] + metric.g2*metric.gab_con[2];
        array_1d<double, 3> g_con_2 = metric.g1*metric.gab_con[2] + metric.g2*metric.gab_con[1];


        //local cartesian coordinates
        double lg1 = norm_2(metric.g1);
        array_1d<double, 3> e1 = metric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        Matrix mG = ZeroMatrix(2, 2);
        mG(0, 0) = inner_prod(e1, g_con_1);
        mG(0, 1) = inner_prod(e1, g_con_2);
        mG(1, 0) = inner_prod(e2, g_con_1);
        mG(1, 1) = inner_prod(e2, g_con_2);

        metric.Q = ZeroMatrix(3, 3);
        metric.Q(0, 0) = pow(mG(0, 0), 2);
        metric.Q(0, 1) = pow(mG(0, 1), 2);
        metric.Q(0, 2) = 2.00*mG(0, 0)*mG(0, 1);

        metric.Q(1, 0) = pow(mG(1, 0), 2);
        metric.Q(1, 1) = pow(mG(1, 1), 2);
        metric.Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

        metric.Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
        metric.Q(2, 1) = 2.00 * mG(0, 1)*mG(1, 1);
        metric.Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));
    }


    void IgaMembraneElement::CalculateStrain(
        Vector& StrainVector,
        Vector& gab,
        Vector& gab0)

    {
        KRATOS_TRY

        StrainVector[0] = 0.5 * (gab[0] - gab0[0]);
        StrainVector[1] = 0.5 * (gab[1] - gab0[1]);
        StrainVector[2] = 0.5 * (gab[2] - gab0[2]); 

        KRATOS_CATCH("")
    }


    void IgaMembraneElement::CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        Vector strain_vector = ZeroVector(3);

        CalculateStrain(strain_vector, rActualMetric.gab, mInitialMetric.gab);
        rThisConstitutiveVariablesMembrane.E = prod(mInitialMetric.Q, strain_vector); 


        //Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.D); //this is an ouput parameter
      
        //rValues.CheckAllParameters();
        mConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure);

        double thickness = this->GetProperties().GetValue(THICKNESS);
         rThisConstitutiveVariablesMembrane.D *= thickness;
        
        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.S = prod(
           trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E); 
             
    }


        void IgaMembraneElement::CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& metric) 
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        if (rB.size1() != mat_size || rB.size2() != mat_size)
            rB.resize(mat_size, mat_size);
        rB = ZeroMatrix(3, mat_size);

        for (int r = 0; r<mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;
            Vector dE_curvilinear = ZeroVector(3);
            // strain
            dE_curvilinear[0] = DN_De(kr, 0)*metric.g1(dirr);
            dE_curvilinear[1] = DN_De(kr, 1)*metric.g2(dirr);
            dE_curvilinear[2] = 0.5*(DN_De(kr, 0)*metric.g2(dirr) + metric.g1(dirr)*DN_De(kr, 1));

            rB(0, r) = mInitialMetric.Q(0, 0)*dE_curvilinear[0] + mInitialMetric.Q(0, 1)*dE_curvilinear[1] + mInitialMetric.Q(0, 2)*dE_curvilinear[2];  //geändert Q->T HIER
            rB(1, r) = mInitialMetric.Q(1, 0)*dE_curvilinear[0] + mInitialMetric.Q(1, 1)*dE_curvilinear[1] + mInitialMetric.Q(1, 2)*dE_curvilinear[2];  //geändert Q->T HIER
            rB(2, r) = mInitialMetric.Q(2, 0)*dE_curvilinear[0] + mInitialMetric.Q(2, 1)*dE_curvilinear[1] + mInitialMetric.Q(2, 2)*dE_curvilinear[2];  //geändert Q->T HIERS
        }
    }

void IgaMembraneElement::CalculateSecondVariationStrain(
        SecondVariations& rSecondVariationsStrain,
        const MetricVariables& rMetric)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        for (int r = 0; r<mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;

            for (int s = 0; s <= r; s++)
            {
                // local node number ks and dof direction dirs
                int ks = s / 3;
                int dirs = s % 3;

                // strain
                Vector ddE_curvilinear = ZeroVector(3);
                if (dirr == dirs)
                {
                    ddE_curvilinear[0] = DN_De(kr, 0)*DN_De(ks, 0);
                    ddE_curvilinear[1] = DN_De(kr, 1)*DN_De(ks, 1);
                    ddE_curvilinear[2] = 0.5*(DN_De(kr, 0)*DN_De(ks, 1) + DN_De(kr, 1)*DN_De(ks, 0));
                }

                rSecondVariationsStrain.B11(r, s) = mInitialMetric.Q(0, 0)*ddE_curvilinear[0] + mInitialMetric.Q(0, 1)*ddE_curvilinear[1] + mInitialMetric.Q(0, 2)*ddE_curvilinear[2]; //HIER
                rSecondVariationsStrain.B22(r, s) = mInitialMetric.Q(1, 0)*ddE_curvilinear[0] + mInitialMetric.Q(1, 1)*ddE_curvilinear[1] + mInitialMetric.Q(1, 2)*ddE_curvilinear[2]; //HIER
                rSecondVariationsStrain.B12(r, s) = mInitialMetric.Q(2, 0)*ddE_curvilinear[0] + mInitialMetric.Q(2, 1)*ddE_curvilinear[1] + mInitialMetric.Q(2, 2)*ddE_curvilinear[2]; //HIER
            }
        }
    }
    
    void IgaMembraneElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight
    )
    {
        KRATOS_TRY
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
        KRATOS_CATCH("")    
    }


    void IgaMembraneElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)

    {
        KRATOS_TRY
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        for (int n = 0; n < mat_size; n++)
        {
            for (unsigned int m = 0; m <= n; m++)
            {
            
                double nm = (SD[0] * SecondVariationsStrain.B11(n, m)
                    + SD[1] * SecondVariationsStrain.B22(n, m)
                    + SD[2] * SecondVariationsStrain.B12(n, m))*rIntegrationWeight;
                //KRATOS_ERROR_IF(nm > 1e-5) << "Problem" << std::endl;
                rLeftHandSideMatrix(n, m) += nm;
                if(n!=m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }
        KRATOS_CATCH("")
    }

    //Transformationsmatrix für Prestress
    void IgaMembraneElement::CalculateTransformationmatrixPrestress(
        const MetricVariables& metric,
        PrestresstransVariables& Prestresstrans
         )
    {
    //define base vector in reference plane ATTENSION: in some cases the vector must be modified (e.g. catenoid t3_unten=[0, 0, 1])
    Vector t3_unten = ZeroVector(3);
    t3_unten[0] = 0;
    t3_unten[1] = 1;
    t3_unten[2] = 0;

    //calculate first base vector in tangential plane
    Vector t1_z = ZeroVector(3);
    MathUtils<double>::CrossProduct(t1_z, t3_unten, metric.g3);

    //calculate second base vector in tangential plane
    Vector t2 = ZeroVector(3);
    MathUtils<double>::CrossProduct(t2, metric.g3, t1_z); 

    //Normalization of the vectors
    Vector t1_n = t1_z/norm_2(t1_z); 
    Vector t2_n = t2/norm_2(t2); 
    Vector t3_n = metric.g3/norm_2(metric.g3); 
    
    array_1d<double, 3> g_con_1 = metric.g1*metric.gab_con[0] + metric.g2*metric.gab_con[2];
    array_1d<double, 3> g_con_2 = metric.g1*metric.gab_con[2] + metric.g2*metric.gab_con[1]; 
 
    //local cartesian coordinates oriented along the 1st base vector in the ref. config.
        double lg1 = norm_2(metric.g1); 
        array_1d<double, 3> e1 = metric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2); 
        array_1d<double, 3> e2 = g_con_2 / lg_con2;
               

    //Transformation matrix from the projected basis T to the local cartesian basis
   
    double eG11 = inner_prod(e1,t1_n);
    double eG12 = inner_prod(e1,t2_n);
    double eG21 = inner_prod(e2,t1_n);
    double eG22 = inner_prod(e2,t2_n);

    Prestresstrans.Tpre = ZeroMatrix(3, 3);
    Prestresstrans.Tpre(0,0) = eG11*eG11;
    Prestresstrans.Tpre(0,1) = eG12*eG12;
    Prestresstrans.Tpre(0,2) = 2.0*eG11*eG12;

    Prestresstrans.Tpre(1,0) = eG21*eG21;
    Prestresstrans.Tpre(1,1) = eG22*eG22;
    Prestresstrans.Tpre(1,2) = 2.0*eG21*eG22;

    Prestresstrans.Tpre(2,0) = eG11*eG21;
    Prestresstrans.Tpre(2,1) = eG12*eG22;
    Prestresstrans.Tpre(2,2) = eG11*eG22+eG12*eG21;
    }

//************************************************************************
//************************************************************************
//************************************************************************
//POSTPROCESSING
void IgaMembraneElement::Calculate(
        const Variable<double>& rVariable,
        double& rOutput, //std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )// override;
    {
        if (rVariable == PRINCIPAL_STRESS_1)
        {
            Vector stresses = ZeroVector(3);
            CalculateStresses(stresses, rCurrentProcessInfo); 

            // Principal normal forces
            double principal_stresses = 0.5*(stresses[0] + stresses[1] + sqrt(pow(stresses[0] - stresses[1], 2)+4*pow(stresses[2], 2))) ;
           
            rOutput = principal_stresses;
        }   
        if(rVariable == PRINCIPAL_STRESS_2)
        {
            Vector stresses = ZeroVector(3);
            CalculateStresses(stresses, rCurrentProcessInfo); 
            double principal_stresses = 0.5*(stresses[0] + stresses[1] - sqrt(pow(stresses[0] - stresses[1], 2)+4*pow(stresses[2], 2))) ;
                   
            rOutput = principal_stresses;
        }  
        else
        {
            //SurfaceBaseDiscreteElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
    
    void IgaMembraneElement::CalculateStresses(
        Vector& rStresses,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rStresses.size() != 3)
            rStresses.resize(3);
        noalias(rStresses) = ZeroVector(3); //resetting LHS
        

        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables(3);
        CalculateConstitutiveVariables(actual_metric, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

        double detF = actual_metric.dA / mInitialMetric.dA;

        Vector prestress_tensor = ZeroVector(3);
        CalculatePresstressTensor(prestress_tensor, actual_metric);
        constitutive_variables.S += prestress_tensor;

        // PK2 normal force in local cartesian e1, e2
        rStresses[0] = constitutive_variables.S[0];
        rStresses[1] = constitutive_variables.S[1];
        rStresses[2] = constitutive_variables.S[2];
    }

    void IgaMembraneElement::CalculatePresstressTensor(
        Vector& rPrestressTensor,
        MetricVariables& rMetric
    )
    {
        rPrestressTensor.resize(3);

        double thickness = this->GetProperties().GetValue(THICKNESS);

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
       
        Vector S_prestress_nichttransformiert = this->GetProperties().GetValue(PRESTRESS);
        PrestresstransVariables prestresstrans_variables(3);
        CalculateTransformationmatrixPrestress(
            actual_metric,
            prestresstrans_variables 
        );
        Vector S_prestress = prod(prestresstrans_variables.Tpre, S_prestress_nichttransformiert); 

        rPrestressTensor = S_prestress * thickness;
    }

//**************************************************************************
}//namespace Kratos