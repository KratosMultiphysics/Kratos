//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_INCOMPRESSIBLE_STRESSES_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_INCOMPRESSIBLE_STRESSES_POTENTIAL_FLOW_ELEMENT_H_INCLUDED

// #define SYMMETRIC_CONSTRAINT_APPLICATION

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "utilities/enrichment_utilities.h"
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

template< int Dim, int NumNodes >
class IncompressibleStressesPotentialFlowElement : public Element
{
public:

    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double,TNumNodes> phis, distances;
        double rho;
        double vol;

        bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
    };



    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressibleStressesPotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressibleStressesPotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    IncompressibleStressesPotentialFlowElement(IndexType NewId = 0) {};

    /**
     * Constructor using an array of nodes
     */
    IncompressibleStressesPotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes):Element(NewId, ThisNodes) {};

    /**
     * Constructor using Geometry
     */
    IncompressibleStressesPotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry):Element(NewId, pGeometry) {};

    /**
     * Constructor using Properties
     */
    IncompressibleStressesPotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Element(NewId, pGeometry, pProperties) {};

    /**
     * Copy Constructor
     */
    IncompressibleStressesPotentialFlowElement(IncompressibleStressesPotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~IncompressibleStressesPotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressibleStressesPotentialFlowElement & operator=(IncompressibleStressesPotentialFlowElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleStressesPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleStressesPotentialFlowElement(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleStressesPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override
    {
        if(this->IsNot(MARKER)) //normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            for (unsigned int i = 0; i < NumNodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(POSITIVE_POTENTIAL).EquationId();

        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = GetGeometry()[i].GetDof(POSITIVE_POTENTIAL).EquationId();
                else
                    rResult[i] = GetGeometry()[i].GetDof(NEGATIVE_POTENTIAL,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(POSITIVE_POTENTIAL).EquationId();
                else
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(NEGATIVE_POTENTIAL,0).EquationId();
            }
        }


    }


    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override
    {
        if(this->IsNot(MARKER)) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            for (unsigned int i = 0; i < NumNodes; i++)
                rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_POTENTIAL);
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_POTENTIAL);
                else
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(NEGATIVE_POTENTIAL);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(POSITIVE_POTENTIAL);
                else
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(NEGATIVE_POTENTIAL);
            }
        }
    }


    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {
        ElementalData<NumNodes,Dim> data;
        array_1d<double,NumNodes> elemental_distance;
        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        //gather nodal data
        for(unsigned int i=0; i<NumNodes; i++)
            data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);

        
        //TEST:
        bool kutta_element = false;
        for(unsigned int i=0; i<NumNodes; ++i)
            if(GetGeometry()[i].Is(STRUCTURE))
            {
                kutta_element = true;
                break;
            }

        if(this->IsNot(MARKER))//normal element (non-wake) - eventually an embedded
        {
            if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
                rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
            if (rRightHandSideVector.size() != NumNodes)
                rRightHandSideVector.resize(NumNodes, false);
            rLeftHandSideMatrix.clear();

            if (this->Is(BOUNDARY)){
                for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
                    elemental_distance[i_node] = GetGeometry()[i_node].GetSolutionStepValue(LEVEL_SET_DISTANCE);

                const Vector& r_elemental_distances=elemental_distance;
                Triangle2D3ModifiedShapeFunctions triangle_shape_functions(pGetGeometry(), r_elemental_distances);
                Matrix positive_side_sh_func;
                ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
                Vector positive_side_weights;
                triangle_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
                    positive_side_sh_func,
                    positive_side_sh_func_gradients,
                    positive_side_weights,
                    GeometryData::GI_GAUSS_2);

                    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
                        MatrixType aux_matrix;
                        bounded_matrix<double,NumNodes,Dim> DN_DX;
                        DN_DX=positive_side_sh_func_gradients(i_gauss);
                        
                        //reading properties and conditions
                        aux_matrix=prod(DN_DX,trans(DN_DX))*positive_side_weights(i_gauss);  // Bt D B

                        noalias(rLeftHandSideMatrix) += aux_matrix;                       
                    }
                noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
            }
            else //if (this->Is(FLUID) || this->IsNotDefined(FLUID)){
                ComputeLHSGaussPointContribution(data.vol,rLeftHandSideMatrix,data);          
            noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
        }
        else //it is a wake element
        {   
            if(this->Is(BOUNDARY)) std::cout<<"Element with both embedded Wake & LevelSet functions:"<<this->Id()<<std::endl;
            GetWakeDistances(data.distances);
            
            //note that the lhs and rhs have double the size!!
            if(rLeftHandSideMatrix.size1() != 2*NumNodes || rLeftHandSideMatrix.size2() != 2*NumNodes)
                rLeftHandSideMatrix.resize(2*NumNodes,2*NumNodes,false);
            if(rRightHandSideVector.size() != 2*NumNodes)
                rRightHandSideVector.resize(2*NumNodes,false);
            rLeftHandSideMatrix.clear();
            
            //subdivide the element
            constexpr unsigned int nvolumes = 3*(Dim-1);
            bounded_matrix<double,NumNodes, Dim > Points;
            array_1d<double,nvolumes> Volumes;
            bounded_matrix<double, nvolumes, NumNodes > GPShapeFunctionValues;
            array_1d<double,nvolumes> PartitionsSign;
            std::vector<Matrix> GradientsValue(nvolumes);
            bounded_matrix<double,nvolumes, 2> NEnriched;
            
            for(unsigned int i=0; i<GradientsValue.size(); ++i)
                GradientsValue[i].resize(2,Dim,false);
           
            
            
            for(unsigned int i = 0; i<NumNodes; ++i)
            {
                const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
                for(unsigned int k = 0; k<Dim; ++k)
                {
                    Points(i, k) = coords[k];
                }
            }
            
            const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(Points,
                                                                                            data.DN_DX,
                                                                                            data.distances,
                                                                                            Volumes, 
                                                                                            GPShapeFunctionValues, 
                                                                                            PartitionsSign, 
                                                                                            GradientsValue, 
                                                                                            NEnriched);
            //compute the lhs and rhs that would correspond to it not being divided
            Matrix lhs_positive = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_negative = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_uncondensed=ZeroMatrix(2*NumNodes+2,2*NumNodes+2);

            Matrix lhs_plus_sigma = ZeroMatrix(NumNodes,Dim);
            Matrix lhs_sigma_plus = ZeroMatrix(Dim,NumNodes);
            Matrix lhs_sigma_sigma = ZeroMatrix(Dim,Dim);
            Matrix lhs_sigma_sigma_no_inv = ZeroMatrix(Dim,Dim);
            Matrix sigma_shape_func = ZeroMatrix(Dim,Dim);
            Matrix lhs_minus_sigma = ZeroMatrix(NumNodes,Dim);
            Matrix lhs_sigma_minus = ZeroMatrix(Dim,NumNodes);

            Matrix lhs_plus_plus = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_plus_minus = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_minus_plus = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_minus_minus = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_plus_rhs = ZeroMatrix(NumNodes,Dim);
            Matrix lhs_minus_rhs = ZeroMatrix(NumNodes,Dim);
            Matrix forcing_sigma_plus =ZeroMatrix(NumNodes,Dim);
            Matrix forcing_sigma_minus = ZeroMatrix(NumNodes,Dim); 
            Matrix forcing_sigma =ZeroMatrix(NumNodes,Dim);

            Matrix K_uu = ZeroMatrix(NumNodes*2,NumNodes*2);

            for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
                elemental_distance[i_node] = GetGeometry()[i_node].GetSolutionStepValue(WAKE_DISTANCE);               
                
            const Vector& r_elemental_distances=elemental_distance;
                
            Triangle2D3ModifiedShapeFunctions triangle_shape_functions(this->pGetGeometry(), r_elemental_distances);
            // Matrix positive_side_interface_sh_func;
            // ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_interface_gradients;
            // Vector positive_side_interface_weights;
            // triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            //     positive_side_interface_sh_func,
            //     positive_side_sh_func_interface_gradients,
            //     positive_side_interface_weights,
            //     GeometryData::GI_GAUSS_1);

            std::vector<Vector> cut_unit_normal;
            triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(cut_unit_normal,GeometryData::GI_GAUSS_1);

            Vector sigma(Dim);
            Vector n(Dim);
            Vector shape_func(NumNodes);
            for (unsigned int i = 0; i<NumNodes; i++)
                shape_func(i) = data.N(i);
            n(0)=cut_unit_normal[0][0];
            n(1)=cut_unit_normal[0][1];
            sigma_shape_func(0,0)=1.0;
            sigma_shape_func(1,1)=1.0;
            double norm_n=sqrt(std::abs(inner_prod(n,n)));
            n(0)=n(0)/norm_n;
            n(1)=n(1)/norm_n;
            // std::cout<<"Cut length: "<< sqrt(std::abs(inner_prod(n,n))); 
            Vector normal_gradient=prod(n,sigma_shape_func);
            
            for (unsigned int i = 0; i<NumNodes; ++i){
                for (unsigned int j= 0; j<Dim; ++j){
                    forcing_sigma(i,j)=0;//-data.N(i)*normal_gradient(j);
                }
            }
            
            double n_parameter=rCurrentProcessInfo[INITIAL_PENALTY];

            lhs_sigma_sigma=n_parameter*sigma_shape_func/data.vol;
            lhs_sigma_sigma_no_inv=-1/n_parameter*sigma_shape_func*data.vol;
            noalias(lhs_sigma_plus) = 1.0/n_parameter*data.vol*prod(sigma_shape_func,trans(data.DN_DX));
            noalias(lhs_sigma_minus) = -1.0/n_parameter*data.vol*prod(sigma_shape_func,trans(data.DN_DX));

            for(unsigned int i=0; i<nsubdivisions; ++i)
            {
                if(PartitionsSign[i] > 0){
                    ComputeLHSGaussPointContribution(Volumes[i],lhs_positive,data); //K++                    
                    
                    noalias(lhs_plus_sigma) += 1.0/n_parameter*Volumes[i]*prod(data.DN_DX,sigma_shape_func);
                    // noalias(lhs_sigma_plus) += 1.0/n_parameter*Volumes[i]*prod(sigma_shape_func,trans(data.DN_DX));
                    
                    noalias(forcing_sigma_plus) += -Volumes[i]*forcing_sigma;
                }
                else{
                    ComputeLHSGaussPointContribution(Volumes[i],lhs_negative,data); //K--

                    noalias(lhs_minus_sigma) += 1.0/n_parameter*Volumes[i]*prod(data.DN_DX,sigma_shape_func);
                    // noalias(lhs_sigma_minus) += -1.0/n_parameter*Volumes[i]*prod(sigma_shape_func,trans(data.DN_DX));
                    
                    noalias(forcing_sigma_minus) += Volumes[i]*forcing_sigma;
                }
            }
            // KRATOS_WATCH(forcing_sigma)
            // KRATOS_WATCH(lhs_sigma_minus)
            // KRATOS_WATCH(lhs_sigma_plus)
            // KRATOS_WATCH(Volumes)
            // KRATOS_WATCH(data.DN_DX)



            noalias(lhs_plus_plus) = prod(lhs_plus_sigma+forcing_sigma,Matrix(prod(lhs_sigma_sigma,lhs_sigma_plus)));
            noalias(lhs_plus_minus) = prod(lhs_plus_sigma+forcing_sigma,Matrix(prod(lhs_sigma_sigma,lhs_sigma_minus)));   
            noalias(lhs_minus_plus) = prod(lhs_minus_sigma-forcing_sigma,Matrix(prod(lhs_sigma_sigma,lhs_sigma_plus)));
            noalias(lhs_minus_minus) = prod(lhs_minus_sigma-forcing_sigma,Matrix(prod(lhs_sigma_sigma,lhs_sigma_minus)));
            noalias(lhs_plus_rhs) = prod(lhs_plus_sigma+forcing_sigma,lhs_sigma_sigma);
            noalias(lhs_minus_rhs) = prod(lhs_minus_sigma-forcing_sigma,lhs_sigma_sigma);

            if(kutta_element)
            {
                for(unsigned int i=0; i<NumNodes; ++i)
                {
                    for(unsigned int j=0; j<NumNodes; ++j)
                    {
                        rLeftHandSideMatrix(i,j)                   =  lhs_positive(i,j); 
                        rLeftHandSideMatrix(i,j+NumNodes)          =  0.0;
                        K_uu(i,j)                                  =  lhs_positive(i,j);
                        K_uu(i,j+NumNodes)                         =  0.0;  
                        
                        rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  lhs_negative(i,j);
                        rLeftHandSideMatrix(i+NumNodes,j)          =  0.0;
                        K_uu(i+NumNodes,j+NumNodes)                =  lhs_negative(i,j);
                        K_uu(i+NumNodes,j)                         =  0.0;
                    }
                }
            }
            else
            {
                for(unsigned int i=0; i<NumNodes; ++i)
                {
                    for(unsigned int j=0; j<NumNodes; ++j)
                    {
                        rLeftHandSideMatrix(i,j)                   =  (1.0-1.0/n_parameter)*lhs_positive(i,j)+lhs_plus_plus(i,j); 
                        rLeftHandSideMatrix(i,j+NumNodes)          =  0.0;
                        K_uu(i,j)                                  =  lhs_positive(i,j);                  
                        K_uu(i,j+NumNodes)                         =  0.0;  
                        lhs_uncondensed(i,j)                       =  lhs_positive(i,j);                  
                        
                        rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  (1.0-1.0/n_parameter)*lhs_negative(i,j)-lhs_minus_minus(i,j);
                        rLeftHandSideMatrix(i+NumNodes,j)          =  0.0;
                        K_uu(i+NumNodes,j+NumNodes)                =  lhs_negative(i,j);
                        K_uu(i+NumNodes,j)                         =  0.0;
                        lhs_uncondensed(i+NumNodes,j+NumNodes)     =  lhs_negative(i,j);
                    }
                    lhs_uncondensed(i,2*NumNodes)                  =  lhs_plus_sigma(i,0)+forcing_sigma(i,0);
                    lhs_uncondensed(i,2*NumNodes+1)                =  lhs_plus_sigma(i,1)+forcing_sigma(i,1);
                    lhs_uncondensed(2*NumNodes,i)                  =  lhs_sigma_plus(0,i);
                    lhs_uncondensed(2*NumNodes+1,i)                =  lhs_sigma_plus(1,i);

                    lhs_uncondensed(i+NumNodes,2*NumNodes)         =  lhs_minus_sigma(i,0)+forcing_sigma(i,0);
                    lhs_uncondensed(i+NumNodes,2*NumNodes+1)       =  lhs_minus_sigma(i,1)+forcing_sigma(i,1);
                    lhs_uncondensed(2*NumNodes,i+NumNodes)         =  lhs_sigma_minus(0,i);
                    lhs_uncondensed(2*NumNodes+1,i+NumNodes)       =  lhs_sigma_minus(1,i);

                    lhs_uncondensed(2*NumNodes,2*NumNodes)         =  lhs_sigma_sigma_no_inv(0,0);
                    lhs_uncondensed(2*NumNodes,2*NumNodes+1)       =  lhs_sigma_sigma_no_inv(0,1);
                    lhs_uncondensed(2*NumNodes+1,2*NumNodes)       =  lhs_sigma_sigma_no_inv(1,0);
                    lhs_uncondensed(2*NumNodes+1,2*NumNodes+1)     =  lhs_sigma_sigma_no_inv(1,1);
                    
                }
                
            
                //side1  -assign constraint only on the NEGATIVE_POTENTIAL dofs
                for(unsigned int i=0; i<NumNodes; ++i)
                {
                    if(data.distances[i]<0)
                    {
                        for(unsigned int j=0; j<NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i,j+NumNodes) = lhs_plus_minus(i,j);
                        }
                    }
                }
                
                //side2 -assign constraint only on the NEGATIVE_POTENTIAL dofs
                for(unsigned int i=0; i<NumNodes; ++i)
                {                            
                    if(data.distances[i]>0)
                    {   
                        for(unsigned int j=0; j<NumNodes; ++j)
                            {
                                rLeftHandSideMatrix(i+NumNodes,j) = lhs_minus_plus(i,j); 
                            }
                    }
                }
            }
            Vector split_element_values(NumNodes*2);
            Vector split_element_values_uncondensed=ZeroVector(NumNodes*2+2);
            Vector split_plus(NumNodes);
            Vector split_minus(NumNodes);
            Vector split_prev_plus(NumNodes);
            Vector split_prev_minus(NumNodes);
            Vector residual_sigma(Dim);
            Vector sigma_aux(Dim);
            Vector residual_sigma_prev(Dim);
            Vector rhs_sigma(NumNodes*2);
            Vector rhs_sigma_plus(NumNodes);
            Vector rhs_sigma_minus(NumNodes);
            Vector rhs_uncondensed(2*NumNodes+2);

            GetValuesOnSplitElement(split_element_values, data.distances);
            
            for (unsigned int i = 0; i<NumNodes; ++i){
                split_plus(i) = split_element_values(i);
                split_minus(i)= split_element_values(i+NumNodes);
                split_element_values_uncondensed(i) = split_element_values(i);
                split_element_values_uncondensed(i+NumNodes)= split_element_values(i+NumNodes);
            }   
            split_element_values_uncondensed(2*NumNodes)=0;
            split_element_values_uncondensed(2*NumNodes+1)=0;
      
            // residual_sigma_prev(0)=this->GetValue(Y1);
            // residual_sigma_prev(1)=this->GetValue(Y2);
            // sigma_aux = -prod(lhs_sigma_plus,split_plus) - prod(lhs_sigma_minus,split_minus) - residual_sigma_prev;
            // sigma = prod(lhs_sigma_sigma,sigma_aux);
            // std::cout<<sigma<< std::endl;
            KRATOS_WATCH(split_element_values);
            KRATOS_WATCH(lhs_negative);
            KRATOS_WATCH(lhs_minus_minus);
            residual_sigma = -prod(lhs_sigma_plus,split_plus) - prod(lhs_sigma_minus,split_minus);//-prod(lhs_sigma_sigma_no_inv,trans(sigma));
            sigma = prod(lhs_sigma_sigma,-residual_sigma);
            KRATOS_WATCH(sigma);
            // this->GetValue(Y1) = residual_sigma(0);
            // this->GetValue(Y2) = residual_sigma(1);

            KRATOS_WATCH(residual_sigma);
            rhs_sigma_plus=prod(lhs_plus_rhs,residual_sigma);//+prod(forcing_sigma_plus,trans(sigma));
            rhs_sigma_minus=prod(lhs_minus_rhs,residual_sigma);//+prod(forcing_sigma_minus,trans(sigma));

            for (unsigned int i = 0; i < NumNodes; ++i){
                rhs_sigma(i) = rhs_sigma_plus(i);                
                rhs_sigma(i+NumNodes) = rhs_sigma_minus(i);
            }
            KRATOS_WATCH(rLeftHandSideMatrix);
            KRATOS_WATCH(rhs_sigma);
            noalias(rRightHandSideVector) = -prod(K_uu,split_element_values)+rhs_sigma;
            noalias(rhs_uncondensed) = -prod(lhs_uncondensed,split_element_values_uncondensed);
            KRATOS_WATCH(rhs_uncondensed);
            KRATOS_WATCH(lhs_uncondensed);
        }
        
    }


    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        //TODO: improve speed
        Matrix tmp;
        CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (this->Is(MARKER) && active == true)
            CheckWakeCondition();      
    }

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {

        KRATOS_TRY

        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "IncompressibleStressesPotentialFlowElement found with Id 0 or negative", "")
        }

        if (this->GetGeometry().Area() <= 0)
        {
            std::cout << "error on IncompressibleStressesPotentialFlowElement -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0", "")
        }

        for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
        {
            if (this->GetGeometry()[i].SolutionStepsDataHas(POSITIVE_POTENTIAL) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing variable POSITIVE_POTENTIAL on node ", this->GetGeometry()[i].Id())
        }

        return 0;

        KRATOS_CATCH("");
    }

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        if(rValues.size() != 1) rValues.resize(1);

        if (rVariable == PRESSURE)
        {
            double p = 0.0;
            p = ComputePressure(rCurrentProcessInfo);
            rValues[0] = p;
        }
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        if(rValues.size() != 1) rValues.resize(1);

        if (rVariable == VELOCITY)
        {
            array_1d<double,3> v(3,0.0);
            array_1d<double,Dim> vaux;
            ComputeVelocity(vaux);
            for(unsigned int k=0; k<Dim; k++) v[k] = vaux[k];
            rValues[0] = v;
        }
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IncompressibleStressesPotentialFlowElement #" << Id();
        return buffer.str();
    }

/// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressibleStressesPotentialFlowElement #" << Id();
    }

/// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }



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
    void GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }

    void ComputeLHSGaussPointContribution(
        const double weight,
        Matrix &lhs,
        const ElementalData<NumNodes, Dim> &data)
    {
        noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
    }

    void ComputeRHSGaussPointContribution(
        const double weight,
        Vector &rhs,
        const ElementalData<NumNodes, Dim> &data)
    {
        array_1d<double, Dim> grad = prod(trans(data.DN_DX), data.phis);
        noalias(rhs) -= weight * prod(data.DN_DX, grad);
    }

    void GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
    {

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] < 0)
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
            else
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
        }
    }

    void ComputeVelocity(array_1d<double,Dim>& velocity)
    {
        velocity.clear();

        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (this->IsNot(MARKER) && active == true)
            ComputeVelocityNormalElement(velocity);
        else if (active == true && this->Is(MARKER))
            ComputeVelocityUpperWakeElement(velocity);
    }

    void ComputeVelocityNormalElement(array_1d<double,Dim>& velocity)
    {
        ElementalData<NumNodes, Dim> data;
        for (unsigned int i = 0; i < NumNodes; i++)
            data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);

        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);    

        noalias(velocity) = -prod(trans(data.DN_DX), data.phis);
    }

    void ComputeVelocityUpperWakeElement(array_1d<double,Dim>& velocity)
    {

        array_1d<double, NumNodes> distances;
        ElementalData<NumNodes, Dim> data;
        GetWakeDistances(distances);

        //taking only positive part
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0)
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
            else
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
        }
        
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        noalias(velocity) = -prod(trans(data.DN_DX), data.phis);
    }

    void ComputeVelocityLowerWakeElement(array_1d<double,Dim>& velocity)
    {
        ElementalData<NumNodes, Dim> data;

        array_1d<double, NumNodes> distances;
        GetWakeDistances(distances);

        //taking only negative part
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] < 0)
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
            else
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
        }
            
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        noalias(velocity) = -prod(trans(data.DN_DX), data.phis);
    }

    void CheckWakeCondition()
    {
        array_1d<double, Dim> upper_wake_velocity;
        ComputeVelocityUpperWakeElement(upper_wake_velocity);
        const double vupnorm = inner_prod(upper_wake_velocity, upper_wake_velocity);

        array_1d<double, Dim> lower_wake_velocity;
        ComputeVelocityLowerWakeElement(lower_wake_velocity);
        const double vlownorm = inner_prod(lower_wake_velocity, lower_wake_velocity);

        if (std::abs(vupnorm - vlownorm) > 0.1)
            std::cout << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << this->Id() <<"    " <<std::abs(vupnorm - vlownorm)<<std::endl;
    }

    double ComputePressure(const ProcessInfo& rCurrentProcessInfo)
    {
        double pressure = 0.0;

        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (active && !this->Is(MARKER))
            pressure = ComputePressureNormalElement(rCurrentProcessInfo);
        else if (active == true && this->Is(MARKER))
            pressure = ComputePressureWakeElement(rCurrentProcessInfo);

        return pressure;
    }

    double ComputePressureNormalElement(const ProcessInfo& rCurrentProcessInfo)
    {
        double pressure = 0.0;
        const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

        KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
            << "Error on element -> " << this->Id() << "\n"
            << "vinfinity_norm2 must be larger than zero." << std::endl;

        array_1d<double, Dim> v;
        ComputeVelocityNormalElement(v);

        pressure = (vinfinity_norm2 - inner_prod(v, v)) / vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

        return pressure;
    }

    double ComputePressureWakeElement(const ProcessInfo& rCurrentProcessInfo)
    {
        double pressure = 0.0;
        const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

        KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
            << "Error on element -> " << this->Id() << "\n"
            << "vinfinity_norm2 must be larger than zero." << std::endl;

        array_1d<double, Dim> v;
        ComputeVelocityLowerWakeElement(v);

        pressure = (vinfinity_norm2 - inner_prod(v, v)) / vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

        return pressure;
    }

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
    }

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

}; // Class IncompressibleStressesPotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
