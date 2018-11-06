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

#if !defined(KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED

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
class CompressiblePotentialFlowElement : public Element
{
public:

    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double,TNumNodes> phis, distances;
        double density;
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
    /// Pointer definition of CompressiblePotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(CompressiblePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    CompressiblePotentialFlowElement(IndexType NewId = 0) {};

    /**
     * Constructor using an array of nodes
     */
    CompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes):Element(NewId, ThisNodes) {};

    /**
     * Constructor using Geometry
     */
    CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry):Element(NewId, pGeometry) {};

    /**
     * Constructor using Properties
     */
    CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Element(NewId, pGeometry, pProperties) {};

    /**
     * Copy Constructor
     */
    CompressiblePotentialFlowElement(CompressiblePotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~CompressiblePotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CompressiblePotentialFlowElement & operator=(CompressiblePotentialFlowElement const& rOther)
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
        return Element::Pointer(new CompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
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
        return Element::Pointer(new CompressiblePotentialFlowElement(NewId, pGeom, pProperties));
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
        return Element::Pointer(new CompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
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
        array_1d<double,NumNodes> elemental_distance;

        //Matrix to compute the residual
        MatrixType rLaplacianMatrix;

        ElementalData<NumNodes,Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
        array_1d<double,Dim> velocity;
        double density = GetProperties().GetValue(DENSITY);
        double derivative = 0;
        Vector DNV = ZeroVector(NumNodes);

        //gather nodal data
        for(unsigned int i=0; i<NumNodes; i++)
            data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);

        
        //TEST:
        bool kutta_element = false;
        if (this->Is(INTERFACE))
            kutta_element = true;

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
            if(rLaplacianMatrix.size1() != NumNodes || rLaplacianMatrix.size2() != NumNodes)
                rLaplacianMatrix.resize(NumNodes,NumNodes,false);
            rLaplacianMatrix.clear();

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
                    GeometryData::GI_GAUSS_3);
                this->ComputeDensity(density,derivative,velocity,rCurrentProcessInfo);
                for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
                    bounded_matrix<double,NumNodes,Dim> DN_DX;
                    DN_DX=positive_side_sh_func_gradients(i_gauss);
                    DNV=prod(DN_DX, velocity);

                    noalias(rLeftHandSideMatrix) += density*(prod(DN_DX,trans(DN_DX)))*positive_side_weights(i_gauss); 
                    noalias(rLeftHandSideMatrix) += derivative*outer_prod(DNV, trans(DNV))*positive_side_weights(i_gauss);     
                    noalias(rLaplacianMatrix) += density*prod(DN_DX, trans(DN_DX))*positive_side_weights(i_gauss);
                }

                noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
            }
            else {
                if (this->Is(FLUID) || this->IsNotDefined(FLUID)){
                        this->ComputeDensity(density,derivative,velocity,rCurrentProcessInfo);
                        DNV = prod(data.DN_DX, velocity);
                        //std::cout << "density normal = " << density  << std::endl;

                        noalias(rLeftHandSideMatrix) += data.vol*density*prod(data.DN_DX, trans(data.DN_DX));
                        noalias(rLeftHandSideMatrix) += data.vol*derivative*outer_prod(DNV, trans(DNV));

                        noalias(rLaplacianMatrix) += data.vol*density*prod(data.DN_DX, trans(data.DN_DX));

                        noalias(rRightHandSideVector) = -prod(rLaplacianMatrix, data.phis); 
                }               
            }        
           
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


            if(rLaplacianMatrix.size1() != 2*NumNodes || rLaplacianMatrix.size2() != 2*NumNodes)
                rLaplacianMatrix.resize(2*NumNodes,2*NumNodes,false);
            rLaplacianMatrix.clear(); 
            
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

            Matrix laplacian_positive = ZeroMatrix(NumNodes,NumNodes);
            Matrix laplacian_negative = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_penalty_positive = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_penalty_negative = ZeroMatrix(NumNodes,NumNodes);
            bounded_matrix<double, 2, 1 > n;
            n(0,0)=0;
            n(1,0)=1;

            Matrix test=prod(data.DN_DX,n);
            // for(unsigned int i=0;i<test.size1();i++){
            //     for(unsigned int j=0;j<test.size2();j++){
            //         std::cout<<test(i,j)<<", ";
            //     }
            //     std::cout<<std::endl;
            // }
            // std::cout<<std::endl;
            for(unsigned int i=0; i<nsubdivisions; ++i)
            {
                this->ComputeDensityOnSplitElement(density,derivative,velocity,PartitionsSign[i],rCurrentProcessInfo);
                DNV = prod(data.DN_DX, velocity);
                if(PartitionsSign[i] > 0){
                    noalias(lhs_positive) += Volumes[i]*density*prod(data.DN_DX, trans(data.DN_DX));
                    noalias(lhs_positive) += Volumes[i]*derivative*outer_prod(DNV, trans(DNV));

                    noalias(laplacian_positive) += Volumes[i]*density*prod(data.DN_DX, trans(data.DN_DX)); 

                    noalias(lhs_penalty_positive) += Volumes[i] * prod(test,trans(test));
                }
                else{
                    noalias(lhs_negative) += Volumes[i]*density*prod(data.DN_DX, trans(data.DN_DX));
                    noalias(lhs_negative) += Volumes[i]*derivative*outer_prod(DNV, trans(DNV));

                    noalias(laplacian_negative) += Volumes[i]*density*prod(data.DN_DX, trans(data.DN_DX));

                    noalias(lhs_penalty_negative) += Volumes[i] * prod(test,trans(test));
                }
            }
            
            // double penalty = rCurrentProcessInfo[INITIAL_PENALTY];

            //also next version works - NON SYMMETRIC - but it does not require a penalty
//                 array_1d<double,Dim> n = prod(data.DN_DX,data.distances); //rCurrentProcessInfo[VELOCITY]; 
//                 n /= norm_2(n);
//                 bounded_matrix<double,Dim,Dim> nn = outer_prod(n,n);
//                 bounded_matrix<double,NumNodes,Dim> tmp = prod(data.DN_DX,nn);
//                 bounded_matrix<double,NumNodes,NumNodes> constraint = data.vol*prod(tmp, trans(data.DN_DX));
//                                 
//                 bounded_matrix<double,Dim,Dim> P = IdentityMatrix(Dim,Dim) - nn;
//                 noalias(tmp) = prod(data.DN_DX,P);
//                 bounded_matrix<double,NumNodes,NumNodes> tangent_constraint = /*1e3**/data.vol*prod(tmp, trans(data.DN_DX));
                if(kutta_element)
                {
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        for(unsigned int j=0; j<NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i,j)            =  lhs_positive(i,j); 
                            rLeftHandSideMatrix(i,j+NumNodes)   =  0.0; 
                            
                            rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  lhs_negative(i,j); 
                            rLeftHandSideMatrix(i+NumNodes,j)          =  0.0;

                            rLaplacianMatrix(i,j)            =  laplacian_positive(i,j); 
                            rLaplacianMatrix(i,j+NumNodes)   =  0.0; 
                            
                            rLaplacianMatrix(i+NumNodes,j+NumNodes) =  laplacian_negative(i,j); 
                            rLaplacianMatrix(i+NumNodes,j) = 0.0;
                        }
                    }
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        if(data.distances[i]<0)
                        {
                            std::cout << "LOWER WAKE NODE = " << GetGeometry()[i].Id()  << std::endl;
                            for(unsigned int j=0; j<NumNodes; ++j)
                            {
                                rLeftHandSideMatrix(i,j)          = lhs_positive(i,j); 
                                rLeftHandSideMatrix(i,j+NumNodes) = -lhs_positive(i,j); 

                                rLaplacianMatrix(i,j)          = laplacian_positive(i,j); 
                                rLaplacianMatrix(i,j+NumNodes) = -laplacian_positive(i,j);
                            }
                        }
                    }
                    
                    //side2 -assign constraint only on the NEGATIVE_POTENTIAL dofs
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {                            
                        if(data.distances[i]>0)
                        {   
                            std::cout << "UPPER WAKE NODE = " << GetGeometry()[i].Id()  << std::endl;
                            for(unsigned int j=0; j<NumNodes; ++j)
                            {
                                rLeftHandSideMatrix(i+NumNodes,j+NumNodes) = lhs_negative(i,j);
                                //rLeftHandSideMatrix(i+NumNodes,j) = -lhs_negative(i,j);

                                rLaplacianMatrix(i+NumNodes,j+NumNodes) = laplacian_negative(i,j);
                                //rLaplacianMatrix(i+NumNodes,j) = -laplacian_negative(i,j);
                            }
                        }
                    }
                }
                else
                {
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        for(unsigned int j=0; j<NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i,j)                   =  lhs_positive(i,j); 
                            rLeftHandSideMatrix(i,j+NumNodes)          =  0.0; 
                            
                            rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  lhs_negative(i,j); 
                            rLeftHandSideMatrix(i+NumNodes,j)          =  0.0;

                            rLaplacianMatrix(i,j)            =  laplacian_positive(i,j); 
                            rLaplacianMatrix(i,j+NumNodes)   =  0.0; 
                            
                            rLaplacianMatrix(i+NumNodes,j+NumNodes) =  laplacian_negative(i,j); 
                            rLaplacianMatrix(i+NumNodes,j)          =  0.0; 
                        }
                    }
                    
                
                    //side1  -assign constraint only on the NEGATIVE_POTENTIAL dofs
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        if(data.distances[i]<0)
                        {
                            for(unsigned int j=0; j<NumNodes; ++j)
                            {
                                rLeftHandSideMatrix(i,j)          = lhs_positive(i,j); 
                                rLeftHandSideMatrix(i,j+NumNodes) = -lhs_positive(i,j); 

                                rLaplacianMatrix(i,j)          = laplacian_positive(i,j); 
                                rLaplacianMatrix(i,j+NumNodes) = -laplacian_positive(i,j);
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
                                rLeftHandSideMatrix(i+NumNodes,j+NumNodes) = lhs_negative(i,j);
                                rLeftHandSideMatrix(i+NumNodes,j) = -lhs_negative(i,j);

                                rLaplacianMatrix(i+NumNodes,j+NumNodes) = laplacian_negative(i,j);
                                rLaplacianMatrix(i+NumNodes,j) = -laplacian_negative(i,j);
                            }
                        }
                    }
                }
            Vector split_element_values(NumNodes*2);
            GetValuesOnSplitElement(split_element_values, data.distances);
            noalias(rRightHandSideVector) = -prod(rLaplacianMatrix,split_element_values);
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
            KRATOS_THROW_ERROR(std::logic_error, "CompressiblePotentialFlowElement found with Id 0 or negative", "")
        }

        if (this->GetGeometry().Area() <= 0)
        {
            std::cout << "error on CompressiblePotentialFlowElement -> " << this->Id() << std::endl;
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

        if (rVariable == PRESSURE)//actually it computes pressure coefficient cp
        {
            double cp = 0.0;

            bool active = true;
            if ((this)->IsDefined(ACTIVE))
                active = (this)->Is(ACTIVE);

            if(active && !this->Is(MARKER))//normal element
            {
                const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
                const double gamma = rCurrentProcessInfo[LAMBDA];
                const double a = rCurrentProcessInfo[SOUND_VELOCITY];

                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                //gather nodal data
                for(unsigned int i=0; i<NumNodes; i++)
                {
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                }

                const array_1d<double,Dim> v = prod(trans(data.DN_DX), data.phis);

                const double v_norm2 = inner_prod(v,v);

                const double base = 1 + (gamma -1)*vinfinity_norm2*(1-v_norm2/vinfinity_norm2)/(2*a*a);

                const double exponent = gamma/(gamma -1);

                cp = 2*a*a*(pow(base,exponent)-1)/(gamma*vinfinity_norm2);
                //double cp1 = (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

                // double cp1 = 2*a*a*(pow(base,exponent)-1)/(gamma*vinfinity_norm2);
                //cp =  (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

                // std::cout << "cp =" << cp << std::endl;
                // std::cout << "Mach =" << vinfinity_norm2/(a*a) << std::endl;
                // std::cout << "cp1 =" << cp1 << std::endl;
            }
            else if(this->Is(MARKER) && active==true)//wake element
            {
                const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
                const double gamma = rCurrentProcessInfo[LAMBDA];
                const double a = rCurrentProcessInfo[SOUND_VELOCITY];

                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                array_1d<double,NumNodes> distances;
                GetWakeDistances(distances);

                //taking only positive part
                // for (unsigned int i = 0; i < NumNodes; i++)
                // {
                //     if(distances[i] > 0)
                //         data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                //     else
                //         data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                // }

                //negative part - sign is opposite to the previous case
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    if(distances[i] < 0)                    
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);                    
                    else                    
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);                    
                }

                const array_1d<double,Dim> v = prod(trans(data.DN_DX), data.phis);

                const double v_norm2 = inner_prod(v,v);
                const double base = 1 + (gamma -1)*vinfinity_norm2*(1-v_norm2/vinfinity_norm2)/(2*a*a);
                const double exponent = gamma/(gamma -1);

                cp = 2*a*a*(pow(base,exponent)-1)/(gamma*vinfinity_norm2);                
                
                //cp =  (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));
            }

            rValues[0] = cp;

        }
        if (rVariable == DENSITY)
        {
            double density = 0.0;

            bool active = true;
            if ((this)->IsDefined(ACTIVE))
                active = (this)->Is(ACTIVE);

            if(active && !this->Is(MARKER))//normal element
            {
                const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
                //const double densityinfinity = rCurrentProcessInfo[DENSITY];
                const double densityinfinity = GetProperties().GetValue(DENSITY);
                const double gamma = rCurrentProcessInfo[LAMBDA];
                const double a = rCurrentProcessInfo[SOUND_VELOCITY];

                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
                
                //gather nodal data
                for(unsigned int i=0; i<NumNodes; i++)
                {
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                }

                //compute local velocity
                const array_1d<double,Dim> v = prod(trans(data.DN_DX), data.phis);

                const double v_norm2 = inner_prod(v,v);
                const double base = 1 + (gamma -1)*vinfinity_norm2*(1-v_norm2/vinfinity_norm2)/(2*a*a);
                const double exponent = 1/(gamma -1);

                density = densityinfinity*pow(base,exponent); 
            }
            else if(this->Is(MARKER) && active==true)//wake element
            {
                const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
                //const double densityinfinity = rCurrentProcessInfo[DENSITY];
                const double densityinfinity = GetProperties().GetValue(DENSITY);
                const double gamma = rCurrentProcessInfo[LAMBDA];
                const double a = rCurrentProcessInfo[SOUND_VELOCITY];

                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                array_1d<double,NumNodes> distances;
                GetWakeDistances(distances);

                //taking only positive part
                // for (unsigned int i = 0; i < NumNodes; i++)
                // {
                //     if(distances[i] > 0)
                //         data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                //     else
                //         data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                // }

                //negative part - sign is opposite to the previous case
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    if(distances[i] < 0)                    
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);                    
                    else                    
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);                    
                }

                const array_1d<double,Dim> v = prod(trans(data.DN_DX), data.phis);

                const double v_norm2 = inner_prod(v,v);
                const double base = 1 + (gamma -1)*vinfinity_norm2*(1-v_norm2/vinfinity_norm2)/(2*a*a);
                const double exponent = 1/(gamma -1);

                density = densityinfinity*pow(base,exponent);          
            }

            rValues[0] = density;

        }
        
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        if(rValues.size() != 1) rValues.resize(1);

        if (rVariable == VELOCITY)
        {
            bool active = true;
            if ((this)->IsDefined(ACTIVE))
                active = (this)->Is(ACTIVE);

            array_1d<double,3> v(3,0.0);
            if(this->IsNot(MARKER) && active==true)
            {
                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                //gather nodal data
                for(unsigned int i=0; i<NumNodes; i++)
                {
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                }

                array_1d<double,Dim> vaux = -prod(trans(data.DN_DX), data.phis);
                
                for(unsigned int k=0; k<Dim; k++) v[k] = vaux[k];
            }
            else if(this->Is(MARKER) && active==true)
            {
                ElementalData<NumNodes,Dim> data;
                
                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                array_1d<double,NumNodes> distances;
                GetWakeDistances(distances);

                //taking only positive part
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    if(distances[i] > 0)
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                    else
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                }

                array_1d<double,Dim> vaux = -prod(trans(data.DN_DX), data.phis);
                double vupnorm = inner_prod(vaux,vaux);

                //taking only negative part
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    if(distances[i] < 0)
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                    else
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                }

                array_1d<double,Dim> vtest = prod(trans(data.DN_DX), data.phis);
                double vdownnorm = inner_prod(vtest,vtest);

                if( abs(vupnorm - vdownnorm) > 0.1)
                {
                   std::cout << "WAKE CONDITION NOT FULFILLED ELEMENT = " << this->Id()  << std::endl; 
                }
                
                for(unsigned int k=0; k<Dim; k++) v[k] = vaux[k];
            }



            rValues[0] = v;
        }
        else if (rVariable == NORMAL)
        {
            rValues[0] = this->GetValue(NORMAL);
        }
        else if (rVariable == VELOCITY_LAPLACIAN)
        {
            bool active = true;
            if ((this)->IsDefined(ACTIVE))
                active = (this)->Is(ACTIVE);

            array_1d<double,3> v = ZeroVector();
            if(this->Is(MARKER) && active==true)
            {
                ElementalData<NumNodes,Dim> data;
                
                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                array_1d<double,NumNodes> distances;
                GetWakeDistances(distances);

                //taking only negative part
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    if(distances[i] < 0)
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                    else
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                }

                array_1d<double,Dim> vaux = prod(trans(data.DN_DX), data.phis);
                
                for(unsigned int k=0; k<Dim; k++) v[k] = vaux[k];
            }



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
        buffer << "CompressiblePotentialFlowElement #" << Id();
        return buffer.str();
    }

/// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CompressiblePotentialFlowElement #" << Id();
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
    void ComputeDensityOnSplitElement(  double& density, 
                                        double& derivative, 
                                        array_1d<double,Dim>& velocity, 
                                        double& rPartitionsSign, 
                                        const ProcessInfo& rCurrentProcessInfo)
    {
        //This function computes density, velocity and the derivative of the density w.r.t. the square of the local velocity
        //for a split element

        const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        //const double densityinfinity = rCurrentProcessInfo[DENSITY];
        const double densityinfinity = GetProperties().GetValue(DENSITY);
        const double gamma = rCurrentProcessInfo[LAMBDA];
        const double a = rCurrentProcessInfo[SOUND_VELOCITY];

        const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

        ElementalData<NumNodes,Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        array_1d<double,NumNodes> distances;
        GetWakeDistances(distances);

        if(rPartitionsSign > 0)//taking positive part
        {
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                else
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
            }
        }
        else if(rPartitionsSign < 0)//taking negative part
        {
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)                    
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);                    
                else                    
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);                    
            }
        }       

        velocity = prod(trans(data.DN_DX), data.phis);

        const double v_norm2 = inner_prod(velocity,velocity);
        const double base = 1 + (gamma -1)*vinfinity_norm2*(1-v_norm2/vinfinity_norm2)/(2*a*a);

        density = densityinfinity*pow(base,1/(gamma -1));
        derivative = -densityinfinity*pow(base,(2 - gamma)/(gamma -1))/(2*a*a);
    }

    void ComputeDensity(double& density, 
                        double& derivative, 
                        array_1d<double,Dim>& velocity, 
                        const ProcessInfo& rCurrentProcessInfo)
    {
        //This function computes density, velocity and the derivative of the density w.r.t. the square of the local velocity
        //for a normal element
        
        const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        //const double densityinfinity = rCurrentProcessInfo[DENSITY];
        const double densityinfinity = GetProperties().GetValue(DENSITY);
        const double gamma = rCurrentProcessInfo[LAMBDA];
        const double a = rCurrentProcessInfo[SOUND_VELOCITY];

        const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);
        
        // Check that all required variables have been read correctly
        if(vinfinity_norm2 == 0 || densityinfinity == 0 || gamma == 0 || a==0 )
            KRATOS_ERROR << "One input variable is 0. Check if the variables are correctly read.";

        

        ElementalData<NumNodes,Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
        
        //gather nodal data
        for(unsigned int i=0; i<NumNodes; i++)
        {
            data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
        }
        
        //compute local velocity
        velocity = prod(trans(data.DN_DX), data.phis);

        double v_norm2 = inner_prod(velocity,velocity);
        if(v_norm2/a > 0.94)
        {
            //std::cout << "Local mach larger than 0.94. Using density correction" << std::endl;
            v_norm2 = 0.94*a;
        }

        const double base = 1 + (gamma -1)*vinfinity_norm2*(1-v_norm2/vinfinity_norm2)/(2*a*a);
        if(base > 0)
        {
            density = densityinfinity*pow(base,1/(gamma -1));
            derivative = -densityinfinity*pow(base,(2 - gamma)/(gamma -1))/(2*a*a); 
        }
        else
        {
            std::cout << "Base smaller than 0. Using density correction" << std::endl;
            std::cout << "vinfinity_norm2 =" << vinfinity_norm2 << std::endl;
            std::cout << "v_norm2 =" << v_norm2 << std::endl;
            std::cout << "base =" << base << std::endl;
            density = densityinfinity*0.00001;
            derivative = -densityinfinity*pow(densityinfinity*0.00001,(2 - gamma)/(gamma -1))/(2*a*a); 
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
        // calculate shape functions
        
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
            std::cout << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << this->Id() << std::endl;
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

}; // Class CompressiblePotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
