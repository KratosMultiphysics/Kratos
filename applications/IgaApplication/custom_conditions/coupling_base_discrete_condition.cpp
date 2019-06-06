//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "custom_conditions/coupling_base_discrete_condition.h"


namespace Kratos
{

    /**
    * Initialization of condition
    * Initialize sets the basis vectors of the two patches for the undeformed system:
    * mg1_0_master, mg2_0_master, mg3_0_master
    * mg1_0_slave, mg2_0_slave, mg3_0_slave
    */
    void CouplingBaseDiscreteCondition::Initialize()

    {
        KRATOS_TRY
        const unsigned int number_of_points = GetGeometry().size();
        const unsigned int working_space_dimension = 3;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();
        const unsigned int local_space_dimension = 2;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().LocalSpaceDimension();

                                                     //calculate basis vectors for MASTER Patch
        Matrix DN_De_Master = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        Matrix JMaster = ZeroMatrix(working_space_dimension, local_space_dimension);
        JacobianElement(DN_De_Master, JMaster, true);

        array_1d<double, 3> g1_Master;
        array_1d<double, 3> g2_Master;
        array_1d<double, 3> g3_Master;

        g1_Master[0] = JMaster(0, 0);
        g2_Master[0] = JMaster(0, 1);
        g1_Master[1] = JMaster(1, 0);
        g2_Master[1] = JMaster(1, 1);
        g1_Master[2] = JMaster(2, 0);
        g2_Master[2] = JMaster(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(g3_Master, g1_Master, g2_Master);
        g3_Master = g3_Master / norm_2(g3_Master);

        mg1_0_master = g1_Master;
        mg2_0_master = g2_Master;
        mg3_0_master = g3_Master;

        //calculate basis vectors for SLAVE Patch
        array_1d<double, 3> g1_Slave;
        array_1d<double, 3> g2_Slave;
        array_1d<double, 3> g3_Slave;
        const Matrix DN_De_Slave = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE);
        Matrix JSlave = ZeroMatrix(working_space_dimension, local_space_dimension);
        JacobianElement(DN_De_Slave, JSlave, false);

        g1_Slave[0] = JSlave(0, 0);
        g2_Slave[0] = JSlave(0, 1);
        g1_Slave[1] = JSlave(1, 0);
        g2_Slave[1] = JSlave(1, 1);
        g1_Slave[2] = JSlave(2, 0);
        g2_Slave[2] = JSlave(2, 1);

        MathUtils<double>::CrossProduct(g3_Slave, g1_Slave, g2_Slave);

        g3_Slave = g3_Slave / norm_2(g3_Slave);

        mg1_0_slave = g1_Slave;
        mg2_0_slave = g2_Slave;
        mg3_0_slave = g3_Slave;

        KRATOS_CATCH("")
    }

    /**
    * ONLY NEEDED FOR COUPLING OF MORE PATCHES
    * MappingGeometricToParameterMasterElement calculates the J-tilde for the mapping from
    Geometric to Parameter Space. This paramater is needed for all condition
    integrations on edges.
    * For coupling condition which contain more patches is needed a function which
    calculates this mapping and the Jacobian for only one patch. As integration
    is done on the edge of the master element the mapping parameter only has to be calculated on the master element
    *
    * @param[in] DN_De_Master derivatives of shape functions of master patch.
    * @param[out] JGeometricToParameter Mapping parameter for Geometric Space to
    Parameter Space
    *
    * @see JacobianElement
    */
    void CouplingBaseDiscreteCondition::MappingGeometricToParameterMasterElement(const Matrix& DN_De_Master,
        const array_1d<double, 2>& Tangents,
        double& JGeometricToParameter)
    {
        Matrix J;
        JacobianElement(DN_De_Master, J, true);

        array_1d<double, 3> g1;
        array_1d<double, 3> g2;
        array_1d<double, 3> g3;

        g1[0] = J(0, 0);
        g2[0] = J(0, 1);
        g1[1] = J(1, 0);
        g2[1] = J(1, 1);
        g1[2] = J(2, 0);
        g2[2] = J(2, 1);
        //basis vector g3
        //CrossProduct(g3, g1, g2);

        array_1d<double, 3> temp = g1 * Tangents[0] + g2 * Tangents[1];
        // g1*localDerivatives[0] + g2*localDerivatives[1];//(g1 / norm_2(g1))*localDerivatives[0] + (g2 / norm_2(g2))*localDerivatives[1];
        JGeometricToParameter = norm_2(temp);

    }

    /**
    * ONLY NEEDED FOR COUPLING BETWEEN PATCHES
    * MappingGeometricToParameterMaster calculates the J-tilde for the mapping from
    Geometry to Parameter Space. This paramater is needed for all continuity condition
    integrations on curves. As integration is evaluated on the master curve the J-tilde
    is not needed for the slave.
    * The function needs the KRATOS-variables SHAPE_FUNCTIONS_LOCAL_DERIVATIVES and the
    TANGENTS.
    *
    * @param[out] JGeometricToParameter Mapping parameter for Geometric Space to
    Parameter Space
    *
    * @see JacobianElement
    */
    void CouplingBaseDiscreteCondition::MappingGeometricToParameterOnMasterCurve(double& JGeometricToParameter)
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Vector& Tangents = this->GetValue(TANGENTS);

        Matrix J;
        JacobianElement(DN_De, J, true);

        array_1d<double, 3> g1;
        array_1d<double, 3> g2;
        array_1d<double, 3> g3;

        g1[0] = J(0, 0);
        g2[0] = J(0, 1);
        g1[1] = J(1, 0);
        g2[1] = J(1, 1);
        g1[2] = J(2, 0);
        g2[2] = J(2, 1);
        //basis vector g3

        array_1d<double, 3> temp = g1 * Tangents[0] + g2 * Tangents[1];
        // g1*localDerivatives[0] + g2*localDerivatives[1];//(g1 / norm_2(g1))*localDerivatives[0] + (g2 / norm_2(g2))*localDerivatives[1];
        JGeometricToParameter = norm_2(temp);
    }

    /**
    * ONLY NEEDED FOR COUPLING OF 2 PATCHES
    * JacobianElement calculates Jacobian for the mapping between Geometry Space and
    Parameter Space. The flag Master asks wether the first points of the geometry
    shall be used (TRUE). Or the points that come after that points
    * For coupling condition which contain more patches is needed a function which
    calculates this the Jacobian for only one patch.
    *
    * @param[in] DN_De derivatives of shape functions of master or slave patch. Only ONE PATCH!!
    * @param[out] Jacobian Mapping parameter for Geometric Space to
    *             Parameter Space
    *
    * @see Jacobian
    */
    void CouplingBaseDiscreteCondition::JacobianElement(
        const Matrix& DN_De,
        Matrix& Jacobian, 
        const bool Master)
    {
        const unsigned int number_of_points = GetGeometry().size();
        const unsigned int working_space_dimension = 3;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();
        const unsigned int local_space_dimension = 2;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().LocalSpaceDimension();

        Jacobian.resize(working_space_dimension, local_space_dimension);
        Jacobian.clear();

        int shift = 0;
        if (!Master)
            shift = number_of_points - DN_De.size1();

        KRATOS_WATCH(GetGeometry().size())
            KRATOS_WATCH(DN_De.size1())
            KRATOS_WATCH(DN_De.size2())

        for (unsigned int i = 0; i < DN_De.size1(); i++)
        {
            for (unsigned int k = 0; k < working_space_dimension; k++)
            {
                for (unsigned int m = 0; m < local_space_dimension; m++)
                {
                    Jacobian(k, m) += (GetGeometry()[i + shift]).Coordinates()[k] * DN_De(i, m);
                }
            }
        }
    }

    void CouplingBaseDiscreteCondition::CaculateRotation(const Matrix &ShapeFunctionDerivatives,
        Vector &Phi_r, Matrix &Phi_rs, array_1d<double, 2> &Phi, array_1d<double, 3> &TrimTangent,
        const Vector &Tangents, const bool Master)
    {
        KRATOS_TRY

            const unsigned int number_of_points = ShapeFunctionDerivatives.size1();
        array_1d<double, 3> g10, g20, g30;

        if (Master)
        {
            g10 = mg1_0_master;
            g20 = mg2_0_master;
            g30 = mg3_0_master;
        }
        else
        {
            g10 = mg1_0_slave;
            g20 = mg2_0_slave;
            g30 = mg3_0_slave;
        }

        Matrix J;
        JacobianElement(ShapeFunctionDerivatives, J, Master);

        array_1d<double, 3> g1, g2, g3;

        g1[0] = J(0, 0);
        g2[0] = J(0, 1);
        g1[1] = J(1, 0);
        g2[1] = J(1, 1);
        g1[2] = J(2, 0);
        g2[2] = J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(g3, g1, g2);
        g3 = g3 / norm_2(g3);

        // t1 normal to trim, t2 tangential to trim
        array_1d<double, 3> t2 = Tangents(0)*g10 + Tangents(1)*g20;  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TrimTangent = t2;
        array_1d<double, 3> t1;
        MathUtils<double>::CrossProduct(t1, t2, g30);
        t2 = t2 / norm_2(t2);
        t1 = t1 / norm_2(t1);

        // computation of the a3 displacement
        array_1d<double, 3> w = g3 - g30;
        array_1d<double, 3> SinusOmegaVector;
        MathUtils<double>::CrossProduct(SinusOmegaVector, g30, w);

        array_1d<double, 2> SinusOmega;
        SinusOmega(0) = inner_prod(SinusOmegaVector, t2);
        SinusOmega(1) = inner_prod(SinusOmegaVector, t1);

        array_1d<double, 3> Omega;
        if (SinusOmega(0) > 1.0)
            SinusOmega(0) = 0.999999;
        if (SinusOmega(1) > 1.0)
            SinusOmega(1) = 0.999999;
        Omega(0) = asin(SinusOmega(0));
        Omega(1) = asin(SinusOmega(1));

        //array_1d<double, 2> Phi;
        Phi(0) = Omega(0);
        Phi(1) = Omega(1);

        //variation of the a3
        array_1d<double, 3> t3 = g3;
        array_1d<double, 3> tilde_t3; //g3
        MathUtils<double>::CrossProduct(tilde_t3, g1, g2);
        double Length_t3 = norm_2(tilde_t3);

        for (unsigned int n = 0; n < number_of_points; n++)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                //variations of the basis vectors
                array_1d<double, 3> a1_r;
                a1_r.clear();
                array_1d<double, 3> a2_r;
                a2_r.clear();

                a1_r(i) = ShapeFunctionDerivatives(n, 0);
                a2_r(i) = ShapeFunctionDerivatives(n, 1);
                //variation of the non normalized local vector
                array_1d<double, 3> tilde_3_r = CrossProduct(a1_r, g2) + CrossProduct(g1, a2_r);
                double line_t3_r = inner_prod(t3, tilde_3_r);
                array_1d<double, 3> t3_r = tilde_3_r / Length_t3 - line_t3_r * t3 / Length_t3;
                array_1d<double, 3> SinusOmega_r;
                MathUtils<double>::CrossProduct(SinusOmega_r, g30, t3_r);
                Phi_r(n * 3 + i) = 1.0 / sqrt(1.0 - pow(SinusOmega(0), 2))*inner_prod(SinusOmega_r, t2);
                // if needed at some point:
                //Phi_r_2(i * 3 + j) = 1.0 / sqrt(1.0 - pow(SinusOmega(1), 2))*inner_prod(SinusOmega_r, t1);
            }
        }
        for (unsigned int n = 0; n < number_of_points; n++)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                //variations of the basis vectors
                array_1d<double, 3> a1_r_n;
                a1_r_n.clear();
                array_1d<double, 3> a2_r_n;
                a2_r_n.clear(); /*[0] = 0.0; a2_r_n[1] = 0.0; a2_r_n[2] = 0.0;*/

                a1_r_n(i) = ShapeFunctionDerivatives(n, 0);
                a2_r_n(i) = ShapeFunctionDerivatives(n, 1);

                //variation of the non normalized local vector
                array_1d<double, 3> tilde_3_r_n = CrossProduct(a1_r_n, g2) + CrossProduct(g1, a2_r_n);
                double line_t3_r_n = inner_prod(t3, tilde_3_r_n);
                array_1d<double, 3> t3_r_n = tilde_3_r_n / Length_t3 - line_t3_r_n * t3 / Length_t3;
                array_1d<double, 3> SinusOmega_r_n;
                MathUtils<double>::CrossProduct(SinusOmega_r_n, g30, t3_r_n);

                for (unsigned int m = 0; m < number_of_points; m++)
                {
                    for (unsigned int j = 0; j < 3; j++)
                    {
                        //variations of the basis vectors
                        array_1d<double, 3> a1_r_m;
                        a1_r_m.clear();
                        array_1d<double, 3> a2_r_m;
                        a2_r_m.clear();

                        a1_r_m(j) = ShapeFunctionDerivatives(m, 0);
                        a2_r_m(j) = ShapeFunctionDerivatives(m, 1);


                        //variation of the non normalized local vector
                        array_1d<double, 3> tilde_3_r_m = CrossProduct(a1_r_m, g2) + CrossProduct(g1, a2_r_m);
                        double line_t3_r_m = inner_prod(t3, tilde_3_r_m);
                        array_1d<double, 3> t3_r_m = tilde_3_r_m / Length_t3 - line_t3_r_m * t3 / Length_t3;
                        array_1d<double, 3> SinusOmega_r_m = CrossProduct(g30, t3_r_m);

                        array_1d<double, 3> tilde_t3_rs = CrossProduct(a1_r_n, a2_r_m) + CrossProduct(a1_r_m, a2_r_n);
                        double line_t3_rs = inner_prod(t3_r_m, tilde_3_r_n) + inner_prod(t3, tilde_t3_rs);
                        array_1d<double, 3> t3_rs = (tilde_t3_rs*Length_t3 - line_t3_r_m * tilde_3_r_n) / pow(Length_t3, 2)
                            - line_t3_rs * t3 / Length_t3 - line_t3_r_n * (t3_r_m * Length_t3 - line_t3_r_m * t3) / pow(Length_t3, 2);
                        array_1d<double, 3> SinusOmega_rs = CrossProduct(g30, t3_rs);

                        Phi_rs(n * 3 + i, m * 3 + j) = inner_prod(SinusOmega_rs, t2) / sqrt(1.0 - pow(SinusOmega(0), 2))
                            + inner_prod(SinusOmega_r_m, t2)*inner_prod(SinusOmega_r_n, t2)*SinusOmega(0) / pow(1.0
                                - pow(SinusOmega(0), 2), 1.5);
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    Vector CouplingBaseDiscreteCondition::CrossProduct(const array_1d<double, 2>& v1, const array_1d<double, 2>& v2)
    {
        Vector v3 = ZeroVector(3);
        MathUtils<double>::CrossProduct(v3, v1, v2);
        return v3;
    }


    void CouplingBaseDiscreteCondition::CaculateRotation2(const Matrix &ShapeFunctionDerivatives,
        Vector &Phi_r, Matrix &Phi_rs, array_1d<double, 2> &Phi, array_1d<double, 3> &TrimTangent, const Vector &Tangents, const bool Master)
    {
        //KRATOS_TRY

        //    const unsigned int number_of_points = ShapeFunctionDerivatives.size1();
        //Vector g10, g20, g30;

        //if (Master)
        //{
        //    g10 = mg1_0_master;
        //    g20 = mg2_0_master;
        //    g30 = mg3_0_master;
        //}
        //else
        //{
        //    g10 = mg1_0_slave;
        //    g20 = mg2_0_slave;
        //    g30 = mg3_0_slave;
        //}

        //Matrix J;
        //JacobianElement(ShapeFunctionDerivatives, J, Master);

        //Vector g1, g2, t3;
        //GetBasisVectors(ShapeFunctionDerivatives, g1, g2, t3);

        //// t1 normal to trim, t2 tangential to trim
        //Vector T2 = Tangents(0)*g10 + Tangents(1)*g20;  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ////TrimTangent = T2;
        //Vector T1 = ZeroVector(3);
        //MathUtils<double>::CrossProduct(T1, T2, g30);
        //T2 = T2 / norm_2(T2);
        //T1 = T1 / norm_2(T1);

        //// computation of the a3 displacement
        //Vector w = t3 - g30;
        //Vector omega = ZeroVector(3);
        //MathUtils<double>::CrossProduct(omega, g30, w);

        //double omega_T2 = asin(inner_prod(omega, T2));
        //double omega_T1 = asin(inner_prod(omega, T1));

        //Phi(0) = omega_T2;
        //Phi(1) = omega_T1;


        //Vector g_3_1 = ZeroVector(3);
        //MathUtils<double>::CrossProduct(g_3_1, g1, g2);
        //double Length_t3 = norm_2(g_3_1);

        //for (unsigned int n = 0; n < number_of_points; n++)
        //{
        //    for (unsigned int i = 0; i < 3; i++)
        //    {
        //        //variations of the basis vectors
        //        Vector t1_r = ZeroVector(3);
        //        Vector t2_r = ZeroVector(3);
        //        t1_r(i) = ShapeFunctionDerivatives(n, 0);
        //        t2_r(i) = ShapeFunctionDerivatives(n, 1);

        //        //variation of the non normalized local vector
        //        Vector t3_r_tilde = CrossProduct(t1_r, g2) + CrossProduct(g1, t2_r);
        //        double t3_r_line = inner_prod(t3, t3_r_tilde);
        //        Vector t3_r = t3_r_tilde / Length_t3 - t3_r_line * t3 / Length_t3;
        //        Vector omega_r = CrossProduct(g30, t3_r);
        //        Phi_r(n * 3 + i) = inner_prod(omega_r, T2) / sqrt(1.0 - pow(omega_T2, 2));
        //        // if needed at some point:
        //        //Phi_r_2(i * 3 + j) = inner_prod(omega_r, T1) / sqrt(1.0 - pow(omega_T1, 2));
        //    }
        //}
        //for (unsigned int n = 0; n < number_of_points; n++)
        //{
        //    for (unsigned int i = 0; i < 3; i++)
        //    {
        //        //variations of the basis vectors
        //        Vector t1_r_n = ZeroVector(3);
        //        Vector t2_r_n = ZeroVector(3);
        //        t1_r_n(i) = ShapeFunctionDerivatives(n, 0);
        //        t2_r_n(i) = ShapeFunctionDerivatives(n, 1);

        //        //variation of the non normalized local vector
        //        Vector t3_r_tilde_n = CrossProduct(t1_r_n, g2) + CrossProduct(g1, t2_r_n);
        //        double t3_r_line_n = inner_prod(t3, t3_r_tilde_n);
        //        Vector t3_r_n = t3_r_tilde_n / Length_t3 - t3_r_line_n * t3 / Length_t3;
        //        Vector omega_r_n = CrossProduct(g30, t3_r_n);

        //        for (unsigned int m = 0; m < number_of_points; m++)
        //        {
        //            for (unsigned int j = 0; j < 3; j++)
        //            {
        //                //variations of the basis vectors
        //                Vector t1_r_m = ZeroVector(3);
        //                Vector t2_r_m = ZeroVector(3);
        //                t1_r_m(j) = ShapeFunctionDerivatives(m, 0);
        //                t2_r_m(j) = ShapeFunctionDerivatives(m, 1);

        //                //variation of the non normalized local vector
        //                Vector t3_r_tilde_m = CrossProduct(t1_r_m, g2) + CrossProduct(g1, t2_r_m);
        //                double t3_r_line_m = inner_prod(t3, t3_r_tilde_m);
        //                Vector t3_r_m = t3_r_tilde_m / Length_t3 - t3_r_line_m * t3 / Length_t3;
        //                Vector omega_r_m = CrossProduct(g30, t3_r_m);

        //                Vector t3_rs_tilde = CrossProduct(t1_r_n, t2_r_m) + CrossProduct(t1_r_m, t2_r_n);
        //                double t3_rs_line = inner_prod(t3_r_m, t3_r_tilde_n) + inner_prod(t3, t3_rs_tilde);
        //                Vector t3_rs = (t3_rs_tilde*Length_t3 - t3_r_line_m * t3_r_tilde_n) / pow(Length_t3, 2)
        //                    - t3_rs_line * t3 / Length_t3 - t3_r_line_n * (t3_r_m * Length_t3 - t3_r_line_m * t3) / pow(Length_t3, 2);
        //                Vector omega_rs = CrossProduct(g30, t3_rs);

        //                Phi_rs(n * 3 + i, m * 3 + j) = inner_prod(omega_rs, T2) / sqrt(1.0 - pow(omega_T2, 2))
        //                    + inner_prod(omega_r_m, T2)*inner_prod(omega_r_n, T2)*omega_T2 / pow(1.0
        //                        - pow(omega_T2, 2), 1.5);
        //            }
        //        }
        //    }
        //}
        //KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    void CouplingBaseDiscreteCondition::CaculateRotationalShapeFunctions(
        Vector &Phi_r, Vector &Phi_r_Lambda, Matrix &Phi_rs, array_1d<double, 2> &Diff_Phi)
    {
        Vector localTrimTangentsMasterVector = this->GetValue(TANGENTS);
        Matrix ShapeFunctionDerivativesMaster = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        int number_of_points_master = ShapeFunctionDerivativesMaster.size1();
        Vector Phi_r_Master = ZeroVector(number_of_points_master * 3);
        Matrix Phi_rs_Master = ZeroMatrix(number_of_points_master * 3, number_of_points_master * 3);
        array_1d<double, 2> Phi_Master;
        array_1d<double, 3> TrimTangentsMaster;
        TrimTangentsMaster.clear();
        CaculateRotation(ShapeFunctionDerivativesMaster, Phi_r_Master, Phi_rs_Master, Phi_Master, TrimTangentsMaster, localTrimTangentsMasterVector, true);


        Vector localTrimTangentsSlaveVector = this->GetValue(TANGENTS_SLAVE);
        Matrix ShapeFunctionDerivativesSlave = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE);
        int number_of_points_slave = ShapeFunctionDerivativesSlave.size1();
        Vector Phi_r_Slave = ZeroVector(number_of_points_slave * 3);
        Matrix Phi_rs_Slave = ZeroMatrix(number_of_points_slave * 3, number_of_points_slave * 3);
        array_1d<double, 2> Phi_Slave;
        array_1d<double, 3> TrimTangentsSlave;
        TrimTangentsSlave.clear();
        CaculateRotation(ShapeFunctionDerivativesSlave, Phi_r_Slave, Phi_rs_Slave, Phi_Slave, TrimTangentsSlave, localTrimTangentsSlaveVector, false);

        bool OppositeDirectionOfTrims = true;
        if (inner_prod(TrimTangentsMaster, TrimTangentsSlave) > 0) // tangents have the same direction (assumption for coordinates system changes)
            OppositeDirectionOfTrims = false;


        if (OppositeDirectionOfTrims)
            Diff_Phi = Phi_Slave + Phi_Master;
        else
            Diff_Phi = Phi_Slave - Phi_Master;

        for (unsigned int i = 0; i < Phi_r_Master.size(); i++)
        {
            Phi_r(i) = Phi_r_Master(i);
            Phi_r_Lambda(i) = Phi_r_Master(i);
        }
        int index = Phi_r_Master.size();
        for (unsigned int i = 0; i < Phi_r_Slave.size(); i++)
        {
            if (OppositeDirectionOfTrims)
                Phi_r(i + index) = Phi_r_Slave(i);
            else
                Phi_r(i + index) = -Phi_r_Slave(i);

            Phi_r_Lambda(i + index) = 0;// Phi_r_Slave(i);
        }

        for (unsigned int i = 0; i < Phi_rs_Master.size1(); i++)
        {
            for (unsigned int j = 0; j < Phi_rs_Master.size2(); j++)
            {
                Phi_rs(i, j) = Phi_rs_Master(i, j);
            }
        }
        int index1 = Phi_rs_Master.size1();
        int index2 = Phi_rs_Master.size2();
        for (unsigned int i = 0; i < Phi_rs_Slave.size1(); i++)
        {
            for (unsigned int j = 0; j < Phi_rs_Slave.size2(); j++)
            {
                if (OppositeDirectionOfTrims)
                    Phi_rs(i + index1, j + index2) = Phi_rs_Slave(i, j);
                else
                    Phi_rs(i + index1, j + index2) = -Phi_rs_Slave(i, j);
            }
        }
    }
} // Namespace Kratos


