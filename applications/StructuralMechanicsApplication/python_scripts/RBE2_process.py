
import numpy as np
import KratosMultiphysics as KM

ALL_DOFS = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z",
            "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return ApplyRbe2Process(Model, settings["Parameters"])

class ApplyRbe2Process(KM.Process):
    """
    @brief 
    Generates an RBE2-type coupling between one master node and a set of slave nodes using LinearMasterSlaveConstraint objects.
    
    @details
    The coupling enforces rigid-body motion between the master node and the slave nodes. 
    The selected degrees of freedom of each slave node are determined by the corresponding motion of the master node.
    
    Expected parameters:
    @param        main_odel_part_name           :      name of the main model part
    @param        master_sub_model_part_name    :      name of the subModelPart with the master node
    @param        slave_sub_model_part_name     :      name of the subModelPart with the slave nodes
    @param        constraint_id                 :      ID for the MasterSlaveConstraint
    @param        constrained_dofs              :      list of degrees of freedom to couple
    """
    def __init__(self, model, settings):
        self.echo_level = settings["echo_level"].GetInt() if settings.Has("echo_level") else 0 
        super().__init__() 
        default_settings = KM.Parameters("""{
            "main_model_part_name"       : "",
            "master_sub_model_part_name" : "",
            "slave_sub_model_part_name"  : "",
            "constraint_id"              : 1,
            "constrained_dofs"           : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z",
                                         "ROTATION_X","ROTATION_Y","ROTATION_Z"]
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["main_model_part_name"].GetString()]
        self.master_mp = model[settings["master_sub_model_part_name"].GetString()]
        self.slave_mp = model[settings["slave_sub_model_part_name"].GetString()]
        self.constraint_id = settings["constraint_id"].GetInt()
        self.constrained_dofs = settings["constrained_dofs"].GetStringArray()

        self.dof_idx = [ALL_DOFS.index(name) for name in self.constrained_dofs]

    def ExecuteInitialize(self):
        if self.master_mp.NumberOfNodes() != 1:
            raise RuntimeError("RBE2 needs exactly 1 master node")
        master_node = next(iter(self.master_mp.Nodes))
        master_vars = [KM.KratosGlobals.GetVariable(name) for name in ALL_DOFS]

        slave_vars = [KM.KratosGlobals.GetVariable(dof) for dof in self.constrained_dofs]
        master_dofs = [master_node.GetDof(v) for v in master_vars]
        x_m = np.array([master_node.X, master_node.Y, master_node.Z])    

        slave_dofs = []
        A_list=[]
        for slave_node in self.slave_mp.Nodes:
            slave_dofs.extend(slave_node.GetDof(v) for v in slave_vars)
            A_list.append(self._build_relation_matrix(slave_node, x_m))

        relation_matrix =KM.Matrix(np.vstack(A_list))
        constant = KM.Vector(len(slave_dofs), 0.0)

        self.model_part.CreateNewMasterSlaveConstraint(
            "LinearMasterSlaveConstraint", self.constraint_id,
            master_dofs, slave_dofs, relation_matrix, constant)

    def _build_relation_matrix(self, slave_node, x_m):      
        """
        builts the relation matrix between master and slave nodes
        expected input: slave node and x_m (vector to master)
        """ 
        A = np.eye(6)
        r = np.array([slave_node.X, slave_node.Y, slave_node.Z]) - x_m
        A[:3, 3:] = -self._skew(r)
        return A[self.dof_idx, :]

          
    @staticmethod
    def _skew(r):
        """
        skew matrix: cross product r x θ 
        describes the rotational part of the slaves will be embedded in the 6x6 relation matrix 
        expected input: r (vector master-slave)
        """
        return np.array([
            [0.0, -r[2], r[1]],
            [r[2], 0.0, -r[0]],
            [-r[1], r[0], 0.0] 
        ])
