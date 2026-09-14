
import numpy as np
import KratosMultiphysics as KM
DOF_MAP = {
    "DISPLACEMENT_X": KM.DISPLACEMENT_X,
    "DISPLACEMENT_Y": KM.DISPLACEMENT_Y,
    "DISPLACEMENT_Z": KM.DISPLACEMENT_Z,
    "ROTATION_X": KM.ROTATION_X,
    "ROTATION_Y": KM.ROTATION_Y,
    "ROTATION_Z": KM.ROTATION_Z,
}
TRANSLATION_NAMES = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")
ROTATION_NAMES = ("ROTATION_X", "ROTATION_Y", "ROTATION_Z")

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return ApplyRbe2Process(Model, settings["Parameters"])

class ApplyRbe2Process(KM.Process):
    """
    Generates a RBE2-coupling with a master and a set of slave nodes (using LinearMasterSlaveConstraint) 
    Expected parameters:
        model_part_name   : name of the model part
        master_sub_model_part : SubModelPart-Name with the master node
        slave_sub_model_part  : SubModelPart-Name with the slave nodes
        constraint_id_start   : MasterSlaveConstraint-Id
        constrained dofs: DOFs that should be coupled
    """
    def __init__(self, model, settings):
        self.echo_level = settings["echo_level"].GetInt() if settings.Has("echo_level") else 0 
        super().__init__() 
        default_settings = KM.Parameters("""{
            "model_part_name"       : "",
            "master_sub_model_part" : "",
            "slave_sub_model_part"  : "",
            "constraint_id_start"   : 1, 
            "constrained_dofs"      : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z",
                                         "ROTATION_X","ROTATION_Y","ROTATION_Z"]
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.master_mp = model[settings["master_sub_model_part"].GetString()]
        self.slave_mp = model[settings["slave_sub_model_part"].GetString()]
        self.constraint_id_start = settings["constraint_id_start"].GetInt()

        requested = settings["constrained_dofs"].GetStringArray()
        invalid = set(requested) - set(DOF_MAP) #DOF must be one of the list above --> otherwise invalid =true
        if invalid:
            raise Exception(f"Invalid input: {invalid}")
        
        self.dof_names = requested #string
        self.dof_vars = [DOF_MAP[name] for name in self.dof_names] #Kratos-Variable (for example KM.DISPLACEMENT_X) 
        self.n_dofs = len(self.dof_vars)
        self.idx_map = {name: i for i, name in enumerate(self.dof_names)} #DOFs --> 0,1,2

        self.couple_rotation = (
            any(n in self.idx_map for n in TRANSLATION_NAMES)
            and any(n in self.idx_map for n in ROTATION_NAMES)
        )
    def ExecuteInitialize(self):
        assert self.master_mp.NumberOfNodes() == 1, "RBE3 needs exactly 1 master node"
        master_node = next(iter(self.master_mp.Nodes))
        x_m = np.array([master_node.X, master_node.Y, master_node.Z])
        master_dofs = [master_node.GetDof(v) for v in self.dof_vars]

        constant = KM.Vector(self.n_dofs, 0.0)
        cid = self.constraint_id_start

        for slave_node in self.slave_mp.Nodes:
            relation_matrix = KM.Matrix(self._build_relation_matrix( slave_node, x_m))
            slave_dofs = [slave_node.GetDof(v) for v in self.dof_vars]

            self.model_part.CreateNewMasterSlaveConstraint(
                "LinearMasterSlaveConstraint", cid,
                master_dofs, slave_dofs, relation_matrix, constant)
            cid += 1 

    def _build_relation_matrix(self, slave_node, x_m):      
        """
        builts the relation matrix between master and slave nodes
        expected input: slave node and x_m (vector to master)
        """ 
        T = np.eye(self.n_dofs)

        if not self.couple_rotation:
            return T
        r = np.array([slave_node.X, slave_node.Y, slave_node.Z]) - x_m
        skew_r = self._skew(r)

        for t_i, t_name in enumerate(TRANSLATION_NAMES):
            if t_name not in self.idx_map:
                continue
            row = self.idx_map[t_name]
            for r_i, r_name in enumerate(ROTATION_NAMES):
                if r_name not in self.idx_map:
                    continue
                col = self.idx_map[r_name]
                T[row, col] = -skew_r[t_i, r_i] 
                #(-ry)·θx + rx·θy + 0·θz  =  rx·θy - ry·θx
        return T
          
    @staticmethod
    def _skew(r):
        """
        skew matrix: cross product r × θ
        describes the rotational part of the slaves
        will be embedded in the 6x6 relation matrix
        expected input: r (vector master-slave)
        """
        return np.array([
            [0.0, -r[2], r[1]],
            [r[2], 0.0, -r[0]],
            [-r[1], r[0], 0.0] 
        ])
