
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
    return ApplyRbe3Process(Model, settings["Parameters"])

class ApplyRbe3Process(KM.Process):
    """
    Generates a RBE3-coupling (using LinearMasterSlaveConstraint) between a reference node and a set of connected/independant nodes.
    
    The movement of the reference node results from the movements of the connected nodes.
    Expected parameters: 
        model_part_name   : name of the modelpart
        connected_sub_model_part : SubModelPart-Name with several independent/connected nodes
        reference_sub_model_part  : SubModelPart-Name with the reference node
        constraint_id_start   : first free MasterSlaveConstraint-Id
        constrained dofs: the coupled DOFs
        weight_variable: variable, that defines the weights of the connected nodes. If empty, uniform weight 1.0 for all nodes.
    """
    def __init__(self, model, settings):
        self.echo_level = settings["echo_level"].GetInt() if settings.Has("echo_level") else 0 
        super().__init__()
        default_settings = KM.Parameters("""{
            "model_part_name"       : "",
            "connected_sub_model_part" : "",
            "reference_sub_model_part"  : "",
            "constraint_id_start"   : 1, 
            "constrained_dofs"      : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z",
                                         "ROTATION_X","ROTATION_Y","ROTATION_Z"], 
            "weight_variable"       : ""
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.reference_mp = model[settings["reference_sub_model_part"].GetString()]
        self.connected_mp = model[settings["connected_sub_model_part"].GetString()]
        self.constraint_id_start = settings["constraint_id_start"].GetInt()

        requested = settings["constrained_dofs"].GetStringArray()
        invalid = set(requested) - set(DOF_MAP) 
        if invalid:
            raise Exception(f"Invalid input: {invalid}")
        
        self.dof_names = requested #string
        self.dof_vars = [DOF_MAP[name] for name in self.dof_names] 
        self.n_dofs = len(self.dof_vars)
        self.idx_map = {name: i for i, name in enumerate(self.dof_names)} #DOFs --> 0,1,2..

        has_translation = any(n in self.idx_map for n in TRANSLATION_NAMES)
        has_rotation = any(n in self.idx_map for n in ROTATION_NAMES)
        self.couple_rotation = has_translation and has_rotation
 
        weight_var_name = settings["weight_variable"].GetString()
        self.weight_var = KM.KratosGlobals.GetVariable(weight_var_name) if weight_var_name else None
        
    def ExecuteInitialize(self):
        """
        Defines the relation between the connected nodes and the reference node with respect to the weigths
        """
        if self.reference_mp.NumberOfNodes() != 1:
            raise RuntimeError("RBE3 needs exactly 1 reference node")
        ref_node = next(iter(self.reference_mp.Nodes))
        x_r = np.array([ref_node.X, ref_node.Y, ref_node.Z])
        ref_dofs = [ref_node.GetDof(v) for v in self.dof_vars]

        connected_nodes = list(self.connected_mp.Nodes)
        n = len(connected_nodes)
        if n < 0:
            raise RuntimeError("RBE3 needs at least 1 connected node")

        A_list = [] 
        w_list = []
        connected_dofs = []
        for node in connected_nodes:
            A_i = self._build_relation_matrix(node, x_r) 
            A_list.append(A_i)
            w_list.append(self._get_weight(node))
            connected_dofs.extend(node.GetDof(v) for v in self.dof_vars)

        K = np.zeros((self.n_dofs, self.n_dofs))
        for A_i, w_i in zip(A_list, w_list):
            K += w_i * A_i.T @ A_i 

        K_inv = np.linalg.pinv(K)

        C_blocks = [K_inv @ (w_i * A_i.T) for A_i, w_i in zip(A_list, w_list)]
        relation_matrix_np = np.hstack(C_blocks) 
 
        relation_matrix = KM.Matrix(relation_matrix_np)
        constant = KM.Vector(self.n_dofs, 0.0)

        self.model_part.CreateNewMasterSlaveConstraint(
            "LinearMasterSlaveConstraint", self.constraint_id_start,
            connected_dofs, ref_dofs, relation_matrix, constant)
        self.constraint_id_start += 1

        if self.echo_level > 0:
            KM.Logger.PrintInfo("ApplyRbe3Process",
                f"RBE3 generates: 1 reference ndoe, {n} Connected nodes, "
                f"{self.n_dofs} DOFs/nodes.")

    def _get_weight(self, node):
        if self.weight_var is not None:
            return node.GetSolutionStepValue(self.weight_var)
        return 1.0

    def _build_relation_matrix(self, node, x_r):  
        A = np.eye(self.n_dofs)
    
        if not self.couple_rotation:
            return A
        
        r = np.array([node.X, node.Y, node.Z]) - x_r
        skew_r = self._skew(r)
    
        for t_i, t_name in enumerate(TRANSLATION_NAMES):
            if t_name not in self.idx_map:
                continue
            row = self.idx_map[t_name]
            for r_i, r_name in enumerate(ROTATION_NAMES):
                if r_name not in self.idx_map:
                    continue
                col = self.idx_map[r_name]
                A[row, col] = -skew_r[t_i, r_i] 
        return A
            

    @staticmethod
    def _skew(r):
        """skew matrix, 
        r: vector from connected to reference node
        """
        return np.array([
            [0.0, -r[2], r[1]],
            [r[2], 0.0, -r[0]],
            [-r[1], r[0], 0.0]  
        ])

      