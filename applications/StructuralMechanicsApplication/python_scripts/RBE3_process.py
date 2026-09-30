
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

ALL_DOFS = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z",
            "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]

KM_ALL_DOFS = [KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z,
            KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z]

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return ApplyRbe3Process(Model, settings["Parameters"])

class ApplyRbe3Process(KM.Process):
    """\
    @brief 
    Creates an RBE3-type coupling between one reference node and multiple connected, independent nodes using LinearMasterSlaveConstraint.
    
    @details 
    The selected degrees of freedom of the reference node are obtained from a weighted combination of the selected degrees of freedom of the connected nodes. 

    The reference-node DOFs act as slave DOFs, while the selected DOFs of the connected nodes act as master DOFs. 
    Similar to a Nastran RBE3 element, the coupling transfers forces and moments without imposing rigid-body motion between the connected nodes.

    Expected parameters: 
        model_part_name             :       name of the modelpart
        reference_sub_model_part    :       name of the reference SubModelPart
        constraint_id_start         :       Constraint-ID
        constrained_dofs_ref        :       to be constrained DOFs of the reference node
        connected_groups            :       list of connected nodes including
                sub_model_part          :       name of the SubModelPart containing the connected independent nodes
                dofs                    :       DOFs of the connected nodes participating in the coupling
                weight                  :       weighing factor assigned to the DOFs of the group
    """
    def __init__(self, model, settings):
        self.echo_level = settings["echo_level"].GetInt() if settings.Has("echo_level") else 0 
        super().__init__()

        group_default_settings = KM.Parameters("""{
            "sub_model_part"  : "",
            "dofs"            : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "weight"          : 1.0
        }""")

        default_settings = KM.Parameters("""{
            "model_part_name"           :   "",
            "reference_sub_model_part"  :   "",
            "constraint_id_start"       :   1, 
            "constrained_dofs_ref"      :   ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z",
                                            "ROTATION_X","ROTATION_Y","ROTATION_Z"], 
            "connected_groups"          :   []
        }""")
        
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.reference_mp = model[settings["reference_sub_model_part"].GetString()]
        self.constraint_id_start = settings["constraint_id_start"].GetInt()
        self.ref_idx = self._to_idx(settings["constrained_dofs_ref"].GetStringArray())

        self.groups = []
        groups = settings["connected_groups"]
        for i in range(groups.size()):
            groups[i].ValidateAndAssignDefaults(group_default_settings)
            g = groups[i]
            self.groups.append({
                "mp": self.model_part.GetSubModelPart(
                    g["sub_model_part"].GetString()),
                "idx": self._to_idx(g["dofs"].GetStringArray()),
                "w": g["weight"].GetDouble()
            })
        if not self.groups:
            raise RuntimeError("RBE3 needs at least 1 group of connected nodes.")
 
        
    def ExecuteInitialize(self):
        """
        Defines the relation between the connected nodes and the reference node with respect to the weigths.
        """

        if self.reference_mp.NumberOfNodes() != 1:
            raise RuntimeError("RBE3 needs exactly 1 reference node")
        
        ref_node = next(iter(self.reference_mp.Nodes))
        for i in self.ref_idx:
            if i >= 3 and not ref_node.HasDofFor(KM_ALL_DOFS[i]):
                raise RuntimeError("Element has no rotational DOF.") #for example a solid element

        x_r = np.array([ref_node.X, ref_node.Y, ref_node.Z])
        n_ref = len(self.ref_idx)

        A_list = [] 
        w_list = []
        reference_dofs = []

        for g in self.groups:
            for node in g["mp"].Nodes:
                if node.Id == ref_node.Id:
                    raise RuntimeError("Reference node must not be a connected node.")
                w = g["w"]     

                if w < 0.0:
                    raise RuntimeError("weights must be >= 0")

                A_list.append(self._build_relation_matrix(node, x_r, g["idx"]))
                w_list.append(w)
                reference_dofs.extend(node.GetDof(KM_ALL_DOFS[i]) for i in g["idx"])

                for i in g["idx"]:
                    if i >= 3 and not node.HasDofFor(KM_ALL_DOFS[i]):
                        raise RuntimeError("Element has no rotational DOF.")

        K = sum(w * A.T @ A for A,w in zip(A_list, w_list))

        if np.linalg.matrix_rank(K) < n_ref:
            raise RuntimeError("Reference dof is not fully determined by connected nodes.")

        C_blocks = np.hstack([w * np.linalg.solve(K, A.T) for A, w in zip(A_list, w_list)])
        relation_matrix = KM.Matrix(C_blocks)

        constant = KM.Vector(n_ref, 0.0)
        connected_dofs = [ref_node.GetDof(KM_ALL_DOFS[i]) for i in self.ref_idx]

        self.model_part.CreateNewMasterSlaveConstraint(
            "LinearMasterSlaveConstraint", self.constraint_id_start,
            reference_dofs, connected_dofs, relation_matrix, constant)
        self.constraint_id_start += 1

        if self.echo_level > 0:
            KM.Logger.PrintInfo("ApplyRbe3Process",
                f"RBE3 generates: 1 reference ndoe, {n} Connected nodes, "
                f"{self.n_dofs} DOFs/nodes.")
            

    def _build_relation_matrix(self, node, x_r, con_idx):  
        """
        Builds the relation matrix between the reference and a connected node.
        """
        A = np.eye(6)
        r = np.array([node.X, node.Y, node.Z]) - x_r
        A[:3,3:] = -self._skew(r)
        return A[np.ix_(con_idx, self.ref_idx)] 
      

    @staticmethod
    def _skew(r):
        """
        Return the skew matrix, 
        r: vector from the reference node to the connected node
        """
        return np.array([
            [0.0, -r[2], r[1]],
            [r[2], 0.0, -r[0]],
            [-r[1], r[0], 0.0]  
        ])

    @staticmethod
    def _to_idx(names):
        """
        Assigns the numbers 0-5 to the DOFs.
        """
        invalid = set(names) - set(ALL_DOFS)
        if invalid:
            raise RuntimeError(f"RBE3: invalid DOFs {invalid}")
        return [ALL_DOFS.index(n) for n in names]

      