import numpy as np
import KratosMultiphysics as KM

ALL_DOFS = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z",
            "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return ApplyRbe3Process(Model, settings["Parameters"])

class ApplyRbe3Process(KM.Process):
    """\
    @brief 
    Creates an RBE3-type coupling between one target node and multiple source, independent nodes using LinearMasterSlaveConstraint.
    
    @details 
    The selected degrees of freedom of the target node are obtained from a weighted combination of the selected degrees of freedom of the source nodes. 

    The target-node DOFs act as slave DOFs, while the selected DOFs of the source nodes act as master DOFs. 
    Similar to a Nastran RBE3 element, the coupling transfers forces and moments without imposing rigid-body motion between the source nodes.

    Expected parameters: 
    @param    main_model_part_name        :       name of the main modelpart
    @param    target_sub_model_part_name  :       name of the target SubModelPart
    @param    constrained_dofs_target     :       to be constrained DOFs of the target node
    @param    source_groups               :       list of groups containing the source nodes
    @param    source_sub_model_part_name         :       name of the SubModelPart containing the source independent nodes
    @param    source_dofs                        :       DOFs of the source nodes participating in the coupling
    @param    weight                      :       weighing factor assigned to the DOFs of the group
    """

    def __init__(self, model, settings):
        self.echo_level = settings["echo_level"].GetInt() if settings.Has("echo_level") else 0 
        super().__init__()

        default_settings = KM.Parameters("""{
            "main_model_part_name"      :   "",
            "target_sub_model_part_name"  :   "",
            "constraint_id"     :   1,
            "constrained_dofs_target"     :   ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z",
                                            "ROTATION_X","ROTATION_Y","ROTATION_Z"], 
            "source_groups"               :   [{
                "source_sub_model_part_name"  : "",
                "source_dofs"            : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"],
                "weight"          : 1.0
            }]
        }""")
        
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["main_model_part_name"].GetString()]
        self.target_mp = model[settings["target_sub_model_part_name"].GetString()]
        self.constraint_id = settings["constraint_id"].GetInt()
        self.constrained_dofs_target = settings["constrained_dofs_target"].GetStringArray()

        self.target_idx = [ALL_DOFS.index(dof) for dof in self.constrained_dofs_target]

        self.groups = []
        groups = settings["source_groups"]
        group_defaults = default_settings["source_groups"][0]

        for i in range(groups.size()):
            g = groups[i]
            g.ValidateAndAssignDefaults(group_defaults)

            self.groups.append({
                "mp": self.model_part.GetSubModelPart(g["source_sub_model_part_name"].GetString()),
                "source_dofs": g["source_dofs"].GetStringArray(),
                "w": g["weight"].GetDouble()
            })

        if not self.groups:
            raise RuntimeError("RBE3 needs at least 1 group of source nodes.")
        
    def ExecuteInitialize(self):
        """
        Defines the relation between the source nodes and the target node with respect to the weigths.
        """
        if self.target_mp.NumberOfNodes() != 1:
            raise RuntimeError("RBE3 needs exactly 1 target node")
        
        target_node = next(iter(self.target_mp.Nodes))

        if any("ROTATION" in dof for dof in self.constrained_dofs_target):
            print("Rotational DOFs are assigned. Check if the element of the target node has rotational DOFs.")             

        x_r = np.array([target_node.X, target_node.Y, target_node.Z])
        n_target = len(self.constrained_dofs_target)

        A_list = [] 
        w_list = []
        s_dofs = []

        if any("ROTATION" in dof for group in self.groups for dof in group["source_dofs"]):
            print("Rotational DOFs are assigned. Check if the element of the target node has rotational DOFs.")             

        for g in self.groups:
            s_dofs_var=[KM.KratosGlobals.GetVariable(dof) for dof in g["source_dofs"]]
            s_idx=[ALL_DOFS.index(dof) for dof in g["source_dofs"]]

            w = g["w"]  
            if w < 0.0:
                raise RuntimeError("weights must be >= 0")

            for node in g["mp"].Nodes:
                if node.Id == target_node.Id:
                    raise RuntimeError("Target node must not be a source node.")
                
                A_list.append(self._build_relation_matrix(node, x_r, s_idx))
                w_list.append(w)
                s_dofs.extend(node.GetDof(s) for s in s_dofs_var)

        K = sum(w * A.T @ A for A,w in zip(A_list, w_list))

        if np.linalg.matrix_rank(K) < n_target:
            raise RuntimeError("Target dof is not fully determined by source nodes.")

        C_blocks = np.hstack([w * np.linalg.solve(K, A.T) for A, w in zip(A_list, w_list)])
        relation_matrix = KM.Matrix(C_blocks)
        constant = KM.Vector(n_target, 0.0)

        t_dofs_var = [KM.KratosGlobals.GetVariable(dof) for dof in self.constrained_dofs_target]
        t_dofs = [target_node.GetDof(t) for t in t_dofs_var]

        self.model_part.CreateNewMasterSlaveConstraint(
            "LinearMasterSlaveConstraint", self.constraint_id, 
            s_dofs, t_dofs, relation_matrix, constant)  

    def _build_relation_matrix(self, node, x_r, con_idx):  
        """
        Builds the relation matrix between the target and a source node.
        """
        A = np.eye(6)
        r = np.array([node.X, node.Y, node.Z]) - x_r
        A[:3,3:] = -self._skew(r)
        return A[np.ix_(con_idx, self.target_idx)]    

    @staticmethod
    def _skew(r):
        """
        Return the skew matrix, 
        r: vector from the target node to the source node
        """
        return np.array([
            [0.0, -r[2], r[1]],
            [r[2], 0.0, -r[0]],
            [-r[1], r[0], 0.0]  
        ])      