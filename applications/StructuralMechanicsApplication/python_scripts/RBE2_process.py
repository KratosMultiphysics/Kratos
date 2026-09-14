
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
    Erzeugt eine RBE2-Starrkoerperkopplung (LinearMasterSlaveConstraint) zwischen einem Master-Knoten und einem SubModelPart voller
    Slave-Knoten.
    Erwartete Parameter:
        model_part_name   : Name der ModelPart, die master/slave enthaelt
        master_sub_model_part : SubModelPart-Name mit genau 1 Knoten (Master)
        slave_sub_model_part  : SubModelPart-Name mit den Slave-Knoten
        constraint_id_start   : erste freie MasterSlaveConstraint-Id
        constrained dofs: hier kann man wählen, was alles übetragen werden soll
    """
    def __init__(self, model, settings):
        self.echo_level = settings["echo_level"].GetInt() if settings.Has("echo_level") else 0 #Ausgaben
        super().__init__() #ruft den Konstruktor der Basisklasse KM.Process auf, 
        # zuerst wird Elternklasse initialisiert, dann die eigene Klasse
        default_settings = KM.Parameters("""{
            "model_part_name"       : "",
            "master_sub_model_part" : "",
            "slave_sub_model_part"  : "",
            "constraint_id_start"   : 1, 
            "constrained_dofs"      : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z",
                                         "ROTATION_X","ROTATION_Y","ROTATION_Z"]
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

# self --> instance of the class
        self.model_part = model[settings["model_part_name"].GetString()]
        self.master_mp = model[settings["master_sub_model_part"].GetString()]
        self.slave_mp = model[settings["slave_sub_model_part"].GetString()]
        self.constraint_id_start = settings["constraint_id_start"].GetInt()

        requested = settings["constrained_dofs"].GetStringArray()
        invalid = set(requested) - set(DOF_MAP) # wenn in Json ein DOF-Name steht, der nicht in DOF_MAP ist, dann ist invalid = True
        if invalid:
            raise Exception(f"Invalid input: {invalid}")
        
        self.dof_names = requested #string
        self.dof_vars = [DOF_MAP[name] for name in self.dof_names] #Kratos-Variable (z.B. KM.DISPLACEMENT_X) 
        self.n_dofs = len(self.dof_vars)
        self.idx_map = {name: i for i, name in enumerate(self.dof_names)} #DOFs bekommen Nummern 0,1,2..

        # Skew-Kopplung nur einbauen, wenn Translation UND Rotation gewählt sind
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
        T = np.eye(self.n_dofs)

        if not self.couple_rotation:
            return T
        r = np.array([slave_node.X, slave_node.Y, slave_node.Z]) - x_m
        skew_r = self._skew(r)
        #KM.Logger.PrintInfo(r)

        for t_i, t_name in enumerate(TRANSLATION_NAMES):
            if t_name not in self.idx_map:
                continue
            row = self.idx_map[t_name]
            for r_i, r_name in enumerate(ROTATION_NAMES):
                if r_name not in self.idx_map:
                    continue
                col = self.idx_map[r_name]
                T[row, col] = -skew_r[t_i, r_i]
        #KM.Logger.PrintInfo(T)
        return T
          
# hier wird die skew matrix definiert; sie ist Teil der 6x6 Matrix und stellt das kreuzprodukt r × θ dar, das die Rotationsanteile der Slave-Knoten beschreibt.
# man kann das nicht über ein Kreuprodukt cross machen, da es noch in die 6x6 matrix eingebettet werden muss
    @staticmethod
    def _skew(r):
        return np.array([
            [0.0, -r[2], r[1]],
            [r[2], 0.0, -r[0]],
            [-r[1], r[0], 0.0]  #mit theta multipliziert: (-ry)·θx + rx·θy + 0·θz  =  rx·θy - ry·θx
        ])
    #jede Zeile drückt eine Komponente (x, y, z) des Kreuzprodukts r × θ als lineare Kombination von θx, θy, θz aus
