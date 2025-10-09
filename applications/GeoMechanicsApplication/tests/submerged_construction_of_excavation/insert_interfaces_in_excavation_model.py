from interface_inserter import InterfaceInserter
from mdpa_reader import MdpaReader
from mdpa_writer import MdpaWriter

def main():
    input_file = "submerged_excavation_gid_project.mdpa"
    lines = []
    with open(input_file, "r") as file:
        lines = file.readlines()
    mesh = MdpaReader().read_mesh(lines)

    interface_inserter_left = InterfaceInserter()
    interface_inserter_left.insert_interface(
        mesh, "Clay_Left", "Diaphragm_Wall", "Wall_Clay_Interface_Left"
    )
    interface_inserter_left.insert_interface(
        mesh, "Sand_Left", "Diaphragm_Wall", "Wall_Sand_Interface_Left"
    )

    interface_inserter_right = InterfaceInserter()
    interface_inserter_right.insert_interface(
        mesh, "Clay_Upper_Right", "Diaphragm_Wall", "Wall_Clay_Interface_Upper_Right"
    )
    interface_inserter_right.insert_interface(
        mesh, "Clay_Middle_Right", "Diaphragm_Wall", "Wall_Clay_Interface_Middle_Right"
    )
    interface_inserter_right.insert_interface(
        mesh, "Clay_Lower_Right", "Diaphragm_Wall", "Wall_Clay_Interface_Lower_Right"
    )
    interface_inserter_right.insert_interface(
        mesh, "Sand_Right", "Diaphragm_Wall", "Wall_Sand_Interface_Right"
    )

    output_lines = MdpaWriter().write_mesh(mesh)
    with open("submerged_excavation_gid_project_with_interfaces.mdpa", "w") as file:
        file.writelines(line + "\n" for line in output_lines)


if __name__ == "__main__":
    main()
