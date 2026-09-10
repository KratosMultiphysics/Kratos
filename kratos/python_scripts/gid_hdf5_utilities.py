import h5py
import hdf5plugin  # Required to read the SZIP compressed datasets safely
import numpy as np

def find_dataset_by_step_and_variable(
    file_path: str, 
    target_step: int, 
    target_variable: str,
    dtype = np.float32
):
    """Scans all groups inside 'Results' to find the specific group where:

    - The group's 'Step' attribute matches 'target_step'
    - The group's 'Name' attribute matches 'target_variable'
    Then returns the data inside it.
    """
    target_step = str(target_step)
    target_variable = str(target_variable)
    try:
        with h5py.File(file_path, "r") as f:
            if "Results" not in f:
                print("Error: 'Results' group not found at the root level.")
                return None

            results_group = f["Results"]

            # Loop through every subgroup inside 'Results'
            for group_name in results_group.keys():
                group = results_group[group_name]

                # Decode the attributes safely
                step_value = group.attrs["Step"].decode("utf-8")
                var_name = group.attrs["Name"].decode("utf-8")

                if step_value == target_step and var_name == target_variable:
                    print(f"Found group: Results/{group_name} containing step_value",step_value," and var_name ",var_name)

                    ids = np.array(group["1"], np.int32)
                    if("4" in group): #vector case
                        data = np.empty((ids.shape[0],3))
                        data[:,0] = np.array(group["2"], dtype=dtype)
                        data[:,1] = np.array(group["3"], dtype=dtype)
                        data[:,2] = np.array(group["4"], dtype=dtype)
                    else: #scalar case
                        data = np.empty((ids.shape[0],1))
                        data[:,0] = group["2"]
                    return ids, data
                       

                    # Find the dataset inside this group and read its contents
                    for child_name in group.keys():
                        child = group[child_name]

                        if isinstance(child, h5py.Dataset):
                            print(
                                f"Reading dataset '{child_name}' from group..."
                            )
                            return (
                                child[:]
                            )  # [:] reads the data array into memory

                    print(
                        f"Warning: Match found, but group 'Results/{group_name}' contains no datasets."
                    )
                    return None

            raise EOFError(f"No match found for Step {target_step} and {target_variable}")
            return None

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


# ==========================================
# Example Usage:
# ==========================================
# path_to_file = "all.h5"
# ids,v = find_dataset_by_step_and_variable(path_to_file, target_step="0.02", target_variable="VELOCITY")
# ids2,p = find_dataset_by_step_and_variable(path_to_file, target_step="0.02", target_variable="PRESSURE")

# print("ids=",ids)
# print("v=",v)
# print("p=",p)

def list_by_variable(file_path: str):
    with h5py.File(file_path, "r", swmr=True) as f:
        results_group = f["Results"]
        var_output = {}

        for group_name in results_group.keys():
            group = results_group[group_name]

            step_value = group.attrs["Step"].decode("utf-8")
            var_name   = group.attrs["Name"].decode("utf-8")

            if var_name not in var_output:
                var_output[var_name] = []

            # Store tuple temporarily
            var_output[var_name].append((step_value, group_name))

        # Now sort and convert to tuple of lists
        for key in var_output:
            # Sort by numeric step_value
            sorted_pairs = sorted(
                var_output[key],
                key=lambda x: float(x[0])
            )

            # Unzip into two parallel lists
            steps, groups = zip(*sorted_pairs)

            # Store as tuple of lists
            var_output[key] = (list(steps), list(groups))

        return var_output
    


##creates a new file from the back
def ExtractFileFromBack(src_path, dst_path, steps_from_back=0):
    with h5py.File(src_path, "r", swmr=True) as src, h5py.File(dst_path, "w") as dst:

        #copy meshes
        #dst.create_group("Meshes")
        src.copy("Meshes", dst)

        # Create the parent group in the destination file
        dst_group = dst.create_group("Results")
        print(dst_group)

        data = list_by_variable(src_path)

        for var in data.keys():
            var_data = data[var]
            print(var_data)
            group_name = var_data[1][-(steps_from_back)]
            print(group_name)

            # Perform a deep copy of everything under the group
            group_to_copy = "Results/"+group_name
            print("group_to_copy-->",group_to_copy)
            src.copy(group_to_copy, dst_group)

            print(f"Copied '{group_to_copy}' to '{dst_path}'")

#aaa = list_by_variable("Pantallas-Pantalla1_interior.h5")
#ExtractFromBack("Pantallas-Pantalla1_interior.h5", "mycopy.h5", steps_from_back=2)