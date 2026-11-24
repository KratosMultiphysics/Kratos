import h5py
import hdf5plugin  # Required to read the SZIP compressed datasets safely
import numpy as np

def _decode_attribute(value):
    """Decode an HDF5 attribute that may be stored either as bytes or as str."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)

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

    Returns (ids, data) on success, None if the file or the 'Results' group is
    missing (or on any other read error), and raises EOFError if no group matches.
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

                # Decode the attributes safely (h5py may return bytes or str depending on version)
                step_value = _decode_attribute(group.attrs["Step"])
                var_name = _decode_attribute(group.attrs["Name"])

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

            raise EOFError(f"No match found for Step {target_step} and {target_variable}")

    except EOFError:
        # 'Not found' is an exceptional (but well-defined) outcome: propagate it to the caller
        raise
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

            step_value = _decode_attribute(group.attrs["Step"])
            var_name   = _decode_attribute(group.attrs["Name"])

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
            # Select the group counting from the back: steps_from_back=0 means the most recent
            # (note that "-0" would silently index the first group, hence the explicit -1).
            # If a variable has fewer steps than requested, clamp to its oldest available step.
            steps_back = min(steps_from_back, len(var_data[1]))
            group_name = var_data[1][-steps_back if steps_back > 0 else -1]
            print(group_name)

            # Perform a deep copy of everything under the group
            group_to_copy = "Results/"+group_name
            print("group_to_copy-->",group_to_copy)
            src.copy(group_to_copy, dst_group)

            print(f"Copied '{group_to_copy}' to '{dst_path}'")

#aaa = list_by_variable("Pantallas-Pantalla1_interior.h5")
#ExtractFromBack("Pantallas-Pantalla1_interior.h5", "mycopy.h5", steps_from_back=2)