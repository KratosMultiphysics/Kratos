import re
import ctypes


def clean_c_buffer(char_buffer):
    try:
        raw = char_buffer.raw.split(b'\x00')[0]
        decoded = raw.decode("utf-8", errors="ignore")
        decoded = re.sub(r"^\s*\d+\s+", "", decoded)
        decoded_clean = re.sub(r"[^\w\s\-\(\)\[\]\.,':=+/]", "", decoded)
        return decoded_clean.strip()
    except Exception:
        return "<?>"

def find_symbol_in_dll(dll_path, dll_lib, symbol_name):
    import pefile
    try:
        pe = pefile.PE(dll_path)
        symbol_name_lower = symbol_name.lower()

        for exp in pe.DIRECTORY_ENTRY_EXPORT.symbols:
            name = exp.name
            if name and name.decode().lower() == symbol_name_lower:
                return getattr(dll_lib, name.decode())

        return None
    except Exception as e:
        print(f"Error reading DLL: {e}")
        return None

def get_model_count(getmodelcount):
    getmodelcount.argtypes = (ctypes.POINTER(ctypes.c_int),)
    getmodelcount.restype = None
    result = ctypes.c_int()
    getmodelcount(ctypes.byref(result))
    return result.value

def get_model_name(getmodelname, model_no):
    BUFFER_SIZE = 256
    char_buffer = ctypes.create_string_buffer(BUFFER_SIZE)
    getmodelname.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, ctypes.c_long)
    getmodelname.restype = None
    getmodelname(ctypes.byref(ctypes.c_int(model_no)), char_buffer, ctypes.c_long(BUFFER_SIZE))
    return clean_c_buffer(char_buffer)

def get_param_count(getparamcount, model_no):
    getparamcount.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
    getparamcount.restype = None
    result = ctypes.c_int()
    getparamcount(ctypes.byref(ctypes.c_int(model_no)), ctypes.byref(result))
    return result.value

def get_param_name(getparamname, model_no, param_no):
    BUFFER_SIZE = 256
    char_buffer = ctypes.create_string_buffer(BUFFER_SIZE)
    getparamname.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, ctypes.c_long)
    getparamname.restype = None
    getparamname(ctypes.byref(ctypes.c_int(model_no)), ctypes.byref(ctypes.c_int(param_no)), char_buffer, ctypes.c_long(BUFFER_SIZE))
    return clean_c_buffer(char_buffer)

def udsm_parser(dll_path):
    dll_lib = ctypes.CDLL(dll_path, winmode=0)

    getmodelcount = find_symbol_in_dll(dll_path, dll_lib, "getmodelcount")
    getmodelname = find_symbol_in_dll(dll_path, dll_lib, "getmodelname")
    getparamcount = find_symbol_in_dll(dll_path, dll_lib, "getparamcount")
    getparamname = find_symbol_in_dll(dll_path, dll_lib, "getparamname")

    if not all([getmodelcount, getmodelname, getparamcount, getparamname]):
        raise ValueError("One or more required symbols not found in the DLL.")

    model_count = get_model_count(getmodelcount)

    model_dict = {"model_name": [], "num_params": [], "param_names": []}
    for model_no in range(1, model_count + 1):
        model_dict["model_name"].append(get_model_name(getmodelname, model_no))
        num_params = get_param_count(getparamcount, model_no)
        model_dict["num_params"].append(num_params)
        param_names = [get_param_name(getparamname, model_no, p) for p in range(1, num_params + 1)]
        model_dict["param_names"].append(param_names)

    return model_dict
