import pickle
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict

model = Kratos.Model()
model_part = model.CreateModelPart("test")

ta1 = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes, Kratos.PRESSURE)
ta2 = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes, Kratos.VELOCITY)

buff_dict = BufferedDict(2)

buff_dict.SetValue("test/test1/ta1", ta1)
buff_dict.SetValue("test/test1/ta2", ta2)

with open("whatever.pkl", "wb") as file_output:
    serializer = Kratos.StreamSerializer()
    serializer.Save("SerializedModel", model)

    # serialize the buffered dict
    data_entry_type: list[tuple[str, str]] = []
    dict_of_str_values = buff_dict.GetMap()
    for k, v in dict_of_str_values.items():
        data_entry_type.append((k, type(v)))
        serializer.Save(f"buff_dict:{k}", v)
    serializer.Save("dict_entries", data_entry_type)

    pickle.dumps(serializer)



with open("watever", "rb") as file_inp:
    serializer: Kratos.StreamSerializer = pickle.load(file_inp)
    serializer.Load("SerializedModel", model)

    data_entry_type: list[tuple[str, str]] = []
    serializer.Load("dict_entries", data_entry_type)

    for k, v_t in data_entry_type:
        if v_t == Kratos.TensorAdaptors.VariableTensorAdaptor:
            ta = Kratos.TensorAdaptors.VariableTensorAdaptor()
        serializer.Load(f"buff_dict:{k}", ta)
        buff_dict.SetValue(k, ta)




