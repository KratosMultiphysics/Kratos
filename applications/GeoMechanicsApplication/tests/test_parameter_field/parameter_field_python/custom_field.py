
def run(input,output):

    new_values = []
    for value, coord in zip(input["values"], input["coordinates"]):

        new_value = value*2*coord[0] + value*3*coord[1]
        new_values.append(new_value)

    output["values"] = new_values


