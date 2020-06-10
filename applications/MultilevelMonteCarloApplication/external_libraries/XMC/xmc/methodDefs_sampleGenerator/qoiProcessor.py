def doNothing(rawSolutions):
    #TODO - accept some sort of number QoI from somewhere
    number_qoi = 1

    list_of_qoi = []
    for _ in range(number_qoi):
        qoi_values = []
        for raw_solution in rawSolutions:
            qoi_values.append(raw_solution)
        list_of_qoi.append(qoi_values)

    return list_of_qoi
