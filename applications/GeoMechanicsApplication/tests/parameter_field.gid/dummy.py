
def run(input,output):

    i=1
    for k,v in input.items():

        output[k] = v["value"]*i
        i+=1

    print(input)
    print(output)

