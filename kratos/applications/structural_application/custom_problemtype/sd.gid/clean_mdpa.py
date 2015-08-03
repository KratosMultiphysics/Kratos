import sys

# for arg in sys.argv:
#    print arg

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

reading_element = False

elem_count = 0
for line in infile.readlines():
    words = line.split()
    if len(words) > 1:
        if words[0] == "Begin":
            if words[1] == "Elements" or words[1] == "Conditions" or words[1] == "NodalData" or words[1] == "ElementalData":
                reading_element = True
                elem_count = 0
                text = []
                text.append(line)
                continue

        if words[0] == "End":
            if words[1] == "Elements" or words[1] == "Conditions" or words[1] == "NodalData" or words[1] == "ElementalData":
                reading_element = False
                text.append(line)
                if elem_count > 0:
                    outfile.write('\n')
                    for t in text:
                        outfile.write(t)
                continue

    if reading_element and len(words) > 0:
        elem_count = elem_count + 1
        text.append(line)
    else:
        if(len(words) > 0):
            if(words[0] == "Begin"):
                outfile.write('\n')
            if(not(words[0] == '\n')):
                outfile.write(line)

infile.close()
outfile.close()

