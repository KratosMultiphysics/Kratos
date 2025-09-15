import os
import glob

def find_and_replace_ref_values(file_name, ref_text_lines):
    print("Reading file " + file_name + "...")
    with open(file_name, "r") as file_input:
        lines = file_input.readlines()

    ref_line_index = 0
    while (ref_line_index < len(ref_text_lines)):
        ref_line = ref_text_lines[ref_line_index]
        s_ref_line = ref_line.strip()
        if (s_ref_line.startswith("Test")):
            ref_test_name = s_ref_line[4:s_ref_line.find(" ")]

            # get ref values
            ref_line_index+=3
            data = []
            while (ref_line_index < len(ref_text_lines)):
                ref_line_data = ref_text_lines[ref_line_index].strip()
                if (ref_line_data == ""):
                    break
                else:
                    name = ref_line_data[:ref_line_data.find("=")].strip()
                    value = ref_line_data[ref_line_data.find("=")+1:]
                    data.append([name, value])
                ref_line_index +=1

            line_index = 0
            while (line_index < len(lines)):
                line = lines[line_index]
                if (line.startswith("KRATOS_TEST_CASE_IN_SUITE(")):
                    test_name = line[line.find("(") + 1:line.find(",")]
                    line_index += 1
                    if (ref_test_name == test_name):
                        while(line_index < len(lines)):
                            data_line = lines[line_index]
                            if (data_line.startswith("}")):
                                break
                            else:
                                for ref_data in data:
                                    if data_line.strip().startswith("ref_" + ref_data[0]):
                                        lines[line_index] = data_line[:data_line.find("=")] + "= " + ref_data[1] + ";\n"
                                        break
                            line_index += 1

                line_index +=1

        ref_line_index += 1

    with open(file_name, "w") as file_output:
        file_output.writelines(lines)



if __name__=="__main__":
    path = "/home/suneth/software/kratos_temp_3/applications/RANSApplication/tests/"
    ref_file_name = "/home/suneth/software/kratos_temp_3/applications/RANSApplication/tests/log"

    with open(ref_file_name, "r") as file_input:
        lines = file_input.readlines()

    files = glob.glob(path + '/**/*.cpp', recursive=True)
    for item in files:
        find_and_replace_ref_values(item, lines)
