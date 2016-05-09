from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import json
import os

def read_external_json(file_name):
  with open(file_name, 'r') as outfile:
      data = json.load(outfile)
  return data
      
def write_external_json(file_name, data):
  with open(file_name, 'w') as outfile:
    json.dump(data, outfile) 


### TODO: Add the excel utilities
#def import_excel_to_json(input_file, output_file):
    #import xlrd
    
    #wb = xlrd.open_workbook(input_file)

    #sheets = {}

    #for s in wb.sheets():
        #values = {}
        #for col in range(s.ncols):
            #values[str(s.cell(0, col))] = []
            #for row in range(1, s.nrows):
                #values[str(s.cell(0, col))].append(s.cell(row,col).value)
        #sheets[s.name] = values
    
    #write_external_json(output_file, sheets)
    
#def import_json_to_excel(input_file, output_file):
    
    #sheets = read_external_json(input_file)
    #import xlwt
    #workbook = xlrd.open_workbook('input.xls')
    #sheet = workbook.sheet_by_index(0)

    #data = [sheet.cell_value(0, col) for col in range(sheet.ncols)]

    ##workbook = xlwt.Workbook()
    ##sheet = workbook.add_sheet('test')

    ##for index, value in enumerate(data):
        ##sheet.write(0, index, value)

    #workbook.save(output_file)
    