import json

def read_external_json(file_name):
  with open(file_name, 'r') as outfile:
      data = json.load(outfile)
  return data

def write_external_json(file_name, data):
  with open(file_name, 'w') as outfile:
    json.dump(data, outfile, indent=2)