#!/usr/bin/env python3

import os, sys, shutil, copy, glob
from SU2 import run
from SU2.util import ordered_bunch 
from SU2.util.ordered_dict import OrderedDict

inf = 1.0e20

# ----------------------------------------------------------------------
#  Configuration Class
# ----------------------------------------------------------------------

class genSampling(ordered_bunch):

    _konfig = {}
    
    def __init__(self,*args,**kwarg):

        # initialize ordered bunch
        super(genSampling,self).__init__(*args,**kwarg)

    def read(self):
        """ reads from a config file """
        if self.cfgOpts:
            try:
                cfg = 'fsi'
                self._konfig['fsi'] = read_config(self.cfgOpts['cfgFiles']['fsi'], self.cfgOpts['cfgPath'])
                cfg = 'fluid'
                self._konfig['fluid'] = read_config(self.cfgOpts['cfgFiles']['fluid'], self.cfgOpts['cfgPath'])
                cfg = 'solid'
                self._konfig['solid'] = read_config(self.cfgOpts['cfgFiles']['solid'], self.cfgOpts['cfgPath'])
            except IOError:
                print('Could not find config file for: %s' % cfg)
            except:
                print('Unexpected error: ', sys.exc_info()[0])
                raise 

    def write(self):
        assert os.path.exists(self.pathDOE['train']) , 'Generate training DOE first'
        assert os.path.exists(self.pathDOE['valid']) , 'Generate validation DOE first'
        
        cfgPathTrain = os.path.join(self.cfgOpts['cfgPath'], 'train')
        if not os.path.isdir(cfgPathTrain):
            os.mkdir(cfgPathTrain)
        if not os.path.isdir(self.outDir['train']):
            os.makedirs(self.outDir['train'])
        write_config(self.pathDOE['train'], self.cfgOpts, self._konfig, cfgPathTrain, self.outDir['train'])
        
        cfgPathValid = os.path.join(self.cfgOpts['cfgPath'], 'valid')
        if not os.path.isdir(cfgPathValid):
            os.mkdir(cfgPathValid)
        if not os.path.isdir(self.outDir['valid']):
            os.makedirs(self.outDir['valid'])
        write_config(self.pathDOE['valid'], self.cfgOpts, self._konfig, cfgPathValid, self.outDir['valid'])

    def run(self):
        training = glob.glob(self.cfgOpts['cfgPath'] + '/train/config_fsi_steady_*')
        for cfg in training:
            run.run_command("mpirun -n 8 -v SU2_CFD {}".format(cfg))
        
        validation = glob.glob(self.cfgOpts['cfgPath'] + '/valid/config_fsi_steady_*')
        for cfg in validation:
            run.run_command("mpirun -n 8 SU2_CFD {}".format(cfg))
#
    def __eq__(self,k):
        return super(genSampling,self).__eq__(k)
    def __ne__(self,k):
        return super(genSampling,self).__ne__(k)
    def __getattr__(self,k):
        try:
            return super(genSampling,self).__getattr__(k)
        except AttributeError:
            raise AttributeError('Config parameter not found')

    def __getitem__(self,k):
        try:
            return super(genSampling,self).__getitem__(k)
        except KeyError:
            raise KeyError('Config parameter not found: %s' % k)


    def local_files(self):
        """ removes path prefix from all *_FILENAME params
        """
        for key, value in self.items():
            if key.split('_')[-1] == 'FILENAME':
                self[key] = os.path.basename(value)

    def __repr__(self):
        #return '<Config> %s' % self._filename
        return self.__str__()

    def __str__(self):
        output = 'Config: %s' % self._filename
        for k,v in self.items():
            output +=  '\n    %s= %s' % (k,v)
        return output
#: class Config

# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_config(filename, configPath):
    """ reads a config file """

    # initialize output dictionary
    data_dict = {}

    input_file = open(configPath + '/' + filename)

    # process each line
    while 1:
        # read the line
        line = input_file.readline()
        if not line:
            break

        # remove line returns
        line = line.strip('\r\n').strip()

        if (len(line) == 0):
            continue
        # make sure it has useful data
        if (line[0] == '%'):
            continue

        # --- Check if there is a line continuation character at the
        # end of the current line or somewhere in between (the rest is ignored then).
        # If yes, read until there is a line without one or an empty line.
        # If there is a statement after a cont. char
        # throw an error. ---*/

        while(line[0].endswith('\\') or len(line.split('\\')) > 1):
            tmp_line = input_file.readline()
            tmp_line = tmp_line.strip()
            assert len(tmp_line.split('=')) <= 1, ('Statement found after line '
                                                   'continuation character in config file %s' % tmp_line)
            if (not tmp_line.startswith('%')):
                line = line.split('\\')[0]
                line += ' ' + tmp_line

        # split across equals sign
        line = line.split("=",1)
        this_param = line[0].strip()
        this_value = line[1].strip()

        assert this_param not in data_dict, ('Config file has multiple specifications of %s' % this_param )
        data_dict[this_param] = this_value

    return data_dict

#: def read_config()



# -------------------------------------------------------------------
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def write_config(doePath,cfgOpts, cfgDict, path, outDir):

    fsiFileName = cfgOpts['cfgFiles']['fsi']
    fsiParamDict = copy.deepcopy(cfgDict['fsi'])
    parVarFileName = cfgOpts['cfgFiles'][cfgOpts['parVar']['cfgFileType']]
    parVarDict = copy.deepcopy(cfgDict[cfgOpts['parVar']['cfgFileType']])
    parVarName = cfgOpts['parVar']['name']

    if (cfgOpts['parVar']['cfgFileType'] == 'fluid'):
        constFileName = cfgOpts['cfgFiles']['solid']
    else:
        constFileName = cfgOpts['cfgFiles']['fluid']

    shutil.copy(cfgOpts['cfgPath'] + '/' + constFileName, path + '/' + constFileName)
    
    doeList = open(doePath + '/listDOE.txt','r')
    i = 0
    while 1:
        # read the line
        lineDOE = doeList.readline()
        if not lineDOE:
            break

        # remove line returns
        lineDOE = lineDOE.strip('\r\n').strip()

        if (len(lineDOE) == 0):
            continue
        # make sure it has useful data
        if (lineDOE[0] == '%'):
            continue
        
        
        output_file_fsi = open(path + '/' + fsiFileName.rstrip('.cfg') + '_{}'.format(i) + '.cfg' ,"w")
        for key,value in fsiParamDict.items():

            if (key == 'CONFIG_LIST'):
                new_value = '(' + path + '/' + parVarFileName.rstrip('.cfg') + '_{}'.format(i) + '.cfg' + ',' + \
                           path + '/' + constFileName + ')'
                output_file_fsi.write(key + " = " + new_value)
            elif (key == 'RESTART_FILENAME'):
                new_value = outDir + '/' + 'restart_{}'.format(i)
                output_file_fsi.write(key + " = " + new_value)
            elif (key == 'VOLUME_FILENAME'):
                new_value = outDir + '/' + 'volume_{}'.format(i)
                output_file_fsi.write(key + " = " + new_value)
            else:
                output_file_fsi.write(key + " = " + value)

            # next line
            output_file_fsi.write("\n")

        output_file_Var = open(path + '/' + parVarFileName.rstrip('.cfg') + '_{}'.format(i) + '.cfg' ,"w")
        for key,value in parVarDict.items():

            if (key == parVarName):
                output_file_Var.write(key + " = " + lineDOE.removesuffix('_0.dat'))
            else:
                output_file_Var.write(key + " = " + value)
            output_file_Var.write("\n")

        i += 1

        output_file_fsi.close()
#        output_file_Var.close()

#: def write_config()


