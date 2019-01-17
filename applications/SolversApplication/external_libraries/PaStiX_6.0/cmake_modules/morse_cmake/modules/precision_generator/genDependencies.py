#!/usr/bin/env python
"""@package Tools

This python script is responsible for precision generation replacements
as well as replacements of any kind in other files.  Different types of
replacements can be defined such that no two sets can conflict.  Multiple
types of replacements can, however, be specified for the same file.

@author Wesley Alvaro
@date 2011-4-8

"""
__author__="alvaro"
__date__ ="$Sep 2, 2010 10:09:19 AM$"
__version__=11.0408

import sys;
import re;
import shlex;
import os;
import shutil;
import imp;
from os import path;
from optparse import OptionParser,OptionGroup;
from datetime import datetime;
from Conversion import KEYWORD,DONE_KEYWORD,REGEX,EXTS,Conversion,check_gen,relpath,subs;

class GenConversion:
  """
  This class works on a single file to create generations
  """
  """Static. Source directory where to find the input files"""
  srcdir = None;
  """Static. What (if any) prefix is specified for the output folder?
  If None, use the file's resident folder."""
  prefix = None;
  """Static. Precisions asked by the user."""
  required_precisions = [];
  """The compiled regular expression for detecting proper files."""
  rex = re.compile(REGEX);

  def run(self, file = None):
    match = None;
    """Local complete """
    if self.srcdir is None:
      lfile = file;
    else:
      lfile = path.join(self.srcdir, file);
    """Reads the file and determines if the file needs generation."""
    for line in open(path.realpath(lfile), 'r'):
      m = self.rex.match(line);
      if m is None: continue;
      else:
        match = m.groups();
        break;

    """This file need to go through precision generation script"""
    if match is None:
      return "";
    else:
      try:
        """['normal','all','mixed'] for example. This(ese) are the replacement types to be used."""
        self.types = match[0].split(',');
        """'z' for example. This is the current file's `type`."""
        self.precision = match[2].lower();
        """['c','d','s'] for example. This is the current file's destination `types`."""
        self.precisions = set(match[3].lower().split());
        """Add the local precision to the list if not present"""
        self.precisions.add(self.precision);
        """Take only the intersection of the available precision and the requested ones"""
        self.precisions.intersection_update(self.required_precisions);
      except:
        raise ValueError(lfile+' : Invalid conversion string');

    filename = list(path.split(lfile))[1];
    result = "";
    for precision in self.precisions:
      """For each destination precision, make the appropriate changes to the file name/data."""
      new_file = self.convert(filename, precision);
      if new_file != filename or self.prefix is not None:
        if self.prefix is None:
          """If no prefix is specified, use the file's current folder."""
          prefix = ''
        else:
          """If a prefix is specified, set it up."""
          prefix = self.prefix;
        """Where the destination file will reside."""
        conversion = path.join(prefix, new_file);
        file_out = relpath(conversion);
        result = result+file+","+precision+","+file_out+";";
    return result;

  def substitute(self, sub_type, data, precision):
    """This operates on a single replacement type.
    @param sub_type The name of the replacement set.
    @param data The content subject for replacments.
    @param precision The target precision for replacements.
    """
    try:
      """Try to select the requested replacements."""
      work = subs[sub_type];
      prec_to = work[0].index(precision);
      prec_from = work[0].index(self.precision);
    except:
      """If requested replacement type does not exist,
      return unaltered contents."""
      return data;
    for i in range(1,len(work)):
      """Requested replacements were found,
      execute replacements for each entry."""
      try:
        search = work[i][prec_from];
        replace = work[i][prec_to];
        if not search: continue;
        replace = replace.replace('\*','*');
        if sub_type != 'tracing' :
          replace = replace.replace('\(','(');
          replace = replace.replace('\)',')');
        data = re.sub(search, replace, data);
      except:
        print('Bad replacement pair ', i, 'in', sub_type);
        continue;
    return data;

  def convert(self, data, precision):
    """Select appropriate replacements for the current file.
    @param data The content subject for the replacements.
    @param precision The target precision for generation.
    """
    global KEYWORD, DONE_KEYWORD;
    try:
      """All files undergo the "all" replacements."""
      data = self.substitute('all', data, precision);
    except: pass;
    for sub_type in self.types:
      """For all the other types of conversion for the current file,
      make the correct replacements."""
      if sub_type == 'all': continue;
      try:
        data = self.substitute(sub_type, data, precision);
      except Exception(e):
        raise ValueError('I encountered an unrecoverable error while working in subtype:',sub_type+'.');
    return data;


def main():
  """Create option parser, set static variables of the converter and manage printing options/order."""
  global REGEX, EXTS;
  """The compiled regular expression for detecting proper files."""
  rex = re.compile(REGEX);
  """Files found to be workable."""
  work = [];

  """Create the options parser for detecting options on the command line."""
  parser = OptionParser(usage="Usage: %prog [options]",version='%prog '+str(__version__));
  group = OptionGroup(parser,"Printing Options","These options control the printing output.");
  group.add_option("-i", "--in-files",  help='Print the filenames of files for precision generation.',      action='store_true', dest='in_print',  default=False);
  group.add_option("-o", "--out-files", help='Print the filenames for the precision generated files.',      action='store_true', dest='out_print', default=False);
  group.add_option("-d", "--debug",     help='Print debugging messages.',                                   action='store_true', dest='debug',     default=False);
  parser.add_option_group(group);
  group = OptionGroup(parser,"Settings","These options specify how the work should be done.");
  group.add_option("-s", "--srcdir",    help='The input source directory.',                                 action='store',      dest='srcdir',    default=None);
  group.add_option("-P", "--prefix",    help='The output directory if different from the input directory.', action='store',      dest='prefix',    default=None);
  group.add_option("-f", "--file",      help='Specify a file(s) on which to operate.',                      action='store',      dest='fileslst', type='string', default="");
  group.add_option("-p", "--prec",      help='Specify a precision(s) on which to operate.',                 action='store',      dest='precslst', type='string', default="");
  group.add_option("-e", "--filetypes", help='Specify file extensions on which to operate when walking.',   action='store',      dest='fileexts', type='string', default="");
  group.add_option("-D", "--dictionnary", help='Specify the dictionnary to use in names conversion.',       action='store',      dest='dictionnary', type='string', default="");
  parser.add_option_group(group);

  (options, args) = parser.parse_args();

  """If file extensions are specified, override defaults."""
  if options.fileexts:
    EXTS = options.fileexts.split();

  """Fill the 'work' array with files found to be operable."""
  # if options.fileslst:
  #   """If files to examine are specified on the command line"""
  #   for file in options.fileslst.split():
  #     check_gen(file, work, rex);
  # else:
  #   """Begin directory walking in the current directory."""
  #   startDir = '.';
  #   for root, dirs, files in os.walk(startDir, True, None):
  #     dirs  = filter(hidden,dirs);
  #     files = filter(hidden,files);
  #     files = filter(valid_extension,files);
  #     for file in files:
  #       check_gen(path.join(root,file), work, rex);

  """Set static options for conversion."""
  GenConversion.srcdir = options.srcdir;
  GenConversion.prefix = options.prefix;
  GenConversion.required_precisions = options.precslst.split();

  c = GenConversion(); """This initializes the variable for static member access."""

  result = ""
  for file in options.fileslst.split():
    """For each valid conversion file found."""
    try:
      """Try creating and executing a converter."""
      result += c.run(file);
    except Exception(e):
      print >> sys.stderr, str(e);
      continue;

  print(result);

if __name__ == "__main__":
    main();
