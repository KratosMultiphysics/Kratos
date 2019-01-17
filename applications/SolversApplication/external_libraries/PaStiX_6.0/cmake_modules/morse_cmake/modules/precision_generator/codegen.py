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
from os import path;
from optparse import OptionParser,OptionGroup;
from datetime import datetime;
from Conversion import KEYWORD,DONE_KEYWORD,REGEX,EXTS,Conversion,check_gen,subs;

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
  group.add_option("-m", "--make",      help='Spew a GNU Make friendly file to standard out.',              action='store_true', dest='make',      default=False);
  group.add_option("-d", "--debug",     help='Print debugging messages.',                                   action='store_true', dest='debug',     default=False);
  parser.add_option_group(group);
  group = OptionGroup(parser,"Operating Mode Options","These options alter the way the program operates on the input/output files.");
  group.add_option("-c", "--clean",     help='Remove the files that are the product of generation.',        action='store_true', dest='out_clean', default=False);
  group.add_option("-T", "--test",      help='Don\'t actually do any work.',                                action='store_true', dest='test',      default=False);
  parser.add_option_group(group);
  group = OptionGroup(parser,"Settings","These options specify how the work should be done.");
  group.add_option("-P", "--prefix",    help='The output directory if different from the input directory.', action='store',      dest='prefix',    default=None);
  group.add_option("-f", "--file",      help='Specify a file(s) on which to operate.',                      action='store',      dest='fileslst', type='string', default="");
  group.add_option("-p", "--prec",      help='Specify a precision(s) on which to operate.',                 action='store',      dest='precslst', type='string', default="");
  group.add_option("-e", "--filetypes", help='Specify file extensions on which to operate when walking.',   action='store',      dest='fileexts', type='string', default="");
  parser.add_option_group(group);

  (options, args) = parser.parse_args();

  """If file extensions are specified, override defaults."""
  if options.fileexts:
    EXTS = options.fileexts.split();

  """Fill the 'work' array with files found to be operable."""
  if options.fileslst:
    """If files to examine are specified on the command line"""
    for file in options.fileslst.split():
      check_gen(file, work, rex);
  else:
    """Begin directory walking in the current directory."""
    startDir = '.';
    for root, dirs, files in os.walk(startDir, True, None):
      dirs  = filter(hidden,dirs);
      files = filter(hidden,files);
      files = filter(valid_extension,files);
      for file in files:
        check_gen(path.join(root,file), work, rex);

  """Set static options for conversion."""
  Conversion.debug  = options.debug;
  Conversion.make   = options.make;
  Conversion.prefix = options.prefix;
  Conversion.required_precisions = options.precslst.split();
  if options.out_print or options.out_clean or options.in_print or options.make or options.test:
    Conversion.test = True;

  if options.make:
    """If the program should be GNU Make friendly."""
    print('## Automatically generated Makefile');
    print('PYTHON ?= python');

  c = Conversion(); """This initializes the variable for static member access."""

  for tuple in work:
    """For each valid conversion file found."""
    try:
      """Try creating and executing a converter."""
      c = Conversion(tuple[0], tuple[1], tuple[2]);
      c.run();
    except Exception(e):
      print >> sys.stderr, str(e);
      continue;

  if options.make:
    """If the program should be GNU Make friendly."""
    print('gen = ',' '+' '.join(c.files_out));
    print('cleangen:');
    print('\trm -f $(gen)');
    print('generate: $(gen)');
    print('.PHONY: cleangen generate');
  if options.in_print:
    """Should we print the input files?"""
    print(' '.join(c.files_in));
  if options.out_print:
    """Should we print the output files?"""
    print(' '.join(c.files_out));
  if options.out_clean:
    """Clean generated files"""
    for file in c.files_out:
      if not path.exists(file): continue;
      os.remove(file);

if __name__ == "__main__":
    main();
