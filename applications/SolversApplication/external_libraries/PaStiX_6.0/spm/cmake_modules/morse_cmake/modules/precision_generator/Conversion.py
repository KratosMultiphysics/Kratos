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
from subs import subs;
from datetime import datetime;

"""Keyword used to signal replacement actions on a file"""
KEYWORD = '@precisions';
"""Replace the above keyword with this one post-replacements"""
DONE_KEYWORD = '@generated';
"""Regular expression for the replacement formatting"""
REGEX = '^.*'+KEYWORD+'\s+((\w+,?)+)\s+(\w+)\s+->\s*((\s\w+)+).*$';
"""Default acceptable extensions for files during directory walking"""
EXTS = ['.c','.cpp','.h','.hpp','.f','.jdf','.f90','.F90','.f77','.F77','.cu','.cuf','.CUF','.jdf'];

def check_gen(file, work, rex):
  """Reads the file and determines if the file needs generation."""
  fd = open(path.realpath(file), 'r');
  lines = fd.readlines();
  fd.close();
  for line in lines:
    m = rex.match(line);
    if m is None: continue;
    work.append((file, m.groups(), ''.join(lines)));

def grep(string,list):
  expr = re.compile(string)
  return filter(expr.search,list)

def hidden(file):
  """Exclude hidden files"""
  return not file.startswith('.');

def valid_extension(file):
  """Exclude non-valid extensions"""
  global EXTS;
  for ext in EXTS:
    if file.endswith(ext):
      return True;
  return False;

def relpath(p):
  """Get the relative path of a file."""
  p = path.realpath(p);
  return p.replace(path.realpath('.')+'/','');

class Conversion:
  """
  This class works on a single file to create generations
  """
  """Static. Is the conversion in debug mode? More verbose."""
  debug = False;
  """Static. Is the conversion in test mode? No real work."""
  test = False;
  """Static. Is the conversion in make mode? Output make commands."""
  make = False;
  """Static. What (if any) prefix is specified for the output folder?
  If None, use the file's resident folder."""
  prefix = None;
  required_precisions = [];
  """Static. A running list of files that are input."""
  files_in = [];
  """Static. A running list of files that are output."""
  files_out = [];
  """Static. A running list of files that are output."""
  dependencies = [];
  def __init__(self, file = None, match = None, content = None):
    """Constructor that takes a file, match, and content.
    @param file The file name of the input.
    @param match The regular expression matches
    @param content The ASCII content of the file.
    """
    if file is None: return;
    self.content = content;
    #file = path.realpath(file);
    rel = relpath(file);
    self.file = list(path.split(file));
    self.date = path.getmtime(file);
    if sys.platform!="win32" and path.samefile(path.join(self.file[0],self.file[1]),sys.argv[0]):
      raise ValueError('Let\'s just forget codegen.py');
    try:
      """['normal','all','mixed'] for example. This(ese) are the replacement types to be used."""
      self.types = match[0].split(',');
      """'z' for example. This is the current file's `type`."""
      self.precision = match[2].lower();
      """['c','d','s'] for example. This is the current file's destination `types`."""
      self.precisions = match[3].lower().split();
      if len(self.required_precisions):
        self.precstmp = [];
        for prec in self.required_precisions:
          if prec in self.precisions or prec == self.precision:
            self.precstmp.append(prec);
        self.precisions = self.precstmp;
    except:
      raise ValueError(path.join(self.file[0],self.file[1])+' : Invalid conversion string');
    self.files_in.append(rel);

  def run(self):
    """Does the appropriate work, if in test mode, this is limited to only converting names."""
    if self.convert_names() and not self.test:
      """If not in test mode, actually make changes and export to disk."""
      self.convert_data();
      self.export_data();

  def convert_names(self):
    """Investigate file name and make appropriate changes."""
    self.names = [];
    self.dates = [];
    self.copy = [];
    self.converted = [];
    load = False;
    if self.debug:
      print( '|'.join(self.types), self.precision, relpath(path.join(self.file[0],self.file[1])) );
    for precision in self.precisions:
      """For each destination precision, make the appropriate changes to the file name/data."""
      new_file = self.convert(self.file[1], precision);
      if self.debug:
        print(precision,':')
      copy = False;
      if new_file != self.file[1] or self.prefix is not None:
        if self.prefix is None:
          """If no prefix is specified, use the file's current folder."""
          prefix = ''
          makeprefix = '';
        else:
          """If a prefix is specified, set it up."""
          prefix = self.prefix;
          makeprefix = '--prefix '+prefix;
          if new_file == self.file[1]:
            copy = True;
        """Where the destination file will reside."""
        conversion = path.join(prefix, new_file);
        file_out = relpath(conversion);
        if self.make:
          """If in GNU Make mode, write the rule to create the file."""
          file_in = relpath(path.join(self.file[0],self.file[1]));
          print(file_out+':',file_in);
          print("\t$(PYTHON)",path.realpath(sys.argv[0]),makeprefix,'-p',precision,"--file",file_in);
        self.names.append(new_file);
        self.files_out.append(file_out);
        self.dependencies.append( (path.join(self.file[0],self.file[1]), precision, file_out) );
        if self.debug:
          print(relpath(conversion), ':')
        try:
          """Try to emulate Make like time based dependencies."""
          date = path.getmtime(conversion);
          diff = self.date - date;
          self.dates.append(diff);
          if self.debug:
            if diff > 0:
              print('Old')
            else:
              print('Current')
            print(diff);
          if diff > 0:
            load = True;
        except:
          if self.debug:
            print('Missing');
          self.dates.append(None);
          load = True;
      elif precision != self.precision :
        """There was no change in the file's name, thus,
        no work can be done without overwriting the original."""
        if self.debug:
          print('<No Change>',':')
        else:
          print >> sys.stderr, new_file, 'had no change for', precision;
        self.names.append(None);
        self.dates.append(None);
      self.copy.append(copy);
    return load;

  def export_data(self):
    """After all of the conversions are complete,
    this will write the output file contents to the disk."""
    for i in range(len(self.names)):
      name = self.names[i];
      data = self.converted[i];
      copy = self.copy[i];
      if copy:
        shutil.copy(self.files_in[i], self.files_out[i]);
        continue;
      if data is None or name is None: continue;
      fd = open(self.files_out[i], 'w');
      fd.write(data);
      fd.close();


  def convert_data(self):
    """Convert the data in the files by making the
    appropriate replacements for each destination precision."""
    for i in range(len(self.precisions)):
      precision = self.precisions[i];
      name = self.names[i];
      date = self.dates[i];
      copy = self.copy[i];
      if name is not None and not copy and (date is None or date > 0):
        self.converted.append(self.convert(self.content, precision));
      else: self.converted.append(None);

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
        print('Bad replacement pair ',i,'in',sub_type);
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
    """Replace the replacement keywork with one that signifies this is an output file,
    to prevent multiple replacement issues if run again."""
    data = re.sub(KEYWORD+' '+','.join(self.types)+'.*', DONE_KEYWORD+' '+precision+' '+datetime.now().ctime(), data);
    return data;
