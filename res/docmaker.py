# A simple tool for generating Markdown documentation for Quaspy.

import os;
import types;

import inspect;

from inspect import signature;
from inspect import getmembers;
from inspect import ismodule;
from inspect import isfunction;
from inspect import isclass;
from inspect import ismethod;

import quaspy;

import quaspy.factoring;
import quaspy.factoring.general;
import quaspy.factoring.general.postprocessing;
import quaspy.factoring.general.postprocessing.ekera;
import quaspy.factoring.general.postprocessing.shor;
import quaspy.factoring.rsa;
import quaspy.factoring.rsa.postprocessing;
import quaspy.factoring.sampling;

import quaspy.logarithmfinding;
import quaspy.logarithmfinding.general;
import quaspy.logarithmfinding.general.postprocessing;
import quaspy.logarithmfinding.sampling;
import quaspy.logarithmfinding.short;
import quaspy.logarithmfinding.short.postprocessing;
import quaspy.logarithmfinding.short.sampling;

import quaspy.math;
import quaspy.math.continued_fractions;
import quaspy.math.crt;
import quaspy.math.groups;
import quaspy.math.kappa;
import quaspy.math.lattices.babai;
import quaspy.math.lattices.cvp;
import quaspy.math.lattices.enumerate;
import quaspy.math.lattices.gram_schmidt;
import quaspy.math.lattices.lagrange;
import quaspy.math.lattices.lll;
import quaspy.math.lattices.svp;
import quaspy.math.matrices;
import quaspy.math.modular;
import quaspy.math.primes;
import quaspy.math.random;
import quaspy.math.vectors;

import quaspy.orderfinding;
import quaspy.orderfinding.general;
import quaspy.orderfinding.general.postprocessing;
import quaspy.orderfinding.general.postprocessing.ekera;
import quaspy.orderfinding.general.sampling;

import quaspy.utils;
import quaspy.utils.timeout;
import quaspy.utils.timer;

def fixup(paragraph):
  paragraph = paragraph.replace("\n", " ").strip();

  while True:
    x = paragraph.replace("  ", " ");

    if x == paragraph:
      break;

    paragraph = x;

  # Citations
  citations = dict();

  citations["[EH17]"] = "https://doi.org/10.1007/978-3-319-59879-6_20";
  citations["[E19p]"] = "https://doi.org/10.48550/arXiv.1905.09084";
  citations["[E20]"] = "https://doi.org/10.1007/s10623-020-00783-2";
  citations["[E21]"] = "https://doi.org/10.1515/jmc-2020-0006";
  citations["[E21b]"] = "https://doi.org/10.1007/s11128-021-03069-1";
  citations["[E23p]"] = "https://doi.org/10.48550/arXiv.2309.01754";
  citations["[E24]"] = "https://doi.org/10.1145/3655026";
  citations["[E24t]"] = \
    "https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf";

  citations["[Shor94]"] = "https://doi.org/10.1109/SFCS.1994.365700";
  citations["[Shor97]"] = "https://doi.org/10.1137/S0097539795293172";

  citations["[Seifert01]"] = "https://doi.org/10.1007/3-540-45353-9_24";

  citations["[Babai86]"] = "https://doi.org/10.1007/BF02579403";
  citations["[LLL82]"] = "https://doi.org/10.1007/BF01457454";

  citations["factoritall repository"] = \
    "https://www.github.com/ekera/factoritall";
  citations["Factoritall repository"] = \
    "https://www.github.com/ekera/factoritall";

  citations["qunundrum repository"] = \
    "https://www.github.com/ekera/qunundrum";
  citations["Qunundrum repository"] = \
    "https://www.github.com/ekera/qunundrum";

  citations["fpLLL"] = "https://github.com/fplll/fplll";

  for key in citations.keys():
    citation = "[" + key + "](" + citations[key] + ")";
    paragraph = paragraph.replace(key, citation);

  paragraph = paragraph.replace("Quaspy", \
                                  "[Quaspy](https://github.com/ekera/quaspy)");

  # Remarks
  if paragraph.startswith("@remark "):
    paragraph = "> " + paragraph.split("@remark ")[1];

  return paragraph;


def get_paragraphs(s):
  s = s.replace("@param ", "\n@param ");
  paragraphs = "\n".join([x.strip() for x in s.split("\n")]);

  paragraphs = paragraphs.split("\n\n");
  paragraphs = [fixup(paragraph) for paragraph in paragraphs];

  return paragraphs;


def get_brief(f):
  if f.__doc__ == None:
    return None;

  paragraphs = get_paragraphs(f.__doc__);

  # Parse the @brief section
  offset = 0;

  if (not paragraphs[offset].startswith("@brief ")):
    raise Exception("Error: Expected the docstring for \"" + str(f) + \
                      "\" to start with @brief.");

  brief = paragraphs[offset].split(" ")[1:];
  brief = " ".join(brief);

  return brief;


def markdown_for_routine(f, path = "docs", c = None):
  if c == None:
    print("  Processing function:", f.__name__);
  else:
    print("    Processing method:", c.__name__ + "." + f.__name__);

  if os.path.exists(path):
    if not os.path.isdir(path):
      raise Exception("Error: The output directory is contaminated by " +
        "unknown activity.");
  else:
    os.mkdir(path);

  if f.__doc__ == None:
    brief = [];
    paragraphs = [];

    offset = 0;
  else:
    paragraphs = get_paragraphs(f.__doc__);

    # Parse the @brief section
    offset = 0;

    if (not paragraphs[offset].startswith("@brief ")):
      raise Exception("Error: Expected the docstring for the routine \"" +
                        str(f) + "\" to start with @brief.");

    brief = paragraphs[offset].split(" ")[1:];
    brief = " ".join(brief);

    offset += 1;

    brief = [brief];

    while offset < len(paragraphs):
      if paragraphs[offset].startswith("@"):
        break;

      if not paragraphs[offset].startswith("[["):
        brief.append(paragraphs[offset]);

      offset += 1;

  # Parse the @param section
  params = dict();
  param_names = [];

  while offset < len(paragraphs):
    if not paragraphs[offset].startswith("@param"):
      break;

    x = paragraphs[offset].split(" ");
    if len(x) < 2:
      raise Exception("Error: Expected parameter name after @param.");
    param_name = x[1];

    if param_name in params.keys():
      raise Exception("Error: The parameter '" + str(param_name) + \
        "' was twice documented.");

    param_description = [" ".join(x[2:])];

    offset += 1;

    while offset < len(paragraphs):
      if paragraphs[offset].startswith("@"):
        break;

      if not paragraphs[offset].startswith("[["):
        param_description.append(paragraphs[offset]);

      offset += 1;

    param_names.append(param_name);
    params[param_name] = param_description;

  # Parse the @return section
  if offset < len(paragraphs):
    if paragraphs[offset].startswith("@return"):
      ret = paragraphs[offset].split(" ")[1:];
      ret = [" ".join(ret)];

    offset += 1;

    while offset < len(paragraphs):
      if paragraphs[offset].startswith("@"):
        break;

      if not paragraphs[offset].startswith("[["):
        ret.append(paragraphs[offset]);

      offset += 1;
  else:
    ret = [];

  # Sanity check
  if offset != len(paragraphs):
    raise Exception("Error: Failed to completely parse the docstring.");

  # Write markdown
  file = open(os.path.join(path, f.__name__ + ".md"), "w");

  # Brief
  prototype_params = [];

  sgn = signature(f);
  sgn_param_names = list(sgn.parameters);
  for i in range(len(sgn_param_names)):
    param = sgn.parameters[sgn_param_names[i]];
    if param.default != inspect._empty:
      prototype_params.append("..");
      break;

    prototype_params.append(sgn_param_names[i]);

  prototype_params = ", ".join(prototype_params);

  if c == None:
    file.write("## Function: <code>" + f.__name__.replace("_", "\\_") + \
                "(" + prototype_params + ")" + "</code>\n");
  else:
    file.write("## Method: <code>" + c.__name__ + "." + \
                f.__name__.replace("_", "\\_") + \
                    "(" + prototype_params + ")" + "</code>\n");
  for x in brief:
    file.write(x + "\n\n");
  if brief == []:
    file.write("\n");

  # Import directive
  file.write("## Import directive\n");
  file.write("```python\n");
  if c == None:
    file.write("from " + str(f.__module__) + " import " + f.__name__ + "\n");
  else:
    file.write("from " + str(f.__module__) + " import " + c.__name__ + "\n");
  file.write("```\n\n");

  # Parent module
  if c == None:
    file.write("## Parent module\n");
    file.write("- [<code>" + str(f.__module__.split(".")[-1]) + \
                "</code>](README.md)\n\n");
  else:
    file.write("## Parent class\n");
    file.write("- [<code>" + str(c.__name__) + "</code>](" + "../" + \
                str(c.__name__) + ".md" + ")\n\n");

  # Prototype
  prototype = "def " + f.__name__ + "(";

  sgn = signature(f);
  sgn_param_names = list(sgn.parameters);

  for i in range(len(sgn_param_names)):
    prototype += sgn_param_names[i];

    param = sgn.parameters[sgn_param_names[i]];
    if param.annotation != inspect._empty:
      if isinstance(param.annotation, types.UnionType) or isinstance(param.annotation, str):
        prototype += " : " + str(param.annotation);
      else:
        prototype += " : " + str(param.annotation.__name__);
    if param.default != inspect._empty:
      prototype += " = " + str(param.default);

    if i + 1 < f.__code__.co_argcount:
      prototype += ",\n";
      for _ in range(4 + len(f.__name__) + 1):
        prototype += " ";
    else:
      prototype += ")";

  file.write("## Prototype\n");
  file.write("```python\n");
  file.write(prototype + "\n");
  file.write("```\n\n");

  # Parameters
  if param_names != []:
    file.write("## Parameters\n");

    file.write("| <b>Name</b> | <b>Description</b> |\n");
    file.write("| ----------- | ------------------ |\n");

    for param_name in param_names:
      param_description = params[param_name];

      description = "<br><br>".join(param_description);

      file.write("| " + param_name + " | " + description + " |\n");

    file.write("\n");

  if (len(sgn_param_names) > 0) and (sgn_param_names[0] == "self"):
      sgn_param_names = sgn_param_names[1:];
  if sgn_param_names != param_names:
    print("** Missing parameters:")
    print(sgn_param_names);
    print(param_names);

  # Return value
  if ret != []:
    file.write("## Return value\n");
    for x in ret:
      file.write(x + "\n\n");

  file.close();


def markdown_for_class(c, path = "docs"):
  print("  Processing class:", c.__name__);

  paragraphs = get_paragraphs(c.__doc__);

  # Parse the @brief section
  offset = 0;

  if (not paragraphs[offset].startswith("@brief ")):
    raise Exception("Error: Expected the docstring for the class \"" + str(c) +
                      "\" to start with @brief.");

  brief = paragraphs[offset].split(" ")[1:];
  brief = " ".join(brief);

  offset += 1;

  brief = [brief];

  while offset < len(paragraphs):
    if paragraphs[offset].startswith("@"):
      break;

    if not paragraphs[offset].startswith("[["):
      brief.append(paragraphs[offset]);

    offset += 1;

  # Methods
  methods = [];

  for function in getmembers(c, isfunction):
    # Only include functions actually defined in the class, not functions
    # defined in classes inherited from by this class.
    if function[0] in c.__dict__.keys():
      origin_module = c.__dict__[function[0]].__module__;
      if origin_module.startswith("quaspy"):
        methods.append(function[1]);

  for method in getmembers(c, ismethod):
    # Only include methods actually defined in the class, not methods
    # defined in classes inherited from by this class.
    if method[0] in c.__dict__.keys():
      origin_module = c.__dict__[method[0]].__module__;
      if origin_module.startswith("quaspy"):
        methods.append(method[1]);

  # Members
  members = [x[1] for x in getmembers(c) if x[0] == "__members__"];

  # Write markdown
  file = open(os.path.join(path, c.__name__ + ".md"), "w");

  # Brief
  file.write("## Class: <code>" + c.__name__ + "</code>\n");
  for x in brief:
    file.write(x + "\n\n");
  if brief == []:
    file.write("\n");

  # Import directive
  file.write("## Import directive\n");
  file.write("```python\n");
  file.write("from " + str(c.__module__) + " import " + c.__name__ + "\n");
  file.write("```\n\n");

  # Parent module
  file.write("## Parent module\n");
  file.write("- [<code>" + str(c.__module__.split(".")[-1]) + \
              "</code>](README.md)\n\n");

  # Methods
  if len(methods) > 0:
    file.write("## Methods\n");
    for m in methods:
      params = [];

      sgn = signature(m);
      sgn_param_names = list(sgn.parameters);
      for i in range(len(sgn_param_names)):
        param = sgn.parameters[sgn_param_names[i]];
        if param.default != inspect._empty:
          params.append("..");
          break;

        params.append(sgn_param_names[i]);

      params = ", ".join(params);

      file.write("- [<code>" + str(m.__name__).replace("__", "\\_\\_") + \
                  "(" + params + ")" + "</code>]" + \
        "(" + c.__name__ + "/" + m.__name__ + ".md)");
      if getattr(m, '__isabstractmethod__', False):
        file.write(" <code>abstract</code>");
      file.write("\n");

      x = get_brief(m);
      if None != x:
        file.write("\n  " + x + "\n\n");
      else:
        print("** Missing @brief documentation...");
        file.write("\n");

      markdown_for_routine(m, os.path.join(path, c.__name__), c);

  # Members
  if len(members) > 0:
    file.write("## Members\n");
    for m in members[0].keys():
      file.write("- <code>" + str(m) + "</code>\n");
    file.write("\n");

  file.close();


def markdown_for_module(module, path = "docs", level = 0):
  print("\nProcessing module: " + str(module.__name__));

  if module.__doc__ == None:
    paragraphs = [];
  else:
    paragraphs = get_paragraphs(module.__doc__);

  # Parse the @brief section
  offset = 0;

  if len(paragraphs) == 0:
    brief = [];
    print("** Missing @brief documentation...");
  else:
    if (paragraphs[offset].startswith("@brief ")):
      brief = paragraphs[offset].split(" ")[1:];
      brief = " ".join(brief);
      brief = [brief];
    else:
      brief = [paragraphs[offset]];
      print("** Missing @brief documentation...");

    offset += 1;

    while offset < len(paragraphs):
      if paragraphs[offset].startswith("@"):
        break;

      if not paragraphs[offset].startswith("[["):
        brief.append(paragraphs[offset]);

      offset += 1;

  submodules = [];

  for submodule in getmembers(module, ismodule):
    subpath = os.path.join(path, submodule[0]);
    submodule = submodule[1];

    if os.path.exists(subpath):
      if not os.path.isdir(subpath):
        raise Exception("Error: The output directory is contaminated by " +
          "unknown activity.");
    else:
      os.mkdir(subpath);

    if not submodule.__name__.startswith(module.__name__):
      continue;

    markdown_for_module(submodule, subpath, level + 1);
    submodules.append(submodule);

  functions = [];

  for f in getmembers(module, isfunction):
    if f[1].__module__ == module.__name__:
      markdown_for_routine(f[1], path);
      functions.append(f[1]);

  classes = [];

  for c in getmembers(module, isclass):
    if c[1].__module__ == module.__name__:
      markdown_for_class(c[1], path);
      classes.append(c[1]);

  # Write markdown
  file = open(os.path.join(path, "README.md"), "w");
  file.write("## Module: <code>" + str(module.__name__.split(".")[-1]) + \
              "</code>\n");
  for x in brief:
    file.write(x + "\n\n");
  if brief == []:
    file.write("\n");

  # Import directive
  file.write("## Import directive\n");
  file.write("```python\n");
  file.write("import " + str(module.__name__) + "\n");
  file.write("```\n\n");

  if level > 0:
    # Parent module
    file.write("## Parent module\n");
    parentname = module.__name__.split(".")[-2];
    file.write("- [<code>" + parentname + "</code>](../README.md)\n\n");

  # Submodules
  if len(submodules) > 0:
    file.write("## Submodules\n");
    for m in submodules:
      if len(m.__name__.split(".")) == 1:
        continue;

      submodulename = m.__name__.split(".")[-1];
      file.write("- [<code>" + submodulename + "</code>]" + \
        "(" + submodulename + "/README.md" + ")\n");

      x = get_brief(m);
      if None != x:
        file.write("\n  " + x + "\n\n");
      else:
        print("** Missing @brief documentation...");
        file.write("\n");

  # Functions
  if len(functions) > 0:
    file.write("## Functions\n");
    for f in functions:
      params = [];

      sgn = signature(f);
      sgn_param_names = list(sgn.parameters);
      for i in range(len(sgn_param_names)):
        param = sgn.parameters[sgn_param_names[i]];
        if param.default != inspect._empty:
          params.append("..");
          break;

        params.append(sgn_param_names[i]);

      params = ", ".join(params);

      file.write("- [<code>" + str(f.__name__) + "(" + params + ")</code>]" + \
        "(" + f.__name__ + ".md)\n");

      x = get_brief(f);
      if None != x:
        file.write("\n  " + x + "\n\n");
      else:
        print("** Missing @brief documentation...");
        file.write("\n");

  # Classes
  if len(classes) > 0:
    file.write("## Classes\n");
    for c in classes:
      file.write("- [<code>" + str(c.__name__) + "</code>]" + \
        "(" + c.__name__ + ".md)\n");

      x = get_brief(c);
      if None != x:
        file.write("\n  " + x + "\n");
      else:
        print("** Missing @brief documentation...");
        file.write("\n");

  file.close();

markdown_for_module(quaspy);
