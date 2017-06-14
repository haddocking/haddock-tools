#!/usr/bin/env python
# coding=utf-8

# To support python 2.5+
from __future__ import with_statement

import argparse
import os
import json

import collections

"""
param_to_json.py

Convert a haddockparam.web in JSON format.
Can also be used to get or modify a parameter value.
Inspired by Sjoerd de Vries parser used in HADDOCK.
"""

__author__    = "Mikael Trellet"
__version__   = "1.0"
__copyright__ = "Copyright 2017, Apache 2"
__email__     = "mikael.trellet@gmail.com"
__credits__   = ['Mikael Trellet']


class HaddockParamWeb(object):

    def __init__(self, filename):
        self.filename = filename
        self.type = self._type()
        self.data = self._parse()
        self._cast_type()


    def _type(self):
        with open(self.filename, 'r') as f:
            s = f.read()
        return s.splitlines()[0].split()[0]


    def _parse(self):
        with open(self.filename, 'r') as f:
            s = f.read()
        s = s.rstrip('\n').rstrip() + ","
        stack = [([], True), ]
        curr, listmode = stack[-1]
        ident = 0
        objectlist = "ObjectList"
        for l in s.splitlines():
            if l.endswith("),"):
                if listmode is objectlist:  # we are parsing an objectlist
                    ll = l[ident + 2 * (curr[2] - 1):]
                    if ll == "),":  # dedent
                        if curr[2] > 0:
                            curr[1] += l[ident:] + "\n"
                            curr[2] -= 1
                            if curr[2] == 0:  # we are at outer level, parse what we have and reset
                                ## Seems to be never reached
                                txt = curr[1][:-2]  # we have removed outer indentation from txt and we are effectively doing a 2nd pass here (suboptimal)
                                print txt, curr[0]
                                curr[3].append(self._parse(curr[0]))
                                curr[:2] = None, None
                        else:  # dedent and leave objectlist mode
                            curr[:] = curr[3]
                            stack.pop()
                            curr, listmode = stack[-1]
                            ident -= 2
                    else:
                        if curr[2] == 0:  # parsing an elemental value, e.g Float(1), at outer level
                            ## Seems to be never reached
                            ind = ll.index("(")
                            typ = ll[2:ind]
                            val = ll[ind + 1:-2]
                            # print typ, val
                            # curr[3].append(spyder.__types__[typ](val))
                        else:  # elemental value at inner level, treat it as continuation
                            curr[1] += l[ident:] + "\n"
                else:  # no objectlist, dedent
                    stack.pop()
                    curr, listmode = stack[-1]
                    ident -= 2
            elif l.endswith(","):  # continuation
                if listmode == True:  # we're in an array, we expect unnamed items
                    v = l[ident:-1]
                    curr.append(v)
                elif listmode is objectlist:  # we're parsing an objectlist, gather it for later
                    curr[1] += l[ident:] + "\n"
                else:  # we're in a class, we expect named items
                    eq = l.index("=")
                    k = l[ident:eq - 1]
                    v = l[eq + 2:-1]
                    curr[k] = v
            else:  # endswith ( => indent
                if listmode == False:  # we're in a class, we expect named items
                    k = l[ident:l.find("=") - 1]
                    if l.find("ObjectList") > -1:  # enter objectlist mode for the inner level
                        new = [None, None, 0, []]
                        listmode = objectlist
                    elif l.endswith("Array ("):  # enter array mode for the inner level
                        new = []
                        listmode = True
                    else:  # enter class mode for the inner level
                        new = {}
                        listmode = False
                    stack.append((new, listmode))  # push the stack
                    curr[k] = new
                elif listmode is objectlist:  # parsing objectlist is just gathering text
                    if curr[2] == 0:  # reset the text if we are at the outer level
                        curr[0] = l[ident:-2]
                        curr[1] = ""
                    curr[1] += l[ident:] + "\n"
                    curr[2] += 1
                    ident -= 2  # to offset the +2 below, we don't want to increase ident read
                else:  # we're in array mode,
                    if l.find("ObjectList") > -1:  # enter objectlist mode for the inner level
                        new = [None, None, 0, []]
                        listmode = "ObjectList"
                    elif l.endswith("Array ("):  # enter array mode for the inner level
                        new = []
                        listmode = True
                    else:  # enter class mode for the inner level
                        new = {}
                        listmode = False
                    stack.append((new, listmode))
                    curr.append(new)
                ident += 2  # increase indentation read
                curr = new
        assert len(stack) == 1  # when we're done, the stack must have been popped
        return curr[0]

    def _change_value(self, key, new_val, dic):
        if hasattr(dic, 'iteritems'):
            for k, v in dic.iteritems():
                if k == key:
                    if type(new_val) != type(v):
                        raise Exception("Old and new values are not of the same type, {} expects {}".format(key, type(v)))
                    else:
                        if hasattr(new_val, "len") and hasattr(v, "len"):
                            if len(new_val) != len(v):
                                raise Exception("Old and new values have different length, {} is {} long".format(key, len(v)))
                            else:
                                print "Changing {} value from {} to {}".format(k, v, new_val)
                                dic[k] = new_val
                                yield dic
                        else:
                            print "Changing {} value from {} to {}".format(k, v, new_val)
                            dic[k] = new_val
                            yield dic
                if isinstance(v, dict):
                    for result in self._change_value(key, new_val, v):
                        yield result
                elif isinstance(v, list):
                    for d in v:
                        for result in self._change_value(key, new_val, d):
                            yield result

    def _get_value(self, key, dic=None):
        # if not dic:
        #     dic = self.data
        if hasattr(dic, 'iteritems'):
            for k, v in dic.iteritems():
                if k == key:
                    yield v
                if isinstance(v, dict):
                    for result in self._get_value(key, v):
                        yield result
                elif isinstance(v, list):
                    for d in v:
                        for result in self._get_value(key, d):
                            yield result

    def _cast_type(self, dic=None, key=None):
        if not dic:
            dic = self.data
        if hasattr(dic, 'iteritems'):
            for k, v in dic.iteritems():
                if isinstance(v, dict):
                    self._cast_type(v, k)
                elif isinstance(v, list):
                    c = 0
                    for s in v:
                        if isinstance(s, dict):
                            # print s
                            self._cast_type(s, k)
                        else:
                            dic[k][c] = eval(s)
                            c += 1
                else:
                    dic[k] = eval(v)
        elif isinstance(dic, list):
            for s in dic:
                self._cast_type(s)
        else:
            dic[key] = eval(dic)

    def update(self, orig_dict, new_dict):
        for key, val in new_dict.iteritems():
            if isinstance(val, collections.Mapping):
                tmp = self.update(orig_dict.get(key, {}), val)
                orig_dict[key] = tmp
            elif isinstance(val, list):
                orig_dict[key] = (orig_dict.get(key, []) + val)
            else:
                orig_dict[key] = new_dict[key]
        return orig_dict

    def cast_type(self, dic=None):
        if not dic:
            dic = self.data
        if hasattr(dic, 'iteritems'):
            for k, v in dic.iteritems():
                if isinstance(v, dict):
                    for result in self._get_value(v):
                        yield result
                elif isinstance(v, list):
                    for d in v:
                        for result in self._get_value(d):
                            yield result
                else:
                    dic[k] = eval(v)

    def change_value(self, key, new_val):
        # Check for the key existence first
        if not haddockparams.get_value(key):
            raise Exception("Key {} not found".format(key))
        # Try to change the value (must be a match between the old and new value types)
        try:
            result = list(haddockparams._change_value(key, new_val, self.data))
            if result:
                print "Updating dictionary..."
                self.data = result[0]
            else:
                raise Exception("An error was not caught during the update of the dictionary")
        except:
            raise

    def get_value(self, key):
        dic = self.data
        value = list(haddockparams._get_value(key, dic))
        if value:
            if len(value) > 1:
                return value
            else:
                return value[0]
        else:
            raise Exception("Key {} not found".format(key))

    def dump_keys(self, d, lvl=0):
        keys = {}
        for k, v in d.iteritems():
            print '%s%s' % (lvl * '  ', k)
            # keys.
            if type(v) == dict:
                self.dump_keys(v, lvl+1)


parser = argparse.ArgumentParser(description="This script parses a HADDOCK parameter file (*.web) and transforms it to "
                                             "JSON format.\n It also allows to change a parameter of the "
                                             "haddockparam.web")
parser.add_argument("web", nargs=1, help="HADDOCK parameter file")
parser.add_argument("-o", "--output", nargs=1, help="Path of JSON output file")
parser.add_argument("-g", "--get", nargs=1, help="Get value of a particular parameter")

args = parser.parse_args()

if os.path.exists(args.web[0]):
    haddockparams = HaddockParamWeb(args.web[0])
    #print haddockparams.data
    if args.get:
        haddockparams.get_value(args.get[0])
    if args.output:
        with open(args.output[0], 'w') as output:
            json.dump(haddockparams.data, output, indent=4, sort_keys=True)

# EXAMPLE - change waterrefine parameter from 400 to 200
print haddockparams.get_value('auto_his')
new_dic = haddockparams.update(haddockparams.data, {'dan1': {'constants': {'stages' : {'hot' : 200}}}})
#haddockparams.change_value('waterrefine', 200)
print haddockparams.get_value('hot')