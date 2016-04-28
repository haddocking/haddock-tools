#!/usr/bin/env python

import re
import argparse
import os

def validate_tbl(restraints, pcs = False):
    ret = ""
    parenthmatch = re.compile('[\(\)]')
    mode = "global"
    lines = restraints.replace('\r','').split("\n")
    lnr = 0
    post = None
    for l in lines:
        lnr +=1
        if l.find("!") > -1: l = l[:l.find("!")]
        l = l.strip()
        if not len(l):
            if mode == "postglobal" and not ret:
                ret += post
                post = None
            mode = "global" 
            continue
        if l.count('"') % 2 > 0:
            raise Exception('Unclosed ": %s' % l)
        if mode in ("global", "postglobal"):
            if l.lower().startswith("assi"):
                # if mode == "postglobal": 
                # ret += post
                if post:
                    ret += post
                mode = "assign"
                selections = []
                if l.find("(") == -1: continue
                post = None
            elif mode == "postglobal" and l.startswith("or"):
                mode = "postassign"
                selections = []
                l = l[len("or"):]
            else:
                raise Exception("Invalid TBL file: Unknown statement (line %d): %s" % (lnr, l))
        matched = True
        while matched:
            matched = False
            if mode in ("assign", "postassign"):
                pos = l.find("(")
                if pos != -1:
                    matched = True
                    l = l[pos+1:]
                    if mode == "postassign": mode = "postsel"
                    else: mode = "sel"
                    lastassign = lnr
                    s = ""
                    level = 1
            if mode in ("sel", "postsel"):
                for match in parenthmatch.finditer(l):
                    if match.group() == "(":
                        level += 1
                    else:
                        level -= 1
                    if level == 0:
                        if mode == "postsel": mode = "postassign"
                        else: mode = "assign"
                        matched = True
                        s += l[:match.start()]
                        selections.append(s)
                        s = None
                        l = l[match.end():]
                        break
                if level > 0: s += l + "\n"
        if mode in ("sel", "postsel"): continue
        if mode == "postassign":
            if len(selections) != postselections:
                raise Exception("Invalid TBL file: wrong number of selections: in-term %d, cross-term %d" % (postselections, len(selections)))
            if not post:
                post = " or"
            else:
                post += " or" 
            for s in selections: post += " (%s)\n" % s
            ret += post
            post = None
            mode = "postglobal"
        if len(l) == 0: continue
        if mode == "assign":
            mode = "numbers"
            if pcs:
                if len(selections) == 5: types = (" %.3f", " %.3f")
                else: raise Exception("Invalid TBL file: wrong number of selections (must be 5 in PCS mode)")
            else:
                if len(selections) == 2: types = (" %.3f", " %.3f", " %.3f")
                elif len(selections) == 4: types = (" %.3f", " %.3f", " %.3f", " %d")
                elif len(selections) == 5: raise Exception("Invalid TBL file: wrong number of selections (can be 5 only in PCS mode)")
                elif len(selections) == 6: types = (" %.3f", " %.3f")
                else: raise Exception("Invalid TBL file: wrong number of selections (must be 2,4 or 6)")
            postselections = len(selections)
            numbers = []
        if mode == "numbers":
            ll = l.split()
            for num in ll:
                if len(numbers) == len(types): break
                numbers.append(float(num))
            if len(numbers) == len(types):
                post = "assign "
                for s in selections: post += "(%s)\n" % s
                post = post[:-len("\n")]
                for n,t in zip(numbers,types): post += t % n
                post += "\n"
                mode = "postglobal"
    if mode == "postglobal":
        ret += post
        post = None
        mode = "global"
    if mode != "global":
        raise Exception("Invalid TBL file: Malformed ASSIGN statement (line %d)" % lastassign)
    if not len(ret.strip()):
        raise Exception("Invalid or empty TBL file")
    if post:
        ret += post
        post = None
    return ret

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script validates a restraint file (*.tbl).\n")
    
    parser.add_argument("file",
      help="TBL file to be validated")

    parser.add_argument("--pcs", action='store_true',
      help="PCS mode")
    
    args = parser.parse_args()

    if os.path.exists(args.file):
        tbldata = open(args.file).read()
        print validate_tbl(tbldata, args.pcs)
    else:
        raise Exception("TBL file %s does not exist, check the path" % args.file) 
