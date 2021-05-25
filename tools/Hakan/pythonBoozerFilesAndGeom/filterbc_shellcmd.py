#!/usr/bin/env python
from __future__ import division
import sys, os
import argparse
from filterbigbcfile import filterbigbcfile

###############################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(description='bcfile filtering')
    parser.add_argument('infile', type=str)
    parser.add_argument('outfile', type=str)
    parser.add_argument('min_Bmn', type=float)
    return parser.parse_known_args(argv[1:])
###############################################################

arg, arg_remains = parse_args(sys.argv)

filterbigbcfile(arg.infile,arg.outfile,arg.min_Bmn)
