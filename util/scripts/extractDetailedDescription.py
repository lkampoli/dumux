#! /usr/bin/python

# Extracts the paragraph with the detailed description of a class from
# a LaTeX file generated by doxygen. This script is used to avoid
# duplication between the API reference documentation and the DuMuX
# handbook.

import sys
import re

while True:
    curLine = sys.stdin.readline().strip()
    if re.search("\\\\[a-z]*section\{Detailed Description\}", curLine):
        # ignore one line
        i = 0
        while i < 1:           
            sys.stdin.readline().strip()
            i += 1

        # ignore all lines until we hit the next paragraph
        i = 0
        while i < 2:
            curLine = sys.stdin.readline().strip()
            if curLine == "":
                i += 1

        # print the next paragraph
        while True:
            curLine = sys.stdin.readline().strip()
            match = re.search("\\\\[a-z]*section|\\\\begin\\{", curLine)
            if match:
                print curLine[:match.start()]
                sys.exit(0)
            print curLine
        
