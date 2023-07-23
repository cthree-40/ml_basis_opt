#!/bin/bash

# input is file name

tail -n 8 $1 | grep -v neo | grep -v end | grep -v S | grep -v H | grep -v "\*" | sed -e "s/1.0000D+00//"
