#! /usr/bash

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./
eosmount ~/eos
.g++ makeRegTree.cc -I `root-config --incdir` -L. -ldictionary_cc `root-config --libs` -lProof  -O2 -pipe -Wall -W -Woverloaded-virtual  `root-config --ldflags` -fPIC -o makeRegTree
