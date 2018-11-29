#!/bin/bash

source activate
conda install muscle

# 
# prefix="muscle3.8.31_i86darwin64"
# if ! [ -x "$(command -v "$prefix")" ]
# then
#   cd ~/Downloads
#   curl "https://www.drive5.com/muscle/downloads3.8.31/"$prefix".tar.gz" -o $prefix".tar.gz"
#   tar xzvf $prefix".tar.gz"
# fi
