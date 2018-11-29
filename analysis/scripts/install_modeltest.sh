#!/bin/bash

if ! [ -x "$(command -v modeltest-ng)" ]
then
  brew install automake libtool
  sudo ln -s /usr/local/bin/glibtoolize /usr/local/bin/libtoolize
  cd ~/Downloads
  git clone --recursive https://github.com/ddarriba/modeltest
  cd modeltest
  mkdir build && cd build
  cmake ..
  make
  cd ../libs/pll-modules
  ./autogen.sh
  ./configure CPPFLAGS="-I./libs/libpll/src" LDFLAGS="-L./libs/libpll/src"
  make
  make install    # as root, otherwise run: sudo make install
  cd ../..
  autoreconf -i
  ./configure
  make
  make install
fi