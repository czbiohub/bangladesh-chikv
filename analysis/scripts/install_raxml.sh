#!/bin/bash

if ! [ -x "$(command -v raxml-ng)" ]
then
  cd ~/Downloads
  mkdir raxml-ng_v0.7.0_macos_x86_64
  cd raxml-ng_v0.7.0_macos_x86_64
  curl "https://github-production-release-asset-2e65be.s3.amazonaws.com/75947982/3f284200-d654-11e8-85d3-c9ab23196daa?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIWNJYAX4CSVEH53A%2F20181126%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20181126T192740Z&X-Amz-Expires=300&X-Amz-Signature=ce3ae30cdb31e989c55a492eeb3cc4eb82f46f6e7b068295e496a53cc7bab8e1&X-Amz-SignedHeaders=host&actor_id=6191248&response-content-disposition=attachment%3B%20filename%3Draxml-ng_v0.7.0_macos_x86_64.zip&response-content-type=application%2Foctet-stream" -o "raxml-ng_v0.7.0_macos_x86_64.zip"
  unzip -a raxml-ng_v0.7.0_macos_x86_64.zip
  sudo cp raxml-ng /usr/local/bin
  sudo chmod +x /usr/local/bin/raxml-ng
fi