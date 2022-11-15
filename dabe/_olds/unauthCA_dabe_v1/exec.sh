#!/bin/sh

echo "\n[Compiling]"
cd ./scripts && make
cd ../
echo "\n[Executing]"
./authorizer ./param/a.param