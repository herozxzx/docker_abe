#!/bin/sh

echo "\n==============Starting Compile & Execute.=============="
echo "[Compiling]"
cd ./scripts && make
cd ../
echo "\n[Executing]"
./authorizer ./param/a.param
echo "=======================================================\n"
