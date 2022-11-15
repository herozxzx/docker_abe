#!/bin/sh

echo "\n==============Starting Compile & Execute.=============="
echo "[Compiling]"
cd /pbc/scripts && make
echo "\n[Executing]"
/pbc/authorizer /pbc/param/a.param
echo "=======================================================\n"
