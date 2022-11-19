#!/bin/sh

echo "\n==============(DABE with option)Starting Compile & Execute.=============="
echo "[Compiling]"
cd /pbc/dabe_option/scripts && make
echo "\n[Executing]"
/pbc/dabe_option/authorizer /pbc/param/a.param
echo "=========================================================================\n"
