#!/bin/sh

echo "\n====================(DABE)Starting Compile & Execute.===================="
echo "[Compiling]"
cd /pbc/dabe/scripts && make
echo "\n[Executing]"
/pbc/dabe/authorizer /pbc/param/a.param
echo "=========================================================================\n"
