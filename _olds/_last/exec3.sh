#!/bin/sh

echo "\n==============(ABE with option)Starting Compile & Execute.==============="
echo "[Compiling]"
cd /pbc/abe_option/scripts && make
echo "\n[Executing]"
/pbc/abe_option/authorizer /pbc/param/a.param
echo "=========================================================================\n"
