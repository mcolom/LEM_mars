#!/usr/bin/bash
set -e # stop on error
make OMP=1 -j4
./lem_bin -o ./ -g ./ -au 1.0 -i input_0.png -r 1.0 -m 0.5 -n 1.0 -e 1.0 -c 1.0 -s 1.0 -ni 100 -lo 0.0 -p 10.0 -lw 0.0  -ca 1.0  -cw 3.0 -cl 0.0 -wc 1.0 -z 1.0 -k 0.0 -q 0.0 -t 0.5 -d 30 -nw 40000 -wm -bw 2.0 -v

