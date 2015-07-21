import argparse
import numpy as np
import sys
import numpy.linalg as la

parser = argparse.ArgumentParser(prog='gen_gold_standard.py',description='generate gold standard for validation')
parser.add_argument('dump_number',type=int)

args = parser.parse_args()

dump_num = args.dump_number

ux_fn = 'ux%d.b_dat'%dump_num
uy_fn = 'uy%d.b_dat'%dump_num
uz_fn = 'uz%d.b_dat'%dump_num

ux = np.fromfile(ux_fn, dtype=np.float32)
uy = np.fromfile(uy_fn, dtype=np.float32)
uz = np.fromfile(uz_fn, dtype=np.float32)

umag = np.sqrt(ux**2+uy**2+uz**2)

filename = 'gold_standard'
np.save(filename,umag)

