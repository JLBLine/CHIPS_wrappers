#!/usr/bin/env python


"""
This was taken from the excellent work of @jlbline
original location: /astro/mwaeor/jline/software/chips_2019/scripts/plotchips_all.py

example use with singularity:

```
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
# chips1D_xx+yy_eor0high_phase1-128T_p0_a3dbf06d_ionosub_30l_src4k_8s_80kHz.png
cp /astro/mwaeor/dev/nfdata/eor0high_phase1-128T_13d68053/ps_metrics/30l_src4k_8s_80kHz/*.dat .

eval singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/chips1D_tsv.py \
    --basedir "./"  \
    --chips_tag "eor0high_phase1-128T_13d68053_30l_src4k_8s_80kHz"  \
    --polarisation "both"  \
    --lowerfreq "166995000.0" \
    --umax "300" \
    --N_chan 384 \
    --num_k_edges "80" \
cp chips1D_xx+yy_eor0high_phase1-128T_p0_a3dbf06d_ionosub_30l_src4k_8s_80kHz.tsv /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/

export group=eor0high_phase2a-128T_410e2cca
export group=eor0high_phase2a-128T_8b6a256b
export name=ionosub_30l_src4k_300it_8s_80kHz_i1000
eval docker run -it --rm -v "/cygnus:/cygnus" -v "${HOME}/src:${HOME}/src" -w "${HOME}/src/MWAEoR-Pipeline" d3vnull0/ssins:latest \
    python templates/chips1D_tsv.py \
        --basedir /cygnus/dev/nfdata/\$group/ps_metrics/ionosub_30l_src4k_300it_8s_80kHz_i1000/ \
        --chips_tag "\${group}_\${name}" \
"""

from chips_wrappers.setup import make_args
from chips_wrappers.ps_methods.chips_data import ChipsDataProducts
import sys
import numpy as np
import pandas as pd

def save_1D(oneD_k_modes, oneD_power_measured, oneD_delta_measured,
            oneD_delta_2sig_noise, pol, delta=False):
    """Save everything to a tsv file"""

    notzero = np.where(oneD_k_modes != 0)
    df = pd.DataFrame({
        'k_modes': oneD_k_modes[notzero],
        'delta': oneD_delta_measured[notzero],
        'sqrt_delta': np.sqrt(oneD_delta_measured[notzero]),
        'power': oneD_power_measured[notzero],
        'noise': oneD_delta_2sig_noise[notzero],
    })
    filename=f"1D_power_{pol}.tsv"
    print(f"saving to {filename}")
    df.to_csv(f"1D_power_{pol}.tsv", index=False, sep=chr(9))

def do_1D_save(chips_data):
    args = chips_data.parser_args

    if args.polarisation == 'xx' or args.polarisation == 'yy':
        pols = [args.polarisation]
    elif args.polarisation == 'both':
        pols = ["xx", "yy"]
    else:
        msg = f'--polarisation={args.polarisation} is not a valid argument{chr(10)}'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)
    ##Read in data and convert to a 1D array for plotting
    for pol in pols:
        oneD_k_modes, oneD_delta_2sig_noise, oneD_power_measured, oneD_delta_measured = chips_data.read_data_and_create_1Darray(pol)
        save_1D(oneD_k_modes, oneD_power_measured, oneD_delta_measured, oneD_delta_2sig_noise, pol=pol, delta=args.plot_delta)

if __name__ == '__main__':

    args = make_args.get_args()
    print(vars(args))
    ptypes = [args.plot_type]
    chips_data = ChipsDataProducts(args)
    print(chips_data.parser_args)
    for ptype in ptypes:
        do_1D_save(chips_data)