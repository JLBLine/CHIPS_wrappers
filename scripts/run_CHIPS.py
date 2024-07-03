#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
# from subprocess import call,check_output
import subprocess as sbp
from os import environ, getcwd
from os.path import exists
from sys import exit

from astropy.io import fits
from numpy import *


def read_env_variables(filename):
    '''Reads the environment variables stored in the
    cluster specific bash script'''
    lines = open(filename,'r').read().split('\n')
    for line in lines:
        if '=' in  line:
            if 'CODEDIR' in line: CODEDIR = line.split('=')[-1]
            elif 'OUTPUTDIR' in line: OUTPUTDIR = line.split('=')[-1]
            elif 'OBSDIR' in line: OBSDIR = line.split('=')[-1]
            elif 'INPUTDIR' in line: INPUTDIR = line.split('=')[-1]
            elif 'BEAMDIR' in line: BEAMDIR = line.split('=')[-1]
            # elif 'PBSDIR' in line: PBSDIR = line.split('=')[-1]
            elif 'PLOTSDIR' in line: PLOTSDIR = line.split('=')[-1]

        else:
            pass

    # return CODEDIR,OUTPUTDIR,OBSDIR,INPUTDIR,BEAMDIR,PBSDIR,PLOTSDIR
    return CODEDIR,OUTPUTDIR,OBSDIR,INPUTDIR,BEAMDIR,PLOTSDIR

def check_env_variables():
    '''Check that all the necessary CHIPS env variables are set, and exit ifread_env_variablesPBSDIR
    missing'''

    ##env variables that are necessary for CHIPS to run
    # env_values = ['CODEDIR', 'OBSDIR', 'INPUTDIR', 'BEAMDIR', 'PBSDIR']
    env_values = ['CODEDIR', 'OBSDIR', 'INPUTDIR', 'BEAMDIR']

    for env in env_values:
        try:
            keyout = environ[env]

        except KeyError:
            exit("CHIPS environment variable {%s} is missing. Things will go wrong without it so exiting now".format(env))

def write_sbatch_exports(outfile, args):
    """Write necessary exports to a text file - for Garrawarla, this is done
    via a module load, so can be skipped"""

    if args.cluster == 'garrawarla':
        pass
    else:
        outfile.write('export CODEDIR={:s}\n'.format(environ['CODEDIR']))
        outfile.write('export OBSDIR={:s}\n'.format(environ['OBSDIR']))
        outfile.write('export BEAMDIR={:s}\n'.format(environ['BEAMDIR']))
        # outfile.write('export PBSDIR={:s}\n\n'.format(environ['PBSDIR']))

    ##Override defaults with optional args (these point to defaults if
    ##nothing was supplied by user)
    outfile.write('export INPUTDIR={:s}\n'.format(args.data_dir))
    outfile.write('export OUTPUTDIR={:s}\n'.format(args.output_dir))

def write_cluster_specifics(outfile, args):
    """Write any necessary module loads and path file funnies that are specific
    to each cluster"""

    if args.cluster == 'ozstar':
        outfile.write('#SBATCH --partition=skylake\n')
        outfile.write('#SBATCH --account=oz048\n\n')
        #outfile.write('module load openblas/0.2.20\n')
        #outfile.write('module load gcc/6.4.0 openmpi/3.0.0\n\n')
        outfile.write('source /fred/oz048/jline/software/chips/module_load_chips.sh\n')
    if args.cluster == 'nt':
        outfile.write('#SBATCH --account=oz048\n\n')
        outfile.write('module use /fred/oz048/achokshi/software/modulefiles\n')
        outfile.write('module load chips/master\n\n')
    elif args.cluster == 'garrawarla':
        outfile.write('#SBATCH --partition=workq\n')
        outfile.write('#SBATCH --account=mwaeor\n\n')
        outfile.write('module use /astro/mwaeor/software/modulefiles\n')
        outfile.write('module load chips/v3.0\n\n')
        
    elif args.cluster == 'garrawarla_cmt':
        outfile.write('#SBATCH --partition=workq\n')
        outfile.write('#SBATCH --account=mwaeor\n\n')
        outfile.write('module use /astro/mwaeor/software/modulefiles\n')
        outfile.write('module load chips/cmt\n\n')
    
        
    else:
        exit('Specified cluster is not recognised: {:s}')

def make_grid_sbatch(obs=None, output_log_dir=None, args=None):

    freq_name = False
    ##check this in check_args and make into a range
    if args.coarse_band_subset:
        coarse_band_subset = args.coarse_band_subset
        freq_name = True
    ##There is a gap in ultra freq coverage due to FM band
    ##so default only process 20 bands
    elif args.band == 'ultra':
        coarse_band_subset = range(1, 21)
    else:
        coarse_band_subset = range(1, 25)

    if args.no_uvfits_check:
        args.low_freq = False
        args.high_freq = False
    ##Test if the uvfits files specified even exist
    else:
        ##Collect the lowest and highest freqs in all uvfits
        all_uvfits_freqs = []
        for coarse_band in coarse_band_subset:
            uvfits_path = "{:s}{:s}/{:s}/{:s}{:02d}.uvfits".format(args.data_dir,
                           obs, args.uvfits_dir, args.uvfits_tag,
                           coarse_band)

            if exists(uvfits_path):

                if args.frequencies_gathered:
                    pass
                else:
                    ##Ok, sometimes the order of uvfits files out of the RTS
                    ##has the frequency decrease with coarse band index, sometimes
                    ##increasing. So take a tally of all the frequencies in this
                    ##subset of coarse bands so we can add it to the name of the
                    ##file, giving no doubt as to what data has gone into it
                    with fits.open(uvfits_path) as hdu:
                        cent_freq = hdu[0].header['CRVAL4']
                        ##subtract one because this is one indexed not zero
                        cent_pix = hdu[0].header['CRPIX4'] - 1
                        freq_res = hdu[0].header['CDELT4']

                        num_freqs = hdu[0].data.data.shape[3]

                        low_freq = cent_freq - cent_pix*freq_res
                        high_freq = cent_freq + (num_freqs - cent_pix - 1)*freq_res

                        all_uvfits_freqs.append(low_freq / 1e+6)
                        all_uvfits_freqs.append(high_freq / 1e+6)

            else:
                exit("{:s} doesn't exist - CHIPS will fail. Check --uvfits_dir "\
                    "and --uvfits_tag make sense (and --data_dir if set). "\
                    "Exiting now".format(uvfits_path))

        if not args.frequencies_gathered:
            args.low_freq = min(all_uvfits_freqs)
            args.high_freq = max(all_uvfits_freqs)

            args.frequencies_gathered = True

    ##if obs contains the uvfits directory, need to strip it off for naming
    ##reasons
    obsname = obs.split('/')[0]

    ##Add in frequency range to name if running with a subset so the user is sure
    ##what freqs are being used. Only do it once, because if gridding multiple
    ##observations we'll keep adding everytime we make a new gridding script
    if freq_name and not args.freq_name_added:
        args.freq_name_added = True
        args.chips_tag += "_freqs_{:.3f}_{:.3f}MHz".format(args.low_freq, args.high_freq)

    outfile = open('{:s}/run_grid_{:s}_{:s}.sh'.format(output_log_dir,obsname,args.chips_tag) ,'w+')

    outfile.write('#!/bin/bash\n')
    outfile.write('#SBATCH --job-name="grid_{:s}_{:s}"\n'.format(obsname,args.chips_tag))
    outfile.write('#SBATCH --export=NONE\n')
    outfile.write('#SBATCH --time=00:10:00\n')
    outfile.write('#SBATCH --nodes=1\n')
    outfile.write('#SBATCH --cpus-per-task=10\n')
    outfile.write('#SBATCH --output=grid_{:s}_{:s}_%A_%a.out\n'.format(obsname,args.chips_tag))
    outfile.write('#SBATCH --error=grid_{:s}_{:s}_%A_%a.err\n'.format(obsname,args.chips_tag))
    outfile.write('#SBATCH --mem=30000\n')
    outfile.write('#SBATCH --array={:d}-{:d}\n'.format(coarse_band_subset[0], coarse_band_subset[-1]))

    ##Write lines specific to each cluster
    write_cluster_specifics(outfile, args)

    ##If the cluster has no 'module load', write out some paths explicitly
    write_sbatch_exports(outfile, args)

    outfile.write('export OMP_NUM_THREADS=10\n')
    outfile.write('GPUBOXN=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")\n\n')

    ##Setup the correct arguments to the CHIPS commands
    if args.band == 'ultra': band_num = 2
    elif args.band == 'low': band_num = 0
    elif args.band == 'high': band_num = 1

    if args.drips:
        command = "${CODEDIR}/gridvisdrips"
    else:
        if args.band == 'ultra':
	    #command = './gridvisultra'
            command = '${CODEDIR}/gridvisdiffbn'
        else:
            if args.cluster == 'garrawarla':
                command = "${CODEDIR}/grid_bn"
            else:
                command = "${CODEDIR}/gridvisdiff"

    cmd = f'srun --mem=30000 --export=ALL {command} '\
          f'{args.data_dir}{obs}/{args.uvfits_dir}/{args.uvfits_tag}${{GPUBOXN}}.uvfits '\
          f'{obsname} {args.chips_tag} {band_num} -f {args.field} '\
          f'-p {args.timeres:.3f} -c {args.freqres:.5f}'

    if args.base_freq:
        cmd += ' -n {:.5f}'.format(args.base_freq)
    elif args.low_freq:
        ##This will have been set if using a subset of frequencies,
        ##need it if running bands 12-24 for example
        cmd += ' -n {:.5f}'.format(args.low_freq*1e+6)

    outfile.write(cmd + '\n')

    return 'run_grid_{:s}_{:s}.sh'.format(obsname, args.chips_tag)


def make_grid_sbatch_singleuvfits(obs=None, output_log_dir=None, args=None,
                                  uvfits_path=False):

    freq_name = False
    # ##check this in check_args and make into a range
    # if args.coarse_band_subset:
    #     coarse_band_subset = args.coarse_band_subset
    #     freq_name = True
    # ##There is a gap in ultra freq coverage due to FM band
    # ##so default only process 20 bands
    # elif args.band == 'ultra':
    #     args.coarse_band_subset = range(1, 21)
    # else:
    #     args.coarse_band_subset = range(1, 25)
    
    if not uvfits_path:
        uvfits_path = "{:s}{:s}/{:s}/{:s}".format(args.data_dir,
                       obs, args.uvfits_dir, args.single_uvfits)

    if args.no_uvfits_check:
        args.low_freq = False
        args.high_freq = False
    ##Test if the uvfits files specified even exist
    else:
        ##Collect the lowest and highest freqs in all uvfits
        all_uvfits_freqs = []
        

        if exists(uvfits_path):

            if args.frequencies_gathered:
                pass
            else:
                ##Ok, sometimes the order of uvfits files out of the RTS
                ##has the frequency decrease with coarse band index, sometimes
                ##increasing. So take a tally of all the frequencies in this
                ##subset of coarse bands so we can add it to the name of the
                ##file, giving no doubt as to what data has gone into it
                with fits.open(uvfits_path) as hdu:
                    cent_freq = hdu[0].header['CRVAL4']
                    ##subtract one because this is one indexed not zero
                    cent_pix = hdu[0].header['CRPIX4'] - 1
                    freq_res = hdu[0].header['CDELT4']

                    num_freqs = hdu[0].data.data.shape[3]

                    low_freq = cent_freq - cent_pix*freq_res
                    high_freq = cent_freq + (num_freqs - cent_pix - 1)*freq_res

                    all_uvfits_freqs.append(low_freq / 1e+6)
                    all_uvfits_freqs.append(high_freq / 1e+6)

        else:
            exit("{:s} doesn't exist - CHIPS will fail. Check --uvfits_dir "\
                "and --uvfits_tag make sense (and --data_dir if set). "\
                "Exiting now".format(uvfits_path))

        if not args.frequencies_gathered:
            args.low_freq = min(all_uvfits_freqs)
            args.high_freq = max(all_uvfits_freqs)

            args.frequencies_gathered = True

    ##if obs contains the uvfits directory, need to strip it off for naming
    ##reasons
    obsname = obs.split('/')[0]

    # ##Add in frequency range to name if running with a subset so the user is sure
    # ##what freqs are being used. Only do it once, because if gridding multiple
    # ##observations we'll keep adding everytime we make a new gridding script
    # if freq_name and not args.freq_name_added:
    #     args.freq_name_added = True
    #     args.chips_tag += "_freqs_{:.3f}_{:.3f}MHz".format(args.low_freq, args.high_freq)

    outfile = open('{:s}/run_grid_{:s}_{:s}.sh'.format(output_log_dir,obsname,args.chips_tag) ,'w+')

    outfile.write('#!/bin/bash\n')
    outfile.write('#SBATCH --job-name="grid_{:s}_{:s}"\n'.format(obsname,args.chips_tag))
    outfile.write('#SBATCH --export=NONE\n')
    outfile.write('#SBATCH --time=01:00:00\n')
    outfile.write('#SBATCH --nodes=1\n')
    outfile.write('#SBATCH --cpus-per-task=36\n')
    outfile.write('#SBATCH --output=grid_{:s}_{:s}_%A_%a.out\n'.format(obsname,args.chips_tag))
    outfile.write('#SBATCH --error=grid_{:s}_{:s}_%A_%a.err\n'.format(obsname,args.chips_tag))
    outfile.write('#SBATCH --mem=30000\n')
    # outfile.write('#SBATCH --array={:d}-{:d}\n'.format(args.coarse_band_subset[0], args.coarse_band_subset[-1]))

    ##Write lines specific to each cluster
    write_cluster_specifics(outfile, args)

    ##If the cluster has no 'module load', write out some paths explicitly
    write_sbatch_exports(outfile, args)

    outfile.write('export OMP_NUM_THREADS=36\n')
    # outfile.write('GPUBOXN=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")\n\n')

    ##Setup the correct arguments to the CHIPS commands
    if args.band == 'ultra': band_num = 2
    elif args.band == 'low': band_num = 0
    elif args.band == 'high': band_num = 1

    if args.drips:
        command = "${CODEDIR}/gridvisdrips"
    else:
        if args.band == 'ultra':
	    #command = './gridvisultra'
            command = '${CODEDIR}/gridvisdiffbn'
        else:
            if args.cluster == 'garrawarla':
                command = "${CODEDIR}/grid_bn"
            else:
                command = "${CODEDIR}/gridvisdiff"

    cmd = f'srun --mem=30000 --export=ALL {command} '\
          f'{uvfits_path} '\
          f'{obsname} {args.chips_tag} {band_num} -f {args.field} '\
          f'-p {args.timeres:.3f} -c {args.freqres:.5f}'

    if args.base_freq:
        cmd += ' -n {:.5f}'.format(args.base_freq)
    elif args.low_freq:
        ##This will have been set if using a subset of frequencies,
        ##need it if running bands 12-24 for example
        cmd += ' -n {:.5f}'.format(args.low_freq*1e+6)

    outfile.write(cmd + '\n')

    return 'run_grid_{:s}_{:s}.sh'.format(obsname, args.chips_tag)

def write_clean_script(args, output_log_dir, no_delete_log=False):

    outfile = open('{:s}/run_clean_{:s}.sh'.format(output_log_dir,args.chips_tag),'w+')

    outfile.write('#!/bin/bash -l\n')
    outfile.write('#SBATCH --job-name="clean_{:s}"\n'.format(args.chips_tag))
    outfile.write('#SBATCH --export=ALL\n')
    outfile.write('#SBATCH --time=4:00:00\n')
    outfile.write('#SBATCH --nodes=1\n')
    outfile.write('#SBATCH --cpus-per-task=1\n')
    outfile.write('#SBATCH --output=clean_{:s}.%j.out\n'.format(args.chips_tag))
    outfile.write('#SBATCH --error=clean_{:s}.%j.err\n'.format(args.chips_tag))
    outfile.write('#SBATCH --mem=30000\n')

    ##Write lines specific to each cluster
    write_cluster_specifics(outfile, args)

    ##If the cluster has no 'module load', write out some paths explicitly
    write_sbatch_exports(outfile, args)

    outfile.write('time rm ${{OUTPUTDIR}}bv_freq*.{:s}.dat\n'.format(args.chips_tag))
    outfile.write('time rm ${{OUTPUTDIR}}noise_freq*.{:s}.dat\n'.format(args.chips_tag))
    outfile.write('time rm ${{OUTPUTDIR}}weights_freq*.{:s}.dat\n'.format(args.chips_tag))
    outfile.write('time rm ${{OUTPUTDIR}}bvdiff_freq*.{:s}.dat\n'.format(args.chips_tag))
    outfile.write('time rm ${{OUTPUTDIR}}noisediff_freq*.{:s}.dat\n'.format(args.chips_tag))

    if no_delete_log:
        pass
    else:
        outfile.write('rm {:s}/grid*{:s}*.out\n'.format(output_log_dir, args.chips_tag))
        outfile.write('rm {:s}/grid*{:s}*.err\n'.format(output_log_dir, args.chips_tag))
        outfile.write('rm {:s}/lssa_{:s}*.out\n'.format(output_log_dir, args.chips_tag))
        outfile.write('rm {:s}/lssa_{:s}*.out\n'.format(output_log_dir, args.chips_tag))
        #outfile.write('rm {:s}/*clean_{:s}*\n'.format(cwd, args.chips_tag))

    return 'run_clean_{:s}.sh'.format(args.chips_tag)

def make_lssa(pol=None, output_log_dir=None, args=None):

    outfile = open('{:s}/run_lssa_{:s}_{:s}.sh'.format(output_log_dir,args.chips_tag, pol),'w+')

    outfile.write('#!/bin/bash -l\n')
    outfile.write('#SBATCH --job-name="lssa_{:s}_{:s}"\n'.format(pol,args.chips_tag))
    outfile.write('#SBATCH --export=NONE\n')
    outfile.write('#SBATCH --time=4:00:00\n')
    outfile.write('#SBATCH --nodes=1\n')
    outfile.write('#SBATCH --cpus-per-task={:d}\n'.format(args.lssa_cores))
    outfile.write('#SBATCH --output=lssa_{:s}_{:s}.%j.out\n'.format(args.chips_tag, pol))
    outfile.write('#SBATCH --error=lssa_{:s}_{:s}.%j.err\n'.format(args.chips_tag, pol))
    outfile.write('#SBATCH --mem=30000\n')

    ##Write lines specific to each cluster
    write_cluster_specifics(outfile, args)

    ##If the cluster has no 'module load', write out some paths explicitly
    write_sbatch_exports(outfile, args)

    outfile.write('export OMP_NUM_THREADS={:d}\n\n'.format(args.lssa_cores))

    ##Setup the correct arguments to the CHIPS commands
    if args.band == 'ultra': band_num = 2
    elif args.band == 'low': band_num = 0
    elif args.band == 'high': band_num = 1

    ##Currently does not support running subsets of freqs from a single uvfits
    ##file so assume full bandwidth
    if args.single_uvfits or args.direct_path_list:
        num_coarse_bands = 24
        full_bandwidth = 30.72e+6
        num_chans = int(full_bandwidth / args.freqres)
    else:
        if args.coarse_band_subset:
            num_coarse_bands = args.coarse_band_subset[-1] - args.coarse_band_subset[0] + 1
        else:
            num_coarse_bands = 24

        full_bandwidth = num_coarse_bands * 1.28e+6
        num_chans = int(full_bandwidth / args.freqres)

    if args.drips:
        command1 = "${CODEDIR}/prepare_lssa_drips"
        command2 = "${CODEDIR}/lssa_fg_drips"
    else:
        if args.band == 'ultra':
            command1 = "${CODEDIR}/prepare_diff"
        else:
            command1 = "${CODEDIR}/prepare_diff"

        if args.cluster == 'garrawarla':
            command2 = "${CODEDIR}/fft_thermal"
        elif args.cluster == 'nt':
            command2 = "${CODEDIR}/lssa_fg_thermal"
        else:
            command2 = "${CODEDIR}/lssa_fg_simple"

    ##Argument in CHIPS is 0 = do kriging, 1= do not
    if args.no_krig:
        krig = 1

    else:
        krig = 0

    ##In the commnand below, 80 is the number of k bins, and 300 is the maximum uv value to grid up to
    cmd1 = "srun --mem=30000 --export=ALL {:s} {:s} {:d} 0 '{:s}' {:s} {:d}"\
               " -c {:.5f} -p {:.3f}".format(command1, args.chips_tag, num_chans,
               pol, args.chips_tag, band_num, args.freqres, args.timeres)
    
    cmd2 = "srun --mem=30000 --export=ALL {:s} {:s} {:d} 80 '{:s}' 300. {:s} {:d}"\
                " {:d} 0 -c {:.5f} -p {:.3f}".format(command2, args.chips_tag,
                num_chans, pol, args.chips_tag, krig, band_num, args.freqres,
                args.timeres)

    if args.base_freq:
        cmd1 += ' -n {:.5f}'.format(args.base_freq)
    elif args.low_freq:
        ##This will have been set if using a subset of frequencies,
        ##need it if running bands 12-24 for example
        cmd1 += ' -n {:.5f}'.format(args.low_freq*1e+6)

    outfile.write(cmd1 + '\n')

    if args.cluster == 'garrawarla_cmt':
        outfile.write("\nexport INPUTDIR=$OUTPUTDIR\n\n")

    outfile.write(cmd2 + '\n')

    return 'run_lssa_{:s}_{:s}.sh'.format(args.chips_tag, pol)

def get_parser():
    parser = argparse.ArgumentParser(description="Setup and/or run jobs to run CHIPS")

    parser.add_argument('--cluster', default='garrawarla',
                        help='Which cluster is being used. Default is garrawarla. ' \
                        'Current clusters: ozstar, nt, garrawarla')
    parser.add_argument('--data_dir', default=False,
                        help='If data does not live in generic /MWA/data dir, ' \
                        'enter the base directory here. Script will search for ' \
                        'uvfits files in /data_dir/obsID/uvfits_dir')
    parser.add_argument('--obs_list', default=None, required=True,
                        help='Text file list of all obs IDs to be processed')
    parser.add_argument('--uvfits_dir', default='/',
                        help='The name of the directory that the uvfits are stored in' \
                        ' - default is to be located be located in */MWA/data/obsID.' \
                        'This arg defaults to "/" in the case that your uvfits files'
                        'are not in a subdir')
    parser.add_argument('--uvfits_tag', default='uvdump_',
                       help="Add different if your uvfits don't follow RTS " \
                       "naming conventions - however they MUST end in *01.uvfits " \
                       "where 01 is band number. Name of uvfits will be interpreted "\
                       "as: '%%s%%02d.uvfits %%(args.uvfits_tag,band_number)' ")
    parser.add_argument('--single_uvfits', default=False,
                       help="Use with hyperdrive outputs. " \
                       "Instead of using `--uvfits_tag` with 24 files, " \
                       "point to a single uvfits file which contains all of the " \
                       "the data to be gridded. Still used in conjunction with  "\
                       "`--data_dir` so multiple obs can be specified with one name.")
    parser.add_argument('--direct_path_list', default=False, required=False,
        help='Instead of using --obs_list, --uvfits_dir, --data_dir and '
        ' --single_uvfits to point to uvfits files , just supply a list of paths '
        'to uvfits files, e.g. first row "/path/to/somename.uvfits" and '
        'second row could be "/path/to/othername.uvfits". These files MUST '
        'match the order in --obs_list')
    
    parser.add_argument('--band', default=None, required=True,
                        help='Which band you are processing: ultra, low, or high')
    parser.add_argument('--output_tag',default=None, required=True,
                        help='Tag to add to output names')
    parser.add_argument('--obs_range', default=None,
                       help='Range of observations within --obs_list to use - ' \
                       'enter 2 values separated by a comma (e.g --obs_range=0,1' \
                       ' to only process the first obs)')

    parser.add_argument('--debug', default=False, action='store_true',
                        help='Add to print out all subprocess commands run by this script')
    
    parser.add_argument('--output_dir', default=False,
                        help='Add this to save to a custom destination')
    parser.add_argument('--no_delete_log', default=False, action='store_true',
                        help='Default is to auto-delete the CHIPS error and output ' \
                        'logs during the clean script - add this to retain' \
                        ' (please consider deleting them ASAP once read)')
    parser.add_argument('--drips', default=False, action='store_true',
                        help='Add to run with DRIPS instead of CHIPS')
    parser.add_argument('--no_run', default=False, action='store_true',
                        help="Don't submit the jobs to the queue")
    parser.add_argument('--base_freq', default=0.0, type=float,
                        help="If using FHD or non-RTS base frequency, add the' \
                        ' base frequency here in Hz")
    parser.add_argument('--timeres', default=8.0, type=float,
                        help="Time resolution of data (s) - defaults to 8s")
    parser.add_argument('--freqres', default=80.0e+3, type=float,
                        help="Frequency resolution of the data (Hz), defaults to 80e+3Hz")
    parser.add_argument('--field', default=3, type=int,
                        help="Limit gridded to specific field. Use field=0 for " \
                        "EoR0 and field=1 for EoR1 (defaults to grid anything)")
    parser.add_argument('--no_uvfits_check', default=False, action='store_true',
                        help="By default script checks all uvfits exist and will ' \
                        ' not run if they are missing. Add option to switch this' \
                        ' function off")
#    parser.add_argument('--no_clean', default=False, action='store_true',
#                        help="By default intermediate ")
    parser.add_argument('--lssa_cores', default=16, type=int,
                        help="Number of cores to use in the lssa step - defaults to 16")
    parser.add_argument('--no_krig', default=False, action='store_true',
                        help="By default, CHIPS applies kriging to account for "\
                        "missing data. Add this to switch kriging off")
    parser.add_argument('--coarse_band_subset', default=False,
                        help="Enter a range of coarse bands (and therefore subset "\
                        "of frequencies) to run CHIPS on as low,high e.g. "\
                        "--coarse_band_subset=0,4 to process bands 1,2,3,4 (zero indexes, end number exclusive). "\
                        "By default, CHIPS is run on all 24 coarse bands for "\
                        "each observation. ")

    parser.add_argument('--only_gridding', default=False, action="store_true",
                        help="Add to ONLY run the gridding step - this is "\
                             "useful if you need to submit 100s of gridding "\
                             "jobs that exceed the maximum allowed number of "\
                             "jobs, or you don't want to spam the cluster with "\
                             "a bazzilion jobs. Run in conjuction with --obs_range "\
                             "as many times as is needed. Drop the --only_gridding "\
                             "on the last set of gridding to then run the PS "\
                             "estimation jobs.")

    return parser

def check_args(args):
    """Take the arguments from the command line parser and make sure we have
    what we need. Update some of the arguments after interrogating them"""

    ##Get optional arguments
    cluster = args.cluster
    debug = args.debug
    uvfits_tag = args.uvfits_tag
    no_delete_log = args.no_delete_log
    drips = args.drips

    base_freq = args.base_freq
    if base_freq: base_freq = float(base_freq)

    ##Test the necessary arguments have been passed
    try:
        obs_list = open(args.obs_list,'r').read().split('\n')
        obs_list = array([obs for obs in obs_list if obs != '' and obs != ' ' ])
    except:
        exit('Cannot read --obs_list=%s, please check file. Exiting now.' %obs_list)
        
    args.obs_list = obs_list
        
        
    if args.direct_path_list:
        try:
            direct_path_list = open(args.direct_path_list,'r').read().split('\n')
            direct_path_list = [path for path in direct_path_list if path != '' and path != ' ' ]
        except:
            exit('Cannot read --direct_path_list=%s, please check file. Exiting now.' %direct_path_list)
        
        args.direct_path_list = direct_path_list
        
    ##If there is an obs range specified, try and convert it to a range
    ##Exit if failure
    if args.obs_range:
        #Try and make a sensible obs range. Exit if failed
        try:
            low,high = map(int,args.obs_range.split(','))
            args.obs_range = range(low,high)
        except:
            exit('Cannot convert --obs_range={:s} into two numbers - please check your input. Exiting'.format(args.obs_range))
    else:
        args.obs_range = range(len(obs_list))

    ##TODO make sure `obs_range` doesn't lie outside length of `obs_list`

    ##Make sure we have the correct amount of dir slashes as expected throughout script
    if args.uvfits_dir[-1] == '/': args.uvfits_dir = args.uvfits_dir[:-1]

    ##Combine arguments into a CHIPS style tag

    if args.uvfits_dir == '':
        args.chips_tag = args.output_tag
    else:
        args.chips_tag = args.uvfits_dir + '_' + args.output_tag

    ##If user specifies a data directory, use that. Otherwise default to the
    ##CHIPS env variable DATADIR
    if args.data_dir:
        DATADIR = args.data_dir
    else:
        DATADIR = environ['DATADIR']

    if DATADIR[-1] == '/': pass
    else: DATADIR += '/'

    args.data_dir = DATADIR

    ##If user specifies an output directory, use that. Otherwise default to the
    ##CHIPS env variable OUTPUTDIR
    if args.output_dir:
        OUTPUTDIR = args.output_dir
    else:
        OUTPUTDIR = environ['OUTPUTDIR']

    if OUTPUTDIR[-1] == '/': pass
    else: OUTPUTDIR += '/'

    args.output_dir = OUTPUTDIR

    if args.coarse_band_subset:
        try:
            low,high = map(int, args.coarse_band_subset.split(','))
            args.coarse_band_subset = range(low + 1, high+1)
        except:
            exit('Cannot convert --coarse_band_subset={:s} into two numbers - please check your input. Exiting.'.format(args.coarse_band_subset))

    ##Flags used to keep track of information gathered from uvfits file headers
    args.freq_name_added = False
    args.frequencies_gathered = False

    if args.single_uvfits:
        if args.single_uvfits[-7:] != ".uvfits": args.single_uvfits += ".uvfits"

    return args

FALSE_RUN = False

def run_sbatch_get_jobID(sbatch_script_name, dependency=False):
    """Run the given `sbatch_script_name` using 'sbatch --parsable', which
    forces sbatch to only return the job number to stdout. Use
    `subprocess.run` to grab this output, and return the job number.
    Optionally add a dependency, which will add `--dependency:afterok:{dependency}`
    to the sbatch command"""


    command_list = ["sbatch", "--parsable"]

    if dependency:
        command_list.append("--dependency=afterok:{:s}".format(dependency))

    command_list.append(sbatch_script_name)

    command = " ".join(command_list)

    if FALSE_RUN:
        print(f"Command run: {command}")
        jobID = 'nothin'
        print(f'Jobs ID is {jobID}')

    else:

        ##Runs the command (which is made up from a list, sigh)
        process_out = sbp.run(command_list, capture_output=True, text=True)
        jobID = process_out.stdout.strip('\n')

        print(f"Command run: {command}")

        if process_out.returncode != 0:
            print('\tcommand returned non-zero returncode')
            print(f'\tstderr report: {process_out.stderr}')

        print(f'Job ID: {jobID}')

    return jobID

##Do the running
if __name__ == '__main__':

    ##Call the argument parser and get the args
    parser = get_parser()
    args = parser.parse_args()

    ##Check the arguments and ensure we have what we need
    args = check_args(args)

    ##At the moment, only garrawarla has a 'module load chips' which does
    ##all the export PATH and fun for you. If not using garrawarla, you need
    ##a certain number of env variables to be set. Check they exist if needed
    if args.cluster == 'garrawarla':
        ##Need some of the env variables to test paths, so call the module
        ##loads and what not on garrawarla

        sbp.call('module use /astro/mwaeor/software/modulefiles',shell=True)
        sbp.call('module load chips/v3.0',shell=True)
        
    elif args.cluster == 'garrawarla_cmt':
        ##I don't understand why, but if chips/v3.0 is loaded, calling
        ##the chips/v3.0 module load doesn't actually change the env variables
        ##So manually set a path here that we need changing
        os.environ['CODEDIR'] = "/astro/mwaeor/software/chips/chips_cmt/bin/"
    
    else:
        ##Check the necessary environment variables are set
        check_env_variables()

    #Make output logs directory
    output_log_dir = 'logs_' + args.output_tag
    if not os.path.exists(output_log_dir):
        os.mkdir(output_log_dir)
    #
    grid_jobs = []
    ##for all requested observations, write gridding sbatch files
    
    ##We are constructing arguemnts from obs and 

    for obs_ind in args.obs_range:
        obs = args.obs_list[obs_ind]
        
        ##If a list of direct paths exist, use those
        if args.direct_path_list:
            print(args.direct_path_list[obs_ind])
            grid_job = make_grid_sbatch_singleuvfits(obs=obs,
                            output_log_dir=output_log_dir,
                            args=args, uvfits_path=args.direct_path_list[obs_ind])
        
        ##Otherwise, piecemeil together a name from other arguments
        else:
        
            if args.single_uvfits:
                grid_job = make_grid_sbatch_singleuvfits(obs=obs,
                                output_log_dir=output_log_dir,
                                args=args)
            else:
                grid_job = make_grid_sbatch(obs=obs, output_log_dir=output_log_dir, args=args)
        grid_jobs.append(grid_job)

    ##Write a script to delete intermediate CHIPS products
    clean_job_name = write_clean_script(args=args,output_log_dir=output_log_dir,
                                        no_delete_log=args.no_delete_log)

    job_xx = make_lssa(pol='xx', output_log_dir=output_log_dir, args=args)
    job_yy = make_lssa(pol='yy', output_log_dir=output_log_dir, args=args)

    ##==================sbatch and do depends=================================##

    if args.no_run:
        pass
    else:
        # ##move to where we wrote out the sbatch jobs
        ##this means all the .out and .err from sbatch get put in the correct
        ##directory
        os.chdir(output_log_dir)

        ##Launch the first grid job and get the jobID to start the long list
        ##of dependencies
        # grid_jobs_IDs = run_sbatch_get_jobID(grid_jobs[0])
        grid_jobs_ID = run_sbatch_get_jobID(grid_jobs[0])


        for grid_job in grid_jobs[1:]:
            grid_jobs_ID = run_sbatch_get_jobID(grid_job, dependency=grid_jobs_ID)
            ##Add this job onto the dependencies for the last job
            
            ##grab the last job ID for a dependency
            # grid_jobs_IDs += ":{:s}".format(grid_jobs_ID)

        ##Don't do the PS estimation if asked to skip
        if args.only_gridding:
            pass
        else:

            ##once gridding is launched, launch the prepare_diff and fft_thermal jobs
            xx_jobs_ID = run_sbatch_get_jobID(job_xx, dependency=grid_jobs_ID)
            yy_jobs_ID = run_sbatch_get_jobID(job_yy, dependency=grid_jobs_ID)

            clean_depend = "{:s}:{:s}".format(xx_jobs_ID, yy_jobs_ID)

            clean_job = run_sbatch_get_jobID(clean_job_name, dependency=clean_depend)
