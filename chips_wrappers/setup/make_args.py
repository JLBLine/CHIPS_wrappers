import sys
import argparse
import shlex
from copy import deepcopy

SCALAR_DENSITY_CORRECTION = 2.6505481727808022

def get_args(argv=None):
    """Parse command line arugments using argparse. Returns the args"""
    
    if argv:
        ##user wants a specific set of args
        pass
    elif len(sys.argv) > 1:
        argv = sys.argv
    else:
        # is being called directly from nextflow
        argv = shlex.split("""${args}""")

    parser = argparse.ArgumentParser()

    file_group = parser.add_argument_group('FILE LOCATIONS')
    file_group.add_argument("--basedir", default='/astro/mwaeor/MWA/output/',
        help="The directory where the CHIPS outputs are located")
    file_group.add_argument("--chips_tag",
        help="Use this when doing a 2D power spectrum (not a ratio). "\
        "Alternatively, use this as the denominator in a ratio-diff 2D PS. "\
        "This is the tag that was passed to the 'lssa_thermal' step, e.g. if you "\
        "have an output called crosspower_xx_0.iter.real_cool_name.dat then "\
        "you should enter --chips_tag=real_cool_name")
    file_group.add_argument("--chips_tag_one",
        help="Use this data as the numerator in a 2D ratio/diff plot. " \
        "Alternatively, use this as the first numerator in a ratio-diff 2D PS. "\
        "This is the string that was passed to the 'lssa_thermal' step, e.g. if you "\
        "have an output called crosspower_xx_0.iter.real_cool_name.dat then "\
        "you should enter --chips_tag_one=real_cool_name")
    file_group.add_argument("--chips_tag_two",
        help="Use this data as the denominator in a 2D ratio/diff plot. " \
        "Alternatively, use this as the second numerator in a ratio-diff 2D PS. "\
        "This is thethat was passed to the 'lssa_thermal' step, e.g. if you "\
        "have an output called crosspower_xx_0.iter.real_cool_name.dat then "\
        "you should enter --chips_tag_two=real_cool_name")
    file_group.add_argument("--outputdir", default='./',
        help="Directory in which to output resultant plots. Default = './'")

    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument("--plot_type",
        help='Which type of plot to make. Options are: 1D, 2D, 2D_Ratio, '\
              '2D_diff, 2D_Ratio_Diff, 1D_comp. Defaults to 2D',default='2D')
    plot_group.add_argument("--polarisation", default='both',
        help='Which polarisation (XX or YY) to plot. Options are: xx, yy, both. Defaults to both')
    plot_group.add_argument("--plot_mode", default='png',
        help="Either 'png' or 'screen'. 'png' saves to .png, 'screen' does plt.show(). " \
        "Default = png")
    plot_group.add_argument("--max_power", default=0.0, type=float,
        help="Maximum power value to show. Defaults to 1e12 for a 2D or 1D" \
             "plot, and 10 for a 2D_Ratio plot. See --max_neg_power if" \
             "doing a 2D difference plot or using --colourscale=pos_and_negs")
    plot_group.add_argument("--min_power", default=0.0, type=float,
        help="Minimum power value to show. Defaults to 1e3 for a 2D or 1D" \
             "plot, and 0.1 for a 2D_Ratio plot. See --min_neg_power if" \
             "doing a 2D difference plot or using --colourscale=pos_and_negs")

    plot_group.add_argument("--max_power_r", default=1e4, type=float,
        help="Maximum power value to show for ratio sections of ratio-diff 2D plots." \
             "Defaults to 1e4 for a 2D_Ratio_Diff plot.")
    plot_group.add_argument("--min_power_r", default=1.0, type=float,
        help="Minimum power value to show for ratio sections of ratio-diff 2D plots." \
             "Defaults to 1 for a 2D_Ratio_Diff plot.")
    plot_group.add_argument("--max_power_rd", default=1e2, type=float,
        help="Maximum power value to show for ratio-diff section of ratio-diff 2D plots." \
             "Defaults to 1e2 for a 2D_Ratio_Diff plot.")
    plot_group.add_argument("--min_power_rd", default=1, type=float,
        help="Minimum power value to show for ratio-diff section of ratio-diff 2D plots." \
             "Defaults to 1 for a 2D_Ratio_Diff plot.")

    plot_group_2D = parser.add_argument_group('2D PLOTTING OPTIONS')
    plot_group_2D.add_argument("--colourscale", default='all_positive',
        help='How to handle negative powers in the colour scale for 2D plot.' \
        "Choices are: negs_are_grey (anything negative is blocked out as grey), " \
        "all_positive (anything below --min_power is same colour), "\
        "pos_and_negs (positive values are orange, negative purple)" \
        "Defaults to --colourscale=all_positive")
    plot_group_2D.add_argument("--max_neg_power", default=False, type=float,
        help="When plotting with --colourscale=pos_and_negs, use this value" \
        "as the upper limit of negative colourbar. Should be larger than "\
        "--min_neg_power. Default = negative --min_power")
    plot_group_2D.add_argument("--min_neg_power", default=False, type=float,
        help="When plotting with --colourscale=pos_and_negs, use this value" \
        "as the upper limit of negative colourbar. Should be smaller than "\
        "--max_neg_power. Default = negative --max_power")

    plot_group_1D = parser.add_argument_group('1D PLOTTING OPTIONS')
    plot_group_1D.add_argument("--no_noise", default=False, action='store_true',
        help="By default, thermal noise is plotted. Add this to switch it off")
    plot_group_1D.add_argument("--wedge_factor", default=-1, type=float,
        help="The scaling factor between k_parrallel and k_perpendicular" \
             "to use during cut. Defaults to the horizon")
    plot_group_1D.add_argument("--plot_delta", default=False, action='store_true',
        help="Plot in unitless Delta rather than P(k)")

    plot_group_1D.add_argument("--ktot_bin_edges", default=False,
        help="Path to text file containing edge values of k-bins for a 1D " \
             "plot. Overrides --low_k_edge, --high_k_edge, --num_k_edges." \
             "Will grid the data to the centre of each pair of bin edges, " \
             "e.g. if you provide 21 bin edges, your data will be gridded " \
             "to 20 bin centres." )

    plot_group_1D.add_argument("--low_k_edge", type=float, default=1e-2,
        help="Lowest k-value to grid the data to (applies to both perp and parra)")
    plot_group_1D.add_argument("--high_k_edge", type=float, default=5,
        help="Highest k-value to grid the data to (applies to both perp and parra)")
    plot_group_1D.add_argument("--num_k_edges", type=int, default=21,
        help="Number of k-bins to grid to between --low_k_edge and --high_k_edge")

    plot_group_1D.add_argument("--kperp_max", default=20, type=float,
        help="Maximum k_perp to use in 1D averaging. Can use this to "\
        "chop off small spatial scales. Defaults to 20")
    plot_group_1D.add_argument("--kperp_min", default=0, type=float,
        help="Minimum k_perp to use in 1D averaging. Can use this to "\
        "chop off large spatial scales. Defaults to no cut")

    plot_group_1D.add_argument("--kparra_min", default=0, type=float,
        help="Minimum k_parra to use in 1D averaging. Can use this to "\
        "chop off things close to the wedge that are leaking")

    plot_group_1D.add_argument("--plot_wedge_cut_2D", default=False,
        action='store_true',
        help="TAdd to plot the 2D power spectra with and without the cuts " \
             "that are being applied before binning to a 1D power spectra")

    plot_group_1D.add_argument("--chips_tag_label", default=False,
        help="When plotting two 1D power spectra on same axes, use this to " \
             "label the the data from --chips_tag")
    plot_group_1D.add_argument("--chips_tag_one_label", default=False,
        help="When plotting two 1D power spectra on same axes, use this to " \
             "label the the data from --chips_tag_one")
    plot_group_1D.add_argument("--chips_tag_two_label", default=False,
        help="When plotting two 1D power spectra on same axes, use this to " \
             "label the the data from --chips_tag_two")


    chips_group = parser.add_argument_group('CHIPS OPTIONS')
    chips_group.add_argument("--N_kperp",type=int, default=80,
        help="The number of kperp bins used in CHIPS 'fft_thermal' command. Default=80")
    chips_group.add_argument("--N_chan",type=int, default=384,
        help="The number of frequency channels used in CHIPS 'prepare_diff' and "\
        "'fft_thermal' commands. Default=384")
    ##TODO make options for low/high band that grab the frequency automagically for you
    chips_group.add_argument("--lowerfreq", default=167.035e6, type=float,
        help="Lowest frequency channel in data (Hz). Default is 167.035e+6")
    chips_group.add_argument("--chan_width", default=80e+3, type=float,
        help="Width of individual spectral channel (Hz). Default = 80e+3")
    # chips_group.add_argument("--deltat", default=8., type=float,
    #     help="Time resolution of data (s). Default = 8")
    chips_group.add_argument("--umax", default=300., type=float,
        help="Maximum u-value used in 'fft_thermal' stage (wavelengths). Default = 300")
    chips_group.add_argument("--density_correction", default=SCALAR_DENSITY_CORRECTION,
        help="Density correction to correct for decoherence. Defaults to factor 2 "
        "based on Barry et al. 2019a (Appendix A). Enter 0 for no correction, "
        "or 'use_fit' to use the fit from Line et al in prep")
    chips_group.add_argument("--num_obs", default=1, type=int,
        help="Total number of observations integrated by CHIPS, used in conjunction "
        "with --density_correction='use_fit'. Defaults to 1.")

    plot_group_1D.add_argument("--RTS_outputs", default=False,
        action='store_true',
        help="This script defaults to scaling the output power spectra " \
        "based on the RTS weighting scheme. If using data with a different " \
        "weighting scheme, add this flag to switch off this scaling")

    chips_group = parser.add_argument_group('COSMOLOGICAL CONSTANTS')
    chips_group.add_argument("--omega_matter", default=0.272, type=float,
        help="Matter density parameter. Default = 0.272")
    chips_group.add_argument("--omega_baryon", default=0.046, type=float,
        help="Baryon density parameter. Default = 0.046")
    chips_group.add_argument("--omega_lambda", default=0.7, type=float,
        help="Dark energy dentisty parameter. Default = 0.7")
    chips_group.add_argument("--hubble" ,default=70.4, type=float,
        help="Hubble param in km/s/Mpc, default=70.4")

    args = parser.parse_args()

    args.Neta = int(args.N_chan/2)
    args.verbose = True

    ##If people like the caps, let them eat cake
    if args.polarisation == 'XX': args.polarisation = 'xx'
    if args.polarisation == 'YY': args.polarisation = 'yy'

    return args

def check_args(args):
    """Given the inputs, check some certain values and make sure things make
    sense. Also auto-populate some values that will be useful later
    
    TODO so so many more sanity checks are needed"""
    
    if args.plot_type in ['2D_ratio_diff', '2D_Ratio_Diff', '2D_ratio',
                          '2D_Ratio', '2D_diff', '2D_Diff', '1D_comp', '1D_Comp']:
        if args.chips_tag_one_label and args.chips_tag_two_label:
            args.out_label1 = deepcopy(args.chips_tag_one_label).replace(' ', '_')
            args.out_label2 = deepcopy(args.chips_tag_two_label).replace(' ', '_')
        else:
            args.out_label1 = deepcopy(args.chips_tag_one).replace(' ', '_')
            args.out_label2 = deepcopy(args.chips_tag_two).replace(' ', '_')
            
    return args
    
    

class FakeArgs(object):
    """When not using argparse, use this object instead which can be called
    from jupyter notebook"""

    def __init__(self):
        self.N_kperp = None
        self.N_chan = None
        self.lowerfreq  = None
        self.chan_width = None
        self.umax = None
        self.omega_matter = None
        self.omega_baryon = None
        self.omega_lambda = None
        self.hubble = None
        self.Neta = None
        self.wedge_factor = -1.0
        self.basedir = None
        self.RTS_outputs = False
        self.max_1D_kmode = 20
        self.kperp_max = 20
        self.kperp_min = 0
        self.kparra_min = 0
        self.plot_wedge_cut_2D = False
        self.plot_type = '2D'
        self.max_power = 1e+13
        self.min_power = 1e+3
        self.num_k_edges = 31
        self.ktot_bin_edges = False
        self.low_k_edge = 1e-3
        self.high_k_edge = 30
        self.density_correction = SCALAR_DENSITY_CORRECTION
        self.verbose = True
        
        self.omega_matter=0.272
        self.omega_baryon=0.046
        self.omega_lambda=0.7
        self.hubble=70.4