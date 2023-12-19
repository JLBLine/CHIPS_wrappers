import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.ma import masked_array
import warnings
import matplotlib as mpl
import sys
from matplotlib.colors import LogNorm, Normalize

def save_or_plot(fig, output_plot_name, plot_mode):
    """Either save the figure to the given file name `output_plot_name`,
    or call plt.show(), based on whether plot_mode == 'png' or
    plot_mode == 'screen' or"""

    if plot_mode == 'png':
        print(f"Saving {output_plot_name}")
        plt.savefig(output_plot_name, bbox_inches='tight')
        plt.close()

    if plot_mode == 'screen':
        print("Printing to screen")
        plt.show()


def plot_wedge_cut(self):
    """Plot the 2D power, what was cut, retained, and which bins were used
    during the gridding to a 1D PS."""

    ##These indexes will represent where power was cut
    wedge_cut_2D = np.where(self.binning_array == -1.0)

    ##Make copies of the data to cut
    power_cut = deepcopy(self.crosspower)
    power_retained = deepcopy(self.crosspower)

    power_retained[wedge_cut_2D] = 0.0

    ##Just subtract off what was cut to see what was retained
    power_cut -= power_retained


    fig, axs = plt.subplots(2,2, figsize=(8,10))

    extent = [self.kper[0], self.kper[-1], self.kpa[0], self.kpa[-1]]

    ##Mask negative power, too much of a hassle otherwise
    twoD_ps_array = masked_array(self.crosspower, self.crosspower <= 0)
    power_cut = masked_array(power_cut, power_cut <= 0)
    power_retained = masked_array(power_retained, power_retained <= 0)

    ##log10 hates all the zeros and negatives etc so mute a tonne
    ##of warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        im0 = axs[0,0].imshow(np.log10(twoD_ps_array), origin='lower',
                                interpolation='none', extent=extent,
                                aspect='auto')

        im1 = axs[0,1].imshow(np.log10(power_cut), origin='lower',
                                interpolation='none', extent=extent,
                                aspect='auto')

        im2 = axs[1,0].imshow(np.log10(power_retained), origin='lower',
                                interpolation='none', extent=extent,
                                aspect='auto')

    ##Mask out bins that aren't used in gridding
    binning_array = masked_array(self.binning_array, self.binning_array == -1)

    im3 = axs[1,1].imshow(binning_array, origin='lower', interpolation='none',
                    extent=extent, aspect='auto')

    axs[0,0].set_title('All power')
    axs[0,1].set_title('Power cut')
    axs[1,0].set_title('Power retained')
    axs[1,1].set_title('Bins applied')

    for ax,im in zip(axs.flatten(), [im0, im1, im2, im3]):

        ax.set_xscale("log")
        ax.set_yscale("log")

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        label = "log10(P(k))"
        if im == im3: label = 'Bin number'
        cb = fig.colorbar(im, cax=cax,extend='both', label=label) # format='%.0e'

    # print(np.nanmax(binning_array))

    axs[1,0].set_xlabel(r"$k_\bot$")
    axs[1,1].set_xlabel(r"$k_\bot$")
    axs[0,0].set_ylabel("$k_\parallel$")
    axs[0,1].set_ylabel("$k_\parallel$")

    plt.tight_layout()
    print("Saving wedgecut_2D.png")
    fig.savefig('wedgecut_2D.png',bbox_inches='tight')
    plt.close()
    

def create_positives_cmap(args, set_under_colour=False):
    """Sets up a colourmap (cmap) with a log10 normalisation based on the input
    arguments. Optionally, sets anything under zero to a specific colour"""
    ##Setup up a bespoke colour map
    cmap = mpl.cm.get_cmap("Spectral_r").copy()

    if args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
        upper = np.ceil(np.log10(args.max_power_r))
        lower = np.ceil(np.log10(args.min_power_r))
    else:
        upper = np.ceil(np.log10(args.max_power))
        lower = np.ceil(np.log10(args.min_power))

    if set_under_colour:
        cmap.set_under(set_under_colour)
        cmap.set_bad(set_under_colour)
        ##Set zero as the lower bound, so that everything under is the desired
        ##colour
        bounds = [0]
    else:
        bounds = []
        cmap.set_bad(cmap(0))

    for bound in np.linspace(lower, upper, cmap.N):
        # bounds.append(10**bound)
        bounds.append(bound)

    # ticks = [10**tick for tick in np.arange(lower, upper + 1, 2)]
    ticks = [tick for tick in np.arange(lower, upper + 1, 2)]
    labels = ["1e+{:d}".format(int(tick)) for tick in ticks]

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    return cmap, norm, ticks, labels


def do_2D_axes_labels(ax, title, polarisation,
                      hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=False):
    """Makes the 2D axis log scale, sticks labels on the axes, sets a title"""

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    if polarisation == 'xx':
        title_pol = 'XX'
    else:
        title_pol = 'YY'

    ax.set_title(f'{title_pol} {title}',size=16)
    if not hide_k_par_label:
        ax.set_xlabel(r'k$_\bot$ ($h$Mpc$^{-1}$)',fontsize=18)
    if not hide_k_perp_label:
        ax.set_ylabel(r'k$_\parallel$ ($h$Mpc$^{-1}$)',fontsize=18)

    ax.set_xscale("log")
    ax.set_yscale("log")

def plot_2D_on_ax(twoD_ps_array, extent, ax, fig, polarisation,
                  args, hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=False):
    """Plot the 2D PS data in `twoD_ps_array`, which covers the kper/k_par
    coords in `extent`, on the axes `ax` on figure `fig`. Uses a single
    cmap"""

    if args.colourscale == "negs_are_grey":
        cmap, norm, ticks, ticklabels = create_positives_cmap(args, "Grey")

    else:
        cmap, norm, ticks, ticklabels = create_positives_cmap(args)

    ##Gotta mask out the negatives or matplotlib just wigs out and shows
    ##incorrect colours everywhere
    twoD_ps_array = masked_array(twoD_ps_array, twoD_ps_array <= 0)

    ##Stop warnings when plotting log10 of zero
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="invalid value encountered in log10")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="divide by zero encountered in log10")


        im = ax.imshow(np.log10(twoD_ps_array), cmap=cmap, origin='lower',
                       norm=norm, aspect='auto', extent=extent,
                       interpolation='none')

    ##Append a smaller axis to plot the colourbar on
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)

    cb = fig.colorbar(im, cax=cax, format='%.0e',extend='both',
                      ticks=ticks)
    cb.ax.set_yticklabels(ticklabels)

    if not hide_cbar_label:
        if args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
            cax.set_ylabel(r'Ratio $\Delta$ P(k) / P(k)',fontsize=14)
        else:
            cax.set_ylabel(r'P(k) mK$^2$ $h^{-3}$ Mpc$^3$',fontsize=14)
    cb.ax.tick_params(labelsize=11)

    if args.colourscale == "negs_are_grey":
        cax.text(1.13, 0.01, '$\leq0$', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes,
                fontsize=11)

    do_2D_axes_labels(ax, 'Crosspower', polarisation, hide_cbar_label, hide_k_par_label, hide_k_perp_label)

def plot_2D_on_ax_two_colour_bars(twoD_ps_array, extent, ax, fig, cax_pos,
                  cax_neg, polarisation, args, title, cmap='PurpOrang',
                  hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=False):
    """Plot the 2D PS data in `twoD_ps_array`, which covers the kper/k_par
    coords in `extent`, on the axes `ax` on figure `fig`. Uses two different
    cmaps, one for positive values, another for negative. Positive colourbar
    is plotted on axes `cax_pos`, negative on `cax_neg`."""

    pos_data = masked_array(twoD_ps_array, twoD_ps_array <= 0)
    neg_data = masked_array(twoD_ps_array, twoD_ps_array >= 0)

    ##==========================================================================
    ##Do the positive plot
    ##==========================================================================
    if cmap == 'PurpOrang':
        # cmap_pos = mpl.cm.get_cmap("Oranges").copy()
        cmap_pos = "Oranges"
    elif cmap == 'BlueRed':
        # cmap_pos = mpl.cm.get_cmap("Blues").copy()
        cmap_pos = "Blues"

    vmax = np.log10(args.max_power)
    vmin = np.log10(args.min_power)

    tick_high = np.floor(vmax)
    tick_low = np.ceil(vmin)

    inc = int(np.ceil((tick_high - tick_low) / 6))
    pos_ticks =  np.arange(tick_low, tick_high + 1, inc)

    # pos_bounds = np.linspace(lower, upper, 100)
    # pos_ticks =  np.arange(lower, upper + 1, 2)
    pos_labels = [f"1e+{int(tick)}" for tick in pos_ticks]

    ##Stop warnings when plotting log10 of zero
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="invalid value encountered in log10")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="divide by zero encountered in log10")

        im_pos = ax.imshow(np.log10(pos_data), cmap=cmap_pos, origin='lower',
                           vmin=vmin, vmax=vmax,
                           aspect='auto', extent=extent, interpolation='none')

    cb_pos = fig.colorbar(im_pos, cax=cax_pos, extend='both', ticks=pos_ticks)
    cb_pos.ax.set_yticklabels(pos_labels)
    cb_pos.ax.tick_params(labelsize=11)

    ##==========================================================================
    ##Do the negative plot
    ##==========================================================================

    if args.max_neg_power:
        neg_upper = abs(args.min_neg_power)
    else:
        neg_upper = args.max_power

    if args.min_neg_power:
        neg_lower = abs(args.max_neg_power)
    else:
        neg_lower = args.min_power

    if cmap == 'PurpOrang':
        # cmap_neg = mpl.cm.get_cmap("Purples_r").copy()
        cmap_neg = "Purples_r"
    elif cmap == 'BlueRed':
        # cmap_neg = mpl.cm.get_cmap("Reds_r").copy()
        cmap_neg = "Reds_r"

    vmin = -np.log10(neg_upper)
    vmax = -np.log10(neg_lower)

    tick_low = np.ceil(np.log10(neg_lower))
    tick_high = np.floor(np.log10(neg_upper))

    inc = int(np.ceil((tick_high - tick_low) / 6))
    neg_ticks =  -np.arange(tick_low, tick_high + 1, inc)

    neg_labels = [f"-1e+{int(abs(tick))}" for tick in neg_ticks]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="invalid value encountered in log10")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="divide by zero encountered in log10")

        im_neg = ax.imshow(-np.log10(np.abs(neg_data)), cmap=cmap_neg, origin='lower',
                           vmin=vmin, vmax=vmax,
                           aspect='auto', extent=extent, interpolation='none')

    cb_neg = fig.colorbar(im_neg, cax=cax_neg, extend='both',
                          ticks=neg_ticks)
    cb_neg.ax.set_yticklabels(neg_labels)
    cb_neg.ax.tick_params(labelsize=11)

    if not hide_cbar_label:
        if args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
            cax_pos.set_ylabel(r'Ratio $\Delta$ P(k) / P(k)',fontsize=14)
            cax_neg.set_ylabel(r'Ratio $\Delta$ P(k) / P(k)',fontsize=14)
        else:
            cax_pos.yaxis.set_label_coords(4.5,0.0)
            cax_pos.set_ylabel(r'P(k) mK$^2$ $h^{-3}$ Mpc$^3$',fontsize=14)
            cax_neg.yaxis.set_label_coords(4.5,0.0)
            cax_neg.set_ylabel(r'P(k) mK$^2$ $h^{-3}$ Mpc$^3$',fontsize=14)

    do_2D_axes_labels(ax, title, polarisation, hide_cbar_label, hide_k_par_label, hide_k_perp_label)

def setup_ax_and_cax_for_double_colourbar(args, fig, num_axes, polarisation):
    """When using two cmaps, sets up an axes and two colour bar axes (cax)
    based on how many axes you are going to plot"""

    if num_axes == 1:
        ax = fig.add_axes([0.16, 0.1, 0.6, 0.84])
        cax_height = 0.42
        cax_width = 0.035
        cax_left = 0.78
        cax_pos = fig.add_axes([cax_left, 0.52, cax_width, cax_height])
        cax_neg = fig.add_axes([cax_left, 0.1, cax_width, cax_height])

    elif num_axes == 2:
        ax_width = 0.32
        ax_height = 0.84
        cax_height = 0.42
        cax_width = 0.02

        if polarisation == 'xx':
            cax_xx_left = 0.42
            ax = fig.add_axes([0.09, 0.1, ax_width, ax_height])
            cax_pos = fig.add_axes([cax_xx_left, 0.52, cax_width, cax_height])
            cax_neg = fig.add_axes([cax_xx_left, 0.1, cax_width, cax_height])

        elif polarisation == 'yy':
            cax_yy_left = 0.88
            ax = fig.add_axes([0.55, 0.1, ax_width, ax_height])
            cax_pos = fig.add_axes([cax_yy_left, 0.52, cax_width, cax_height])
            cax_neg = fig.add_axes([cax_yy_left, 0.1, cax_width, cax_height])

    return ax, cax_pos, cax_neg


def setup_ax_and_cax_for_double_colourbar_ratio_diff(fig, polarisation):
    """When using two cmaps, sets up an axes and two colour bar axes (cax)
    based on how many axes you are going to plot"""

    if polarisation == 'both':
        ax_height = 0.4
        ax_width = 0.2
        cax_height = 0.2
        cax_width = 0.01

        ax_0_0 = fig.add_axes([0.06, 0.54, ax_width, ax_height])
        cax_0_0_pos = fig.add_axes([0.06+ax_width+0.01, 0.54+cax_height, cax_width, cax_height])
        cax_0_0_neg = fig.add_axes([0.06+ax_width+0.01, 0.54, cax_width, cax_height])
        ax_cx_00 = (ax_0_0, cax_0_0_pos, cax_0_0_neg)

        ax_0_1 = fig.add_axes([0.35, 0.54, ax_width, ax_height])
        cax_0_1_pos = fig.add_axes([0.35+ax_width+0.01, 0.54+cax_height, cax_width, cax_height])
        cax_0_1_neg = fig.add_axes([0.35+ax_width+0.01, 0.54, cax_width, cax_height])
        ax_cx_01 = (ax_0_1, cax_0_1_pos, cax_0_1_neg)

        ax_0_2 = fig.add_axes([0.66, 0.54, ax_width, ax_height])
        cax_0_2_pos = fig.add_axes([0.66+ax_width+0.01, 0.54+cax_height, cax_width, cax_height])
        cax_0_2_neg = fig.add_axes([0.66+ax_width+0.01, 0.54, cax_width, cax_height])
        ax_cx_02 = (ax_0_2, cax_0_2_pos, cax_0_2_neg)

        ax_1_0 = fig.add_axes([0.06, 0.08, ax_width, ax_height])
        cax_1_0_pos = fig.add_axes([0.06+ax_width+0.01, 0.08+cax_height, cax_width, cax_height])
        cax_1_0_neg = fig.add_axes([0.06+ax_width+0.01, 0.08, cax_width, cax_height])
        ax_cx_10 = (ax_1_0, cax_1_0_pos, cax_1_0_neg)

        ax_1_1 = fig.add_axes([0.35, 0.08, ax_width, ax_height])
        cax_1_1_pos = fig.add_axes([0.35+ax_width+0.01, 0.08+cax_height, cax_width, cax_height])
        cax_1_1_neg = fig.add_axes([0.35+ax_width+0.01, 0.08, cax_width, cax_height])
        ax_cx_11 = (ax_1_1, cax_1_1_pos, cax_1_1_neg)

        ax_1_2 = fig.add_axes([0.66, 0.08, ax_width, ax_height])
        cax_1_2_pos = fig.add_axes([0.66+ax_width+0.01, 0.08+cax_height, cax_width, cax_height])
        cax_1_2_neg = fig.add_axes([0.66+ax_width+0.01, 0.08, cax_width, cax_height])
        ax_cx_12 = (ax_1_2, cax_1_2_pos, cax_1_2_neg)

        ax_tup = (ax_cx_00, ax_cx_01, ax_cx_02, ax_cx_10, ax_cx_11, ax_cx_12)

    else:
        ax_height = 0.8
        ax_width = 0.2
        cax_height = 0.4
        cax_width = 0.01

        ax_0 = fig.add_axes([0.06, 0.1, ax_width, ax_height])
        cax_0_pos = fig.add_axes([0.06+ax_width+0.01, 0.1+cax_height, cax_width, cax_height])
        cax_0_neg = fig.add_axes([0.06+ax_width+0.01, 0.1, cax_width, cax_height])
        ax_cx_0 = (ax_0, cax_0_pos, cax_0_neg)

        ax_1 = fig.add_axes([0.35, 0.1, ax_width, ax_height])
        cax_1_pos = fig.add_axes([0.35+ax_width+0.01, 0.1+cax_height, cax_width, cax_height])
        cax_1_neg = fig.add_axes([0.35+ax_width+0.01, 0.1, cax_width, cax_height])
        ax_cx_1 = (ax_1, cax_1_pos, cax_1_neg)

        ax_2 = fig.add_axes([0.66, 0.1, ax_width, ax_height])
        cax_2_pos = fig.add_axes([0.66+ax_width+0.01, 0.1+cax_height, cax_width, cax_height])
        cax_2_neg = fig.add_axes([0.66+ax_width+0.01, 0.1, cax_width, cax_height])
        ax_cx_2 = (ax_2, cax_2_pos, cax_2_neg)

        ax_tup = (ax_cx_0, ax_cx_1, ax_cx_2)

    return ax_tup


def do_2D_plot(chips_data):
    """Given the user supplied `args`, plot a 2D power spectrum"""

    if chips_data.parser_args.max_power == 0.0:
        chips_data.parser_args.max_power = 1e+12
    if chips_data.parser_args.min_power == 0.0:
        chips_data.parser_args.min_power = 1e+3

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        # ##Read in data and convert to a 2D array for plotting
        # crosspower, weights = read_in_data(args, chips_data.parser_args.polarisation)
        # twoD_ps_array, extent = convert_to_2D_PS_array(crosspower, weights)

        twoD_ps_array, extent = chips_data.read_data_and_create_2Darray(chips_data.parser_args.polarisation)

        if chips_data.parser_args.colourscale == 'pos_and_negs':
            fig = plt.figure(figsize=(6,7))
            num_axes = 1
            ax, cax_pos, cax_neg = setup_ax_and_cax_for_double_colourbar(chips_data.parser_args,
                                                fig, num_axes, chips_data.parser_args.polarisation)

            plot_2D_on_ax_two_colour_bars(twoD_ps_array, extent, ax, fig,
                                cax_pos, cax_neg, chips_data.parser_args.polarisation, chips_data.parser_args,
                                'Crosspower')
        else:

            fig, ax = plt.subplots(1,1,figsize=(6,7))
            plot_2D_on_ax(twoD_ps_array, extent, ax, fig, chips_data.parser_args.polarisation, chips_data.parser_args)
            plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{chips_data.parser_args.polarisation}_{chips_data.parser_args.chips_tag}_crosspower.png"

    elif chips_data.parser_args.polarisation == 'both':

        twoD_ps_array_xx, extent_xx = chips_data.read_data_and_create_2Darray('xx')
        twoD_ps_array_yy, extent_yy = chips_data.read_data_and_create_2Darray('yy')

        if chips_data.parser_args.colourscale == 'pos_and_negs':

            fig = plt.figure(figsize=(11,7))

            num_axes = 2

            ax_xx, cax_pos_xx, cax_neg_xx = setup_ax_and_cax_for_double_colourbar(chips_data.parser_args,
                                                                   fig, num_axes, 'xx')
            plot_2D_on_ax_two_colour_bars(twoD_ps_array_xx, extent_xx,
                                ax_xx, fig, cax_pos_xx, cax_neg_xx, 'xx', chips_data.parser_args,
                                'Crosspower', hide_cbar_label=True)

            ax_yy, cax_pos_yy, cax_neg_yy = setup_ax_and_cax_for_double_colourbar(chips_data.parser_args,
                                                                   fig, num_axes, 'yy')
            plot_2D_on_ax_two_colour_bars(twoD_ps_array_yy, extent_yy,
                                ax_yy, fig, cax_pos_yy, cax_neg_yy, 'yy', chips_data.parser_args,
                                'Crosspower', hide_k_perp_label=True)

        else:
            fig, axs = plt.subplots(1,2,figsize=(12,7))
            plot_2D_on_ax(twoD_ps_array_xx, extent_xx, axs[0], fig, 'xx', chips_data.parser_args,
                          hide_cbar_label=True)
            plot_2D_on_ax(twoD_ps_array_yy, extent_yy, axs[1], fig, 'yy', chips_data.parser_args,
                          hide_k_perp_label=True)

            plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_{chips_data.parser_args.chips_tag}_crosspower.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument\n'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def make_2D_ratio(chips_data, polarisation):
    """Using input user arguments `args`, make a 2D Ratio array for the given
    `polarisation`"""
    twoD_ps_array_numer, extent_numer = chips_data.read_data_and_create_2Darray(polarisation,
                                                      chips_data.parser_args.chips_tag_one)

    twoD_ps_array_denom, extent_denom = chips_data.read_data_and_create_2Darray(polarisation,
                                                      chips_data.parser_args.chips_tag_two)
    ratio = twoD_ps_array_numer / twoD_ps_array_denom

    return ratio, extent_numer


def plot_2D_ratio_on_ax(twoD_ps_ratio_array, extent, ax, fig, polarisation,
                  args, hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=False):
    """Plot the 2D ratio PS data in `twoD_ps_ratio_array`, which covers the kper/k_par
    coords in `extent`, on the axes `ax` on figure `fig`"""

    if args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
        norm = LogNorm(args.min_power_rd,args.max_power_rd)
    else:
        norm = LogNorm(args.min_power,args.max_power)

    cmap = mpl.cm.get_cmap('RdBu').copy()
    # cmap = mpl.cm.get_cmap('RdBu')

    # if set_under_colour:
    # cmap.set_under('Grey')

    cmap.set_bad('Grey')

    im = ax.imshow(twoD_ps_ratio_array, cmap=cmap, origin='lower',
                   aspect='auto', extent=extent,
                   norm=norm, interpolation='none')
                   # vmin=args.min_power,vmax=args.max_power,

    ##Append a smaller axis to plot the colourbar on
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)

    # cb = fig.colorbar(im, cax=cax, format='%.0e',extend='both')
    cb = fig.colorbar(im, cax=cax, extend='both')

    if not hide_cbar_label:
        if args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
            cax.set_ylabel('Ratio Difference',fontsize=14)
        else:
            cax.set_ylabel('Ratio P(k) / P(k)',fontsize=14)

    cb.ax.tick_params(labelsize=11)

    if args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
        title = 'Ratio-Difference (R2-R1)'
    else:
        title = 'Ratio'

    do_2D_axes_labels(ax, title, polarisation, hide_cbar_label, hide_k_par_label, hide_k_perp_label)

def do_2D_ratio_plot(chips_data):
    """Given the user supplied `chips_data.parser_args`, plot a 2D power spectrum ratio"""

    if chips_data.parser_args.max_power == 0.0:
        chips_data.parser_args.max_power = 10
    if chips_data.parser_args.min_power == 0.0:
        chips_data.parser_args.min_power = -10

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_ratio_array, extent =  make_2D_ratio(chips_data, chips_data.parser_args.polarisation)


        fig, ax = plt.subplots(1,1,figsize=(6,7))
        plot_2D_ratio_on_ax(twoD_ps_ratio_array,
                            extent, ax, fig, chips_data.parser_args.polarisation, chips_data.parser_args)
        plt.tight_layout()
        
        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{chips_data.parser_args.polarisation}_" \
                            f"{chips_data.parser_args.out_label1}_" \
                            f"{chips_data.parser_args.out_label2}_ratio.png"
            

    elif chips_data.parser_args.polarisation == 'both':
        fig, axs = plt.subplots(1,2,figsize=(12,7))

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_ratio_array_xx, extent =  make_2D_ratio(chips_data, 'xx')

        plot_2D_ratio_on_ax(twoD_ps_ratio_array_xx,
                            extent, axs[0], fig, 'xx', chips_data.parser_args, hide_cbar_label=True)

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_ratio_array_yy, extent =  make_2D_ratio(chips_data, 'yy')

        plot_2D_ratio_on_ax(twoD_ps_ratio_array_yy,
                            extent, axs[1], fig, 'yy', chips_data.parser_args, hide_k_perp_label=True)

        plt.tight_layout()
        
        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_" \
                        f"{chips_data.parser_args.out_label1}_" \
                        f"{chips_data.parser_args.out_label2}_ratio.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument\n'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def make_2D_diff(chips_data, polarisation):
    """Using input user arguments `args`, make a 2D Ratio array for the given
    `polarisation`"""

    twoD_ps_array_one, extent_one = chips_data.read_data_and_create_2Darray(polarisation,
                                                    chips_data.parser_args.chips_tag_one)

    twoD_ps_array_two, extent_two = chips_data.read_data_and_create_2Darray(polarisation,
                                                    chips_data.parser_args.chips_tag_two)

    difference = twoD_ps_array_one - twoD_ps_array_two

    return difference, extent_one

def do_2D_diff_plot(chips_data):
    """Given the user supplied `args`, plot a 2D power spectrum ratio"""
    
    args = chips_data.parser_args

    if chips_data.parser_args.max_power == 0.0:
        chips_data.parser_args.max_power = 1e12
    if chips_data.parser_args.min_power == 0.0:
        chips_data.parser_args.min_power = 1e3

    if args.chips_tag_one_label and args.chips_tag_two_label:
        title = f'Difference\n{chips_data.parser_args.chips_tag_one_label} - {chips_data.parser_args.chips_tag_two_label}'
    else:
        title = 'Difference'


    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_diff_array, extent =  make_2D_diff(chips_data, chips_data.parser_args.polarisation)

        fig = plt.figure(figsize=(6,7))
        num_axes = 1
        ax, cax_pos, cax_neg = setup_ax_and_cax_for_double_colourbar(chips_data.parser_args,
                                            fig, num_axes, chips_data.parser_args.polarisation)

        plot_2D_on_ax_two_colour_bars(twoD_ps_diff_array, extent, ax, fig,
                            cax_pos, cax_neg, chips_data.parser_args.polarisation, args,
                            title, cmap='BlueRed')
        
        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{chips_data.parser_args.polarisation}_" \
                        f"{chips_data.parser_args.out_label1}_" \
                        f"{chips_data.parser_args.out_label2}_diff.png"

    elif chips_data.parser_args.polarisation == 'both':

        fig = plt.figure(figsize=(11,7))

        num_axes = 2

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_diff_array_xx, extent_xx =  make_2D_diff(chips_data, 'xx')

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_diff_array_yy, extent_yy =  make_2D_diff(chips_data, 'yy')

        ax_xx, cax_pos_xx, cax_neg_xx = setup_ax_and_cax_for_double_colourbar(args,
                                                               fig, num_axes, 'xx')
        plot_2D_on_ax_two_colour_bars(twoD_ps_diff_array_xx, extent_xx,
                            ax_xx, fig, cax_pos_xx, cax_neg_xx, 'xx', args,
                            title,cmap='BlueRed', hide_cbar_label=True)

        ax_yy, cax_pos_yy, cax_neg_yy = setup_ax_and_cax_for_double_colourbar(args,
                                                               fig, num_axes, 'yy')
        plot_2D_on_ax_two_colour_bars(twoD_ps_diff_array_yy, extent_yy,
                            ax_yy, fig, cax_pos_yy, cax_neg_yy, 'yy', args,
                            title,cmap='BlueRed', hide_k_perp_label=True)

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_" \
                           f"{chips_data.parser_args.out_label1}_" \
                           f"{chips_data.parser_args.out_label2}_diff.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument\n'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def make_2D_ratio_diff(chips_data, polarisation):
    """Using input user arguments `args`, make a 2D Ratio array for the given
    `polarisation`"""
    
    args = chips_data.parser_args

    twoD_ps_array, extent = chips_data.read_data_and_create_2Darray(polarisation,
                                                      args.chips_tag)

    twoD_ps_array_one, extent_one = chips_data.read_data_and_create_2Darray(polarisation,
                                                      args.chips_tag_one)

    twoD_ps_array_two, extent_two = chips_data.read_data_and_create_2Darray(polarisation,
                                                      args.chips_tag_two)

    frac_1 = (twoD_ps_array_one - twoD_ps_array)/twoD_ps_array
    frac_2 = (twoD_ps_array_two - twoD_ps_array)/twoD_ps_array
    ratio_diff = frac_2 - frac_1

    return frac_1, frac_2, ratio_diff, extent

def do_2D_ratio_diff_plot(chips_data):
    """Given the user supplied `args`, plot a 2D power spectrum ratio"""

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        pol = chips_data.parser_args.polarisation

        ##Read in data and convert to a 2D array for plotting
        frac_1, frac_2, ratio_diff, extent =  make_2D_ratio_diff(chips_data, pol)

        
        if chips_data.parser_args.chips_tag_one_label and chips_data.parser_args.chips_tag_two_label and chips_data.parser_args.chips_tag_label:
        
            chips_tag = deepcopy(chips_data.parser_args.chips_tag_label)
            chips_tag_one = deepcopy(chips_data.parser_args.chips_tag_one_label)
            chips_tag_two = deepcopy(chips_data.parser_args.chips_tag_two_label)
            
        else:
            
            chips_tag = deepcopy(chips_data.parser_args.chips_tag)
            chips_tag_one = deepcopy(chips_data.parser_args.chips_tag_one)
            chips_tag_two = deepcopy(chips_data.parser_args.chips_tag_two)
        
        label1 = rf' R1 = $\frac{{{chips_tag_one.replace("_", "ˍ")+" - "+chips_tag.replace("_", "ˍ")}}}{{{chips_tag.replace("_", "ˍ")}}}$'
        label2 = rf' R2 = $\frac{{{chips_tag_two.replace("_", "ˍ")+" - "+chips_tag.replace("_", "ˍ")}}}{{{chips_tag.replace("_", "ˍ")}}}$'
        
        if chips_data.parser_args.colourscale == 'pos_and_negs':
            fig = plt.figure(figsize=(18,7))
            ax_cx_0, ax_cx_1, ax_cx_2 = setup_ax_and_cax_for_double_colourbar_ratio_diff(fig, pol)

            plot_2D_on_ax_two_colour_bars(frac_1, extent,
                                ax_cx_0[0], fig, ax_cx_0[1], ax_cx_0[2], pol, chips_data.parser_args,
                                label1, hide_cbar_label=True, hide_k_par_label=True)

            plot_2D_on_ax_two_colour_bars(frac_2, extent,
                                ax_cx_1[0], fig, ax_cx_1[1], ax_cx_1[2], pol, chips_data.parser_args,
                                label2, hide_cbar_label=False, hide_k_par_label=True, hide_k_perp_label=True)

            plot_2D_on_ax_two_colour_bars(ratio_diff, extent, 
                                ax_cx_2[0], fig, ax_cx_2[1], ax_cx_2[2], pol, chips_data.parser_args,
                                'Ratio-Difference (R2-R1)', cmap='BlueRed', hide_cbar_label=False, hide_k_par_label=True, hide_k_perp_label=True)

            output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{pol}_"\
                            f"{chips_data.parser_args.chips_tag}_"\
                            f"{chips_data.parser_args.out_label1}_"\
                            f"{chips_data.parser_args.out_label2}_ratio_diff.png"

        else:
            fig, axs = plt.subplots(1,3,figsize=(18,7))
            plot_2D_on_ax(frac_1, extent, axs[0], fig, pol, chips_data.parser_args, hide_cbar_label=True, hide_k_par_label=False, hide_k_perp_label=False)
            do_2D_axes_labels(axs[0], label1, pol, hide_cbar_label=False,hide_k_par_label=False, hide_k_perp_label=False)

            plot_2D_on_ax(frac_2, extent, axs[1], fig, pol, chips_data.parser_args, hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=True)
            do_2D_axes_labels(axs[1], label2, pol, hide_cbar_label=False,hide_k_par_label=True, hide_k_perp_label=False)

            plot_2D_ratio_on_ax(ratio_diff, extent, axs[2], fig, pol, chips_data.parser_args, hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=True)
            plt.tight_layout()
            
            output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{pol}_"\
                            f"{chips_data.parser_args.chips_tag}_"\
                            f"{chips_data.parser_args.out_label1}_"\
                            f"{chips_data.parser_args.out_label2}_ratio_diff.png"

    elif chips_data.parser_args.polarisation == 'both':

        ##Read in data and convert to a 2D array for plotting
        frac_1_xx, frac_2_xx, ratio_diff_xx, extent_xx =  make_2D_ratio_diff(chips_data, 'xx')
        frac_1_yy, frac_2_yy, ratio_diff_yy, extent_yy =  make_2D_ratio_diff(chips_data, 'yy')

        if chips_data.parser_args.chips_tag_one_label and chips_data.parser_args.chips_tag_two_label and chips_data.parser_args.chips_tag_label:
        
            chips_tag = deepcopy(chips_data.parser_args.chips_tag_label)
            chips_tag_one = deepcopy(chips_data.parser_args.chips_tag_one_label)
            chips_tag_two = deepcopy(chips_data.parser_args.chips_tag_two_label)
            
        else:
            
            chips_tag = deepcopy(chips_data.parser_args.chips_tag)
            chips_tag_one = deepcopy(chips_data.parser_args.chips_tag_one)
            chips_tag_two = deepcopy(chips_data.parser_args.chips_tag_two)
        
        label1 = rf' R1 = $\frac{{{chips_tag_one.replace("_", "ˍ")+" - "+chips_tag.replace("_", "ˍ")}}}{{{chips_tag.replace("_", "ˍ")}}}$'
        label2 = rf' R2 = $\frac{{{chips_tag_two.replace("_", "ˍ")+" - "+chips_tag.replace("_", "ˍ")}}}{{{chips_tag.replace("_", "ˍ")}}}$'

        if chips_data.parser_args.colourscale == 'pos_and_negs':
            fig = plt.figure(figsize=(18,13))
            ax_cx_00, ax_cx_01, ax_cx_02, ax_cx_10, ax_cx_11, ax_cx_12 = setup_ax_and_cax_for_double_colourbar_ratio_diff(fig, 'both')

            plot_2D_on_ax_two_colour_bars(frac_1_xx, extent_xx,
                                ax_cx_00[0], fig, ax_cx_00[1], ax_cx_00[2], 'xx', chips_data.parser_args,
                                label1, hide_cbar_label=True, hide_k_par_label=True)

            plot_2D_on_ax_two_colour_bars(frac_2_xx, extent_xx,
                                ax_cx_01[0], fig, ax_cx_01[1], ax_cx_01[2], 'xx', chips_data.parser_args,
                                label2, hide_cbar_label=False, hide_k_par_label=True, hide_k_perp_label=True)

            plot_2D_on_ax_two_colour_bars(ratio_diff_xx, extent_xx, 
                                ax_cx_02[0], fig, ax_cx_02[1], ax_cx_02[2], 'xx', chips_data.parser_args,
                                'Ratio-Difference (R2-R1)', cmap='BlueRed', hide_cbar_label=False, hide_k_par_label=True, hide_k_perp_label=True)

            plot_2D_on_ax_two_colour_bars(frac_1_yy, extent_yy,
                                ax_cx_10[0], fig, ax_cx_10[1], ax_cx_10[2], 'yy', chips_data.parser_args,
                                label1, hide_cbar_label=True)

            plot_2D_on_ax_two_colour_bars(frac_2_yy, extent_yy,
                                ax_cx_11[0], fig, ax_cx_11[1], ax_cx_11[2], 'yy', chips_data.parser_args,
                                label2, hide_cbar_label=False, hide_k_perp_label=True)

            plot_2D_on_ax_two_colour_bars(ratio_diff_yy, extent_yy, 
                                ax_cx_12[0], fig, ax_cx_12[1], ax_cx_12[2], 'yy', chips_data.parser_args,
                                'Ratio-Difference (R2-R1)', cmap='BlueRed', hide_cbar_label=False, hide_k_perp_label=True)

            # plt.tight_layout()
            
            output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_" \
                            f"{chips_data.parser_args.out_label1}_" \
                            f"{chips_data.parser_args.out_label2}_ratio_diff.png"
        else:
            fig, axs = plt.subplots(2,3,figsize=(18,13))

            plot_2D_on_ax(frac_1_xx, extent_xx, axs[0, 0], fig, 'xx', chips_data.parser_args, hide_cbar_label=True, hide_k_par_label=True, hide_k_perp_label=False)
            do_2D_axes_labels(axs[0, 0], label1, 'xx', hide_cbar_label=True,hide_k_par_label=True, hide_k_perp_label=False)

            plot_2D_on_ax(frac_2_xx, extent_xx, axs[0, 1], fig, 'xx', chips_data.parser_args, hide_cbar_label=False, hide_k_par_label=True, hide_k_perp_label=True)
            do_2D_axes_labels(axs[0, 1], label2, 'xx', hide_cbar_label=False,hide_k_par_label=True, hide_k_perp_label=True)

            plot_2D_ratio_on_ax(ratio_diff_xx, extent_xx, axs[0, 2], fig, 'xx', chips_data.parser_args, hide_cbar_label=False, hide_k_par_label=True, hide_k_perp_label=True)

            plot_2D_on_ax(frac_1_yy, extent_yy, axs[1, 0], fig, 'yy', chips_data.parser_args, hide_cbar_label=True, hide_k_par_label=False, hide_k_perp_label=False)
            do_2D_axes_labels(axs[1, 0], label1, 'yy', hide_cbar_label=True,hide_k_par_label=False, hide_k_perp_label=False)

            plot_2D_on_ax(frac_2_yy, extent_yy, axs[1, 1], fig, 'yy', chips_data.parser_args, hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=True)
            do_2D_axes_labels(axs[1, 1], label2, 'yy', hide_cbar_label=False,hide_k_par_label=False, hide_k_perp_label=True)

            plot_2D_ratio_on_ax(ratio_diff_yy, extent_yy, axs[1, 2], fig, 'yy', chips_data.parser_args, hide_cbar_label=False, hide_k_par_label=False, hide_k_perp_label=True)

            plt.tight_layout()
            
            output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_" \
                            f"{chips_data.parser_args.out_label1}_" \
                            f"{chips_data.parser_args.out_label2}_ratio_diff.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument\n'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)