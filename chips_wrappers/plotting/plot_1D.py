
from chips_wrappers.plotting.plot_2D import save_or_plot
import matplotlib.pyplot as plt
import numpy as np
import sys


def plot_1D_on_ax(ax, oneD_k_modes, oneD_power_measured, oneD_delta_measured,
                  oneD_delta_2sig_noise,
                  min_power, max_power,
                  pol, colour_power='C0', marker_power='x',
                  label='', delta=False, noise=True):
    """Plot the power and two sigma noise on the given axes"""

    notzero = np.where((oneD_k_modes != 0) & (np.isnan(oneD_delta_measured) == False))
    plot_k_modes = oneD_k_modes[notzero]
    plot_delta = oneD_delta_measured[notzero]
    plot_power = oneD_power_measured[notzero]
    plot_noise = oneD_delta_2sig_noise[notzero]

    if pol == 'xx':
        pol_label = f'XX {label}'
    else:
        pol_label = f'YY {label}'

    if delta:
        ax.plot(plot_k_modes, plot_delta, drawstyle='steps-mid',
               color=colour_power, marker=marker_power, mfc='none', ms=4,
               label=f'{pol_label}Power')
        ax.set_ylabel(r'$\Delta$ (mK$^2$)', fontsize=16)

        if noise:
            ax.plot(plot_k_modes, plot_noise, linestyle='--',
                color=colour_power, label=f'{pol_label}Noise')
    else:
        ax.plot(plot_k_modes, plot_power, drawstyle='steps-mid',
               color=colour_power, marker=marker_power, mfc='none', ms=4,
               label=f'{pol_label}Power')

        if noise:
            ax.plot(plot_k_modes, (plot_noise*2*np.pi**2)/plot_k_modes**3, linestyle='--',
                color=colour_power, label=f'{pol_label}Noise')

        ax.set_ylabel(r'P(k) mK$^2$ $h^{-3}$ Mpc$^3$',fontsize=16)

    np.savez_compressed(f"1D_power_{pol}.npz", k_modes=plot_k_modes, power=plot_power,
                                        delta=plot_delta, delta_error=plot_noise)



    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlabel(r'k ($h$Mpc$^{-1}$)',fontsize=16)

    ax.set_xscale("log")
    ax.set_yscale("log")

    if max_power:
        ax.set_ylim(top=max_power)
    if min_power:
        ax.set_ylim(bottom=min_power)

    ax.legend(prop={"size":14}, loc='best')

def do_1D_plot(chips_data):
    """Plot a 1D power spectrum given the input arguments from the parser"""

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        ##Read in data and convert to a 1D array for plotting
        oneD_k_modes, oneD_delta_2sig_noise, oneD_power_measured, oneD_delta_measured =  chips_data.read_data_and_create_1Darray(chips_data.parser_args.polarisation)


        fig, ax = plt.subplots(1,1,figsize=(8,6))

        plot_1D_on_ax(ax, oneD_k_modes, oneD_power_measured,
                      oneD_delta_measured, oneD_delta_2sig_noise,
                      chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                      chips_data.parser_args.polarisation,
                      delta=chips_data.parser_args.plot_delta,
                      noise=not(chips_data.parser_args.no_noise))

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_{chips_data.parser_args.polarisation}_{chips_data.parser_args.chips_tag}.png"

    elif chips_data.parser_args.polarisation == 'both':
        ##Read in data and convert to a 1D array for plotting
        oneD_k_modes_xx, oneD_delta_2sig_noise_xx, oneD_power_measured_xx, oneD_delta_measured_xx  =  chips_data.read_data_and_create_1Darray("xx")

        oneD_k_modes_yy, oneD_delta_2sig_noise_yy, oneD_power_measured_yy, oneD_delta_measured_yy =  chips_data.read_data_and_create_1Darray("yy")


        fig, axs = plt.subplots(2,1,figsize=(10,10))

        plot_1D_on_ax(axs[0], oneD_k_modes_xx, oneD_power_measured_xx,
                              oneD_delta_measured_xx, oneD_delta_2sig_noise_xx,
                              chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                              "xx", delta=chips_data.parser_args.plot_delta,
                              noise=not(chips_data.parser_args.no_noise))

        plot_1D_on_ax(axs[1], oneD_k_modes_yy, oneD_power_measured_yy,
                              oneD_delta_measured_yy, oneD_delta_2sig_noise_yy,
                              chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                              "yy", delta=chips_data.parser_args.plot_delta,
                              noise=not(chips_data.parser_args.no_noise))

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_xx+yy_{chips_data.parser_args.chips_tag}.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument\n'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def plot_1D_comparison(chips_data, pol, ax):
    """Given the axes `ax` and user supplied `chips_data.parser_args`, plot
    a comparison for the given polarisation `pol`"""

    ##Read in data and convert to a 1D array for plotting
    oneD_k_modes_1, oneD_delta_2sig_noise_1, oneD_power_measured_1, oneD_delta_measured_1 =  chips_data.read_data_and_create_1Darray(pol,
                                            chips_tag=chips_data.parser_args.chips_tag_one)


    if chips_data.parser_args.chips_tag_one_label:
        label_one = chips_data.parser_args.chips_tag_one_label + ' '
    else:
        label_one = chips_data.parser_args.chips_tag_one + ' '

    plot_1D_on_ax(ax, oneD_k_modes_1, oneD_power_measured_1,
                  oneD_delta_measured_1, oneD_delta_2sig_noise_1,
                  chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                  pol, label=label_one,
                  delta=chips_data.parser_args.plot_delta,
                  noise=not(chips_data.parser_args.no_noise))


    if chips_data.parser_args.chips_tag_two_label:
        label_two = chips_data.parser_args.chips_tag_two_label + ' '
    else:
        label_two = chips_data.parser_args.chips_tag_two + ' '

    oneD_k_modes_2, oneD_delta_2sig_noise_2, oneD_power_measured_2, oneD_delta_measured_2 =  chips_data.read_data_and_create_1Darray(pol,
                                            chips_tag=chips_data.parser_args.chips_tag_two)

    plot_1D_on_ax(ax, oneD_k_modes_2, oneD_power_measured_2,
                  oneD_delta_measured_2, oneD_delta_2sig_noise_2,
                  chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                  pol, label=label_two,
                  marker_power='o',colour_power='C1', delta=chips_data.parser_args.plot_delta,
                  noise=not(chips_data.parser_args.no_noise))

def do_1D_comparison_plot(chips_data):
    """Given the user supplied `chips_data.parser_args`, plot 2 different
    1D power spectra on the same axes"""

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        fig, ax = plt.subplots(1,1,figsize=(8,6))

        plot_1D_comparison(chips_data, chips_data.parser_args.polarisation, ax)

        plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_{chips_data.parser_args.polarisation}_" \
                           f"{chips_data.parser_args.out_label1}_" \
                           f"{chips_data.parser_args.out_label2}_comparison.png"

    elif chips_data.parser_args.polarisation == 'both':

        fig, axs = plt.subplots(2,1,figsize=(8,12))

        plot_1D_comparison(chips_data, 'xx', axs[0])
        plot_1D_comparison(chips_data, 'yy', axs[1])

        plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_xx+yy_" \
                           f"{chips_data.parser_args.out_label1}_" \
                           f"{chips_data.parser_args.out_label2}_comparison.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument\n'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)