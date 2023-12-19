#!/usr/bin/env python

import sys

from chips_wrappers.setup import make_args
from chips_wrappers.ps_methods.chips_data import ChipsDataProducts
from chips_wrappers.plotting.plot_2D import do_2D_plot, do_2D_ratio_plot, do_2D_diff_plot, do_2D_ratio_diff_plot
from chips_wrappers.plotting.plot_1D import do_1D_plot, do_1D_comparison_plot

if __name__ == '__main__':
    
    args = make_args.get_args()
    args = make_args.check_args(args)
    chips_data = ChipsDataProducts(args)
    
    if args.plot_type == '2D':
        do_2D_plot(chips_data)

    elif args.plot_type == '2D_Ratio' or args.plot_type == '2D_ratio':
        do_2D_ratio_plot(chips_data)

    elif args.plot_type == '2D_Diff' or args.plot_type == '2D_diff':
        do_2D_diff_plot(chips_data)

    elif args.plot_type == '2D_Ratio_Diff' or args.plot_type == '2D_ratio_diff':
        do_2D_ratio_diff_plot(chips_data)

    elif args.plot_type == '1D':
        do_1D_plot(chips_data)

    elif args.plot_type == '1D_comp':
        do_1D_comparison_plot(chips_data)

    else:
        sys.exit(f"You entered --plot_type={args.plot_type}, which is not"
                 "a recognised type. You can enter either: 2D, 2D_ratio,"
                 " 2D_ratio_diff, 1D, or 1D_comp")
