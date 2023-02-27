#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Felix Höfling, Mahesh Yadav
#
import matplotlib.pyplot as plt


def get_params():
    """ """
    params = {"latex.text.preamble": [r"\usepackage{textcomp}"]}

    return params


def set_plot_defaults():
    """
    set matplotlib defaults
    """

    from pylab import rc, cycler
    from palette import LondonUnderground_Colours

    plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath} \usepackage{textcomp} \usepackage[T1]{fontenc} \usepackage{times} \usepackage{lmodern} \renewcommand{\familydefault}{\sfdefault}" 

    rc("font", **{"family": "sans-serif", "sans-serif": ["lmss"], "size": 10})
    rc("text", usetex=True)
    # rc(
    #     "text.latex",
    #     preamble=( r"\usepackage{textcomp} \usepackage{amsmath}"
    #         #r"\usepackage[T1]{fontenc}",
    #         #r"\usepackage{times}",
    #         #        r'\usepackage[scaled]{helvet} \renewcommand{\familydefault}{\sfdefault}',
    #         #r"\usepackage{lmodern} \renewcommand{\familydefault}{\sfdefault}",
    #         #        r'\usepackage[lite,eucal,subscriptcorrection,slantedGreek,zswash]{mtpro2}',
    #         #        r'\usepackage[mtphrb,mtpcal,subscriptcorrection,slantedGreek,zswash]{mtpro2}',
    #         #        r'\usepackage{txfonts}', # alternative if MathTimePro II fonts are not installed
    #         #        r'\usepackage{lmodern}',  # another alternative font family
    #     ),
    # )

    rc(
        "legend",
        frameon=False,
        loc="best",
        numpoints=1,
        fontsize="small",
        labelspacing=0.2,
        handlelength=0.5,
        handletextpad=0.5,
        borderaxespad=0.5,
    )
    rc("figure", figsize=(3.4, 2.4))  # gnuplot standard is (5, 3.5)
    rc(
        "axes",
        linewidth=0.7,
        prop_cycle=cycler(color=LondonUnderground_Colours.color_cycle),
    )
    rc("axes.formatter", use_mathtext=True)

    rc("lines", linewidth=0.8, markersize=3, markeredgewidth=0)
    rc("xtick", direction="in", top=True)
    rc("ytick", direction="in", right=True)
    rc("ytick.minor",left=True)
    rc("savefig", bbox="tight", pad_inches=0.05, dpi=1200, transparent=False)
    # http://www.scipy.org/Cookbook/Matplotlib/UsingTex
    rc("ps", usedistiller="xpdf")
