import os
from matplotlib import pyplot as plt


def setPlotSettings(font=False):

    ''' Set font and dpi for plotting.
    
    Inputs:
    font = Boolean, whether to set font to LaTeX -- requires LaTeX installation! '''
    
    # Set dpi:
    plt.rcParams['figure.dpi'] = 300

    # Set LaTeX font:
    if font == True:
        os.environ["PATH"] = "/Library/TeX/texbin:" + os.environ["PATH"] # path to LaTeX installation
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern"],
            "text.latex.preamble": r"\usepackage{amsmath}",
            "font.size": 14
        })