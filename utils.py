import sys

import numpy as np
from matplotlib import pyplot as plt


def plot(top, bot):
    fig, axs = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
    axs[0].scatter(top[:, 0], top[:, 1], top[:, 2], linewidth=0, antialiased=False, color="blue")
    axs[0].set_title("Top")
    axs[1].scatter(bot[:, 0], bot[:, 1], bot[:, 2], linewidth=0, antialiased=False, color="red")
    axs[1].set_title("Bottom")
    plt.show()
    return fig, axs


