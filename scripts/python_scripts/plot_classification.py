import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

classification_dir = "classfication_SMP_Repository/"
classification_plots = "classification_plots/"

def plot_classification(classification_dir):
    """
    plot the distribution of metrics in a plot, to show what value of metric
    can lead to
        1. always fail
        2. sometimes fail
        3. never fail
    """
    files = os.listdir(classification_dir)
    if not os.path.exists(classification_plots):
        os.makedirs(classification_plots)

    for f in files:
        if f.endswith(".csv"):
            plt.figure()
            fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(7, 7))
            fig.suptitle(f[:-4])

            data = pd.read_csv(classification_dir + f)
            norms = data["norm"].values
            num_constrs = data["number of constraints"].values
            num_ranges = data["numerical range"].values

            axs[0].plot(np.arange(1, len(data) + 1), np.sort(norms))
            axs[0].set(xlabel="count", ylabel="norm")
            plt.sca(axs[0]); plt.yscale("log")

            axs[1].plot(np.arange(1, len(data) + 1), np.sort(num_constrs))
            axs[1].set(xlabel="count", ylabel="number of constraints")
            plt.sca(axs[1]); plt.yscale("log")

            axs[2].plot(np.arange(1, len(data) + 1), np.sort(num_ranges))
            axs[2].set(xlabel="count", ylabel="numerical range")
            plt.sca(axs[2]); plt.yscale("log")

            plt.tight_layout()
            plt.savefig(classification_plots + f[:-4] + ".png")
    

if __name__ == "__main__":
    plot_classification(classification_dir)