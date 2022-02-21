import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


def plot_epidemic(area_names, R, X=None, C=None, title=None):
    N = R.shape[0]

    fig, axs = plt.subplots(2, N, sharex=True, sharey="row", figsize=(10, 5))

    for i, area in enumerate(area_names):
        axs[0][i].set_ylabel("Case numbers")
        axs[1][i].plot(R[i], color="tab:red")
        axs[1][i].set_xlabel("Days into epidemic")
        axs[1][i].set_ylabel("R")
        axs[0][i].set_title(area)

    if X is not None:
        for i, area in enumerate(area_names):
            axs[0][i].plot(X[i], color="tab:blue", alpha=0.6)

    if C is not None:
        for i, area in enumerate(area_names):
            axs[0][i].plot(C[i], color="tab:orange", alpha=0.6)

    for ax in axs.flat:
        ax.label_outer()

    import matplotlib.patches as mpatches

    handles = [
        mpatches.Patch(color=c, label=k)
        for (k, c) in {
            "Latent Epidemic": "tab:blue",
            "Observed Cases": "tab:orange",
        }.items()
    ]
    axs[0][N - 1].legend(handles=handles, loc="upper left", bbox_to_anchor=(1.05, 1))

    if title is not None:
        plt.suptitle(title)
