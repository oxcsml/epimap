from functools import partial
import os

from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


def map_lowest(func, dct):
    return {
        k: map_lowest(func, v) if isinstance(v, dict) else func(v)
        for k, v in dct.items()
    }


def swaplevel(dct_of_dct):
    keys = next(iter(dct_of_dct.values())).keys()
    return {in_k: {out_k: v[in_k] for out_k, v in dct_of_dct.items()} for in_k in keys}


def collapse(dct, axis=0, names=None):
    """collapse.
    Collapse nested dictionary of dataframes into MultiIndex DataFrame

    Args:
        dct (dict): nested dictionary of dataframes.
        axis (int): axis on which to concat, default 0.
        names (list, optional): index names to pass to pd.concat, default None.
    """
    if isinstance(next(iter(dct.values())), dict):
        return collapse({k: collapse(v, axis=axis, names=None) for k, v in dct.items()}, axis=axis, names=names)
    else:
        return pd.concat(dct, axis=axis, names=names)
    

class PdfDeck:
    def __init__(self, figs=None, names=None):
        self.figs = figs or []
        self.fignames = names or []

    @classmethod
    def save_as_pdf(cls, figs, fpath):
        return cls(figs).make(fpath)

    def default_figname(self):
        return f"{repr(self).replace(' ', '_')}_figure_{len(self.figs)}"

    def add_figure(self, fig, *, position=None, name=None):
        if position is None:
            self.figs.append(fig)
            self.fignames.append(name or self.default_figname())
        else:
            self.figs.insert(position, fig)

    def make(self, fpath):
        with PdfPages(fpath) as pdf:
            for fig in self.figs:
                pdf.savefig(fig)

    def make_individual(self, folder=None, **savefig_kwds):
        folder = folder or os.cwd()
        for fig, name in zip(self.figs, self.fignames):
            fpath = os.path.join(folder, name + "." + savefig_kwds.get("format", "pdf"))
            fig.savefig(fpath, **savefig_kwds)
