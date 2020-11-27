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


def stack_table(dct, axis=0, stack_func=None):
    "recursively concat dict of dict of ... of dataframes"
    stacker = stack_func or partial(pd.concat, axis=axis)
    return stacker(
        {
            k: stack_table(v)
            if isinstance(next(iter(v.values())), dict)
            else stacker(v)
            for k, v in dct.items()
        }
    )


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
