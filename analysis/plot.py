import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(name):
    with open(f"results/{name}.csv", newline="") as fh:
        return [{k: try_num(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def try_num(v):
    try:
        return float(v)
    except ValueError:
        return v


def group(rows, key):
    out = defaultdict(list)
    for r in rows:
        out[r[key]].append(r)
    return out


def count_min_bound():
    fig, ax = plt.subplots(figsize=(6, 5))
    for stream, rows in group(load("cm_bound"), "stream").items():
        rows.sort(key=lambda r: r["bound"])
        b = [r["bound"] for r in rows]
        # exact estimates are zero error, which a log axis cannot show
        ax.plot(b, [r["max_err"] or None for r in rows], "o-",
                label=f"{stream} (max)")
        ax.plot(b, [r["mean_err"] or None for r in rows], "o--", alpha=.5,
                label=f"{stream} (mean)")
    lim = [min(r["bound"] for r in load("cm_bound")),
           max(r["bound"] for r in load("cm_bound"))]
    ax.plot(lim, lim, "k:", label="guarantee: eps*N")
    ax.set(xscale="log", yscale="log", xlabel="bound eps*N",
           ylabel="observed overestimate",
           title="Count-Min: observed error vs guarantee")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig("figures/cm_bound.png", dpi=150)


def count_min_hash():
    rows = load("cm_hash")
    streams = list(dict.fromkeys(r["stream"] for r in rows))
    by_hash = group(rows, "hash")
    x = range(len(streams))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    width = 0.38
    floor = 3e-2
    for i, (hash_name, hrows) in enumerate(by_hash.items()):
        lookup = {r["stream"]: r["max_over_bound"] for r in hrows}
        vals = [lookup[s] for s in streams]
        pos = [p + i * width for p in x]
        ax.bar(pos, [v if v > 0 else floor for v in vals], width, label=hash_name)
        for px, v in zip(pos, vals):
            if v == 0:
                ax.text(px, floor * 1.1, "exact", ha="center", va="bottom",
                        fontsize=6, rotation=90)
    ax.axhline(1.0, color="k", ls=":", label="guarantee")
    ax.set(yscale="log", ylim=(floor * 0.7, None), ylabel="max error / (eps*N)",
           title="Count-Min: cost of a non-2-universal hash")
    ax.set_xticks([p + width / 2 for p in x])
    ax.set_xticklabels(streams)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/cm_hash.png", dpi=150)


def hyperloglog():
    rows = load("hll")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    big = [r for r in rows if r["n_true"] == 1000000]
    for keys, krows in group(big, "keys").items():
        krows.sort(key=lambda r: r["m"])
        axes[0].plot([r["m"] for r in krows], [r["rmse_rel"] for r in krows],
                     "o-", label=f"{keys} keys")
    ref = sorted({r["m"]: r["predicted_se"] for r in rows}.items())
    axes[0].plot([m for m, _ in ref], [s for _, s in ref], "k:",
                 label="1.04/sqrt(m)")
    axes[0].set(xscale="log", yscale="log", xlabel="registers m",
                ylabel="relative RMSE",
                title="Error vs prediction (n = 1e6)")
    axes[0].legend(fontsize=8)

    p14 = sorted([r for r in rows if r["p"] == 14 and r["keys"] == "random"],
                 key=lambda r: r["n_true"])
    axes[1].plot([r["n_true"] for r in p14], [r["rmse_rel"] for r in p14], "o-",
                 label="corrected")
    axes[1].plot([r["n_true"] for r in p14], [r["rmse_rel_raw"] for r in p14],
                 "s--", label="raw estimator")
    axes[1].axhline(p14[0]["predicted_se"], color="k", ls=":",
                    label="1.04/sqrt(m)")
    axes[1].set(xscale="log", yscale="log", xlabel="true cardinality",
                ylabel="relative RMSE",
                title="What the small-range correction buys (p = 14)")
    axes[1].legend(fontsize=8)

    for keys, krows in group([r for r in rows if r["p"] >= 8], "keys").items():
        pts = sorted((r["n_true"] / r["m"], r["mean_rel_bias"]) for r in krows)
        axes[2].plot([x for x, _ in pts], [y for _, y in pts], "o",
                     alpha=.7, label=f"{keys} keys")
    axes[2].axhline(0, color="k", ls=":")
    axes[2].set(xscale="log", xlabel="load n/m", ylabel="mean relative bias",
                title="Random-oracle assumption: bias vs load")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig("figures/hll.png", dpi=150)


def reservoir():
    rows = load("reservoir_positions")
    p = rows[0]["expected_p"]
    trials = 20000
    sd = (p * (1 - p) / trials) ** 0.5
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([r["position"] for r in rows], [r["inclusion_freq"] for r in rows],
            lw=.7, color="steelblue")
    ax.axhline(p, color="k", label=f"k/n = {p}")
    ax.axhspan(p - 3 * sd, p + 3 * sd, color="k", alpha=.15, label="+/- 3 sigma")
    ax.set(xlabel="stream position", ylabel="inclusion frequency",
           title=f"Reservoir sampling: no positional bias (n=1000, k=50, {trials} trials)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/reservoir.png", dpi=150)


def misra_gries():
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for stream, rows in group(load("misra_gries"), "stream").items():
        rows.sort(key=lambda r: r["k"])
        ax.plot([r["k"] for r in rows], [r["max_underest"] for r in rows], "o-",
                label=stream)
    rows = sorted(load("misra_gries"), key=lambda r: r["k"])
    ax.plot([r["k"] for r in rows], [r["bound"] for r in rows], "k:",
            label="bound N/(k+1)")
    ax.set(xscale="log", yscale="log", xlabel="counters k",
           ylabel="max underestimate",
           title="Misra-Gries: deterministic bound vs observed")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/misra_gries.png", dpi=150)


if __name__ == "__main__":
    count_min_bound()
    count_min_hash()
    hyperloglog()
    reservoir()
    misra_gries()
    print("figures/*.png written")
