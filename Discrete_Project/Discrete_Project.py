"""
CS1005 Discrete Structures - Fall 2026 - Network motifs in email-Eu-core

Group roll numbers (last 4 digits): 0775, 1050, 1052
R = 775 + 1050 + 1052 = 2877
  pair number       = R mod 190 = 27  -> (i, j) = (1, 10)
  random networks   = 20 + R mod 31 = 45
  swap attempts     = (5 + R mod 6) * |E| = 8 * |E|
  random seed       = 2877

Run:  python motifs.py            (reads the SNAP files from the data/ folder)
      python motifs.py --data-dir path/to/files
"""
import argparse
import gzip
import itertools
import math
import os
import random

ROLLS = [775, 1050, 1052]
R = sum(ROLLS)
N_RANDOM = 20 + (R % 31)
SWAP_MULT = 5 + (R % 6)
SEED = R

EDGE_NAME = "email-Eu-core.txt.gz"
DEPT_NAME = "email-Eu-core-department-labels.txt.gz"

# ---------------------------------------------------------------- T3
# bit k of the code <-> k-th edge below
EDGE_BITS = {(0, 1): 0, (0, 2): 1, (1, 0): 2, (1, 2): 3, (2, 0): 4, (2, 1): 5}
CLASS_CODES = [3, 6, 7, 10, 11, 15, 21, 23, 25, 27, 30, 31, 63]
CODE_TO_CLASS = {c: k + 1 for k, c in enumerate(CLASS_CODES)}


def canonical_code(has_edge, nodes):
    """Minimum 6-bit code over all 6 orderings of the 3 nodes.
    has_edge(u, v) -> bool tells whether u->v is present."""
    best = 64
    for order in itertools.permutations(nodes):
        code = 0
        for (i, j), bit in EDGE_BITS.items():
            if has_edge(order[i], order[j]):
                code |= 1 << bit
        if code < best:
            best = code
    return best


def triad_class(a, b, c, out):
    """Class 1..13 of the triad {a,b,c}; `out` maps node -> set of out-neighbours.
    Returns None if the triad is not connected."""
    code = canonical_code(lambda u, v: v in out[u], (a, b, c))
    return CODE_TO_CLASS.get(code)


# ---------------------------------------------------------------- T1
def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def find_file(data_dir, name):
    """Look only in the data folder (next to this script, or ./data in the
    folder you run from). Accepts the .txt.gz or the plain .txt version.
    Nothing is downloaded."""
    here = os.path.dirname(os.path.abspath(__file__))
    folders = [data_dir]
    if not os.path.isabs(data_dir):
        folders.append(os.path.join(here, data_dir))
    for folder in folders:
        for fname in (name, name[:-3]):
            p = os.path.join(folder, fname)
            if os.path.exists(p):
                print("using local file:", p)
                return p
    raise FileNotFoundError(
        f"Could not find '{name}' (or '{name[:-3]}') in the '{data_dir}' folder. "
        f"Put the SNAP files in a folder called '{data_dir}' next to motifs.py.")


def load(data_dir):
    edges = set()
    with _open(find_file(data_dir, EDGE_NAME)) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2 and p[0] != p[1]:          # drop self-loops
                edges.add((int(p[0]), int(p[1])))      # set drops duplicates
    dept = {}
    with _open(find_file(data_dir, DEPT_NAME)) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                dept[int(p[0])] = int(p[1])
    return sorted(edges), dept


def top20(edges, dept):
    internal = {}
    for u, v in edges:
        if u in dept and v in dept and dept[u] == dept[v]:
            internal[dept[u]] = internal.get(dept[u], 0) + 1
    ranked = sorted(internal, key=lambda d: (-internal[d], d))
    return ranked[:20], internal


def pair_from_number(n):
    pairs = [(i, j) for i in range(20) for j in range(i + 1, 20)]
    return pairs[n]


# ---------------------------------------------------------------- T2
def subgraph(edges, dept, d):
    members = {n for n, x in dept.items() if x == d}
    return members, [(u, v) for u, v in edges if u in members and v in members]


# ---------------------------------------------------------------- T4
# Precomputed table: 64 labelled configurations -> class (None if disconnected)
_CODE_CLASS = {}
for _cfg in range(64):
    best = 64
    for perm in itertools.permutations(range(3)):
        code = 0
        for (i, j), bit in EDGE_BITS.items():
            if _cfg >> EDGE_BITS[(perm[i], perm[j])] & 1:
                code |= 1 << bit
        best = min(best, code)
    _CODE_CLASS[_cfg] = CODE_TO_CLASS.get(best)


def count_triads(edges):
    """Triad census N[1..13] (list index 0..12). Each connected triple is counted
    exactly once: it is counted at the smallest-id node among the nodes that are
    adjacent to both others (its 'centres'), so no set of triples is needed."""
    out, nbr = {}, {}
    for u, v in edges:
        out.setdefault(u, set()).add(v)
        nbr.setdefault(u, set()).add(v)
        nbr.setdefault(v, set()).add(u)
        out.setdefault(v, set())
    N = [0] * 13
    for v, nb in nbr.items():
        for a, b in itertools.combinations(sorted(nb), 2):
            # a and b adjacent -> all three are centres; count only at the minimum
            if b in nbr[a] and not (v < a):
                continue
            t = (v, a, b)
            cfg = 0
            for (i, j), bit in EDGE_BITS.items():
                if t[j] in out[t[i]]:
                    cfg |= 1 << bit
            N[_CODE_CLASS[cfg] - 1] += 1
    return N


# ---------------------------------------------------------------- T5
def randomize(edges, attempts, seed):
    """Double-edge swaps. `seed` may be an int or a random.Random instance.
    Rejected attempts still count. Always starts from the given edge list."""
    rng = seed if isinstance(seed, random.Random) else random.Random(seed)
    E = list(edges)
    S = set(E)
    m = len(E)
    for _ in range(attempts):
        i = rng.randrange(m)
        j = rng.randrange(m)
        if i == j:                       # same edge picked twice -> reject
            continue
        a, b = E[i]
        c, d = E[j]
        if a == d or c == b or (a, d) in S or (c, b) in S:
            continue
        S.discard((a, b)); S.discard((c, d))
        E[i], E[j] = (a, d), (c, b)
        S.add((a, d)); S.add((c, b))
    return E


def degrees(edges):
    ind, outd = {}, {}
    for u, v in edges:
        outd[u] = outd.get(u, 0) + 1
        ind[v] = ind.get(v, 0) + 1
    return ind, outd


# ---------------------------------------------------------------- T6
def zscores(real, samples):
    z = []
    n = len(samples)
    for k in range(13):
        vals = [s[k] for s in samples]
        mu = sum(vals) / n
        sd = math.sqrt(sum((x - mu) ** 2 for x in vals) / n)   # population std
        z.append(0.0 if sd == 0 else (real[k] - mu) / sd)
    return z


def analyse(name, edges, rng):
    real = count_triads(edges)
    ind0, out0 = degrees(edges)
    attempts = SWAP_MULT * len(edges)
    samples = []
    for _ in range(N_RANDOM):
        g = randomize(edges, attempts, rng)
        ind, outd = degrees(g)
        assert ind == ind0 and outd == out0, "degree sequence changed!"
        assert len(set(g)) == len(g) and all(u != v for u, v in g), "not simple"
        samples.append(count_triads(g))
    print(f"{name}: degree-preservation assertions passed for {N_RANDOM} graphs")
    return real, zscores(real, samples)


# ---------------------------------------------------------------- T7
def plot(zA, zB, labelA, labelB, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = list(range(1, 14))
    w = 0.38
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar([i - w / 2 for i in x], zA, w, label=labelA)
    ax.bar([i + w / 2 for i in x], zB, w, label=labelB)
    for y in (2, -2):
        ax.axhline(y, ls="--", color="k", lw=1)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(x)
    ax.set_xlabel("Triad class")
    ax.set_ylabel("Z-score")
    ax.set_title(f"Motif Z-scores, departments {labelA} vs {labelB}  (R = {R})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    # ---- T1
    edges, dept = load(args.data_dir)
    nodes = {x for e in edges for x in e}
    print(f"nodes in edges: {len(nodes)}, simple edges: {len(edges)}")
    top, internal = top20(edges, dept)
    print("top-20 departments (index: dept id, internal edges):")
    for k, d in enumerate(top):
        print(f"  {k:2d}: dept {d:3d}  {internal[d]}")
    pair = R % 190
    i, j = pair_from_number(pair)
    A, B = top[i], top[j]
    print(f"R = {R}; pair number = {pair} -> ({i}, {j}); "
          f"Department A = {A}, Department B = {B}")
    print(f"random networks = {N_RANDOM}; swap attempts = {SWAP_MULT} x |E|; seed = {SEED}")

    # ---- T2
    results = {}
    for name, d in (("A", A), ("B", B)):
        members, sub = subgraph(edges, dept, d)
        used = {x for e in sub for x in e}
        print(f"Dept {name} ({d}): {len(members)} members, "
              f"{len(used)} with at least one internal edge, {len(sub)} edges")
        results[name] = sub

    # ---- T3 hand example (cycle a->b->c->a): code 25 -> class 9
    demo = {0: {1}, 1: {2}, 2: {0}}
    print("hand example: 0->1,1->2,2->0 gives class", triad_class(0, 1, 2, demo),
          "(expected 9, code 25)")

    # ---- T4 / T6
    rng = random.Random(SEED)
    final = {}
    for name in ("A", "B"):
        final[name] = analyse(name, results[name], rng)

    print("\nclass | N_A | Z_A | N_B | Z_B")
    for k in range(13):
        print(f"{k+1:5d} | {final['A'][0][k]:5d} | {final['A'][1][k]:7.2f} | "
              f"{final['B'][0][k]:5d} | {final['B'][1][k]:7.2f}")
    for name in ("A", "B"):
        mot = [k + 1 for k, z in enumerate(final[name][1]) if z > 2]
        anti = [k + 1 for k, z in enumerate(final[name][1]) if z < -2]
        print(f"Dept {name}: motifs {mot}, anti-motifs {anti}")

    os.makedirs(args.out_dir, exist_ok=True)
    plot(final["A"][1], final["B"][1], f"A (dept {A})", f"B (dept {B})",
         os.path.join(args.out_dir, f"zscores_R{R}.png"))

    # optional sanity check against networkx (different class order!) - checking only
    try:
        import networkx as nx
        for name in ("A", "B"):
            G = nx.DiGraph(results[name])
            tc = nx.triadic_census(G)
            print(f"[check] Dept {name}: networkx connected-triad total =",
                  sum(v for k, v in tc.items() if k not in ("003", "012", "102")),
                  "vs mine", sum(final[name][0]))
    except ImportError:
        pass


if __name__ == "__main__":
    main()