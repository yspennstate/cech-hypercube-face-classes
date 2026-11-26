# Checks the face retractions R_F of the paper on small cubes.
# For every coordinate p-face F of Q_n, at r = 2p-2 (even case) and r = 2p-3 (odd case):
#  (1) R_F is simplicial N(Q_n,r) -> N(Q_p,r) on every maximal generator of Adams-Shukla-Singh Lemma 2.2,
#  (2) R_F restricted to F is the projection pi onto Q_p,
#  (3) every other p-face F' is sent into a single simplex of N(Q_p,r).
# Control: the plain clamp (pi with no correction) passes (1) and fails (3).
# Usage: python3 check_face_retractions.py [output_file]. Standard library only; under a minute.
import itertools, random, sys
pc = lambda x: bin(x).count("1")
def gens(n, r):
    V = range(1 << n); k = r // 2
    ball = lambda c: frozenset(w for w in V if pc(c ^ w) <= k)
    if r % 2 == 0: G = {ball(c) for c in V}
    else: G = {ball(u) | ball(u ^ (1 << i)) for u in V for i in range(n)}
    G = list(G); return [g for g in G if not any(g < h for h in G)]
def faces(n, p):
    for free in itertools.combinations(range(n), p):
        fixed = [i for i in range(n) if i not in free]
        for vals in itertools.product((0, 1), repeat=n - p):
            yield free, tuple(zip(fixed, vals))
def inF(x, face): return all((x >> i) & 1 == b for i, b in face[1])
def proj(x, free): return sum(((x >> i) & 1) << j for j, i in enumerate(free))
def make_R(face, mode, rng=None):
    free = face[0]
    def R(x):
        y = proj(x, free)
        if inF(x, face) or mode == "clamp": return y
        if mode == "even": return y if pc(y) % 2 == 0 else y ^ 1                    # flip first free coordinate
        if mode == "even_rand": return y if pc(y) % 2 == 0 else y ^ (1 << rng.randrange(len(free)))
        if mode == "odd":                                                          # marked bits 0,1,2 of the target word
            t = [(y >> j) & 1 for j in range(3)]
            if t[0] == t[1] == t[2]: return y
            maj = 1 if sum(t) >= 2 else 0
            j = [j for j in range(3) if t[j] != maj][0]; return y ^ (1 << j)
        raise ValueError(mode)
    return R
def into(S, T): return any(S <= g for g in T)
def run(p, n, r, mode, rng=None):
    src = gens(n, r); tgt = gens(p, r); bad1 = bad3 = 0; nf = 0
    Fs = list(faces(n, p))
    for F in Fs:
        nf += 1
        R = make_R(F, mode, rng); Rmap = {x: R(x) for x in range(1 << n)}
        assert all(Rmap[x] == proj(x, F[0]) for x in range(1 << n) if inF(x, F))
        assert {Rmap[x] for x in range(1 << n) if inF(x, F)} == set(range(1 << p))
        bad1 += sum(1 for s in src if not into(frozenset(Rmap[x] for x in s), tgt))
        for G in Fs:
            if G == F: continue
            if not into(frozenset(Rmap[x] for x in range(1 << n) if inF(x, G)), tgt): bad3 += 1
    return nf, bad1, bad3
rng = random.Random(7)
out = []
# target sanity: N(Q_p, r) not a simplex at Q_p, a simplex at Q_{p-1}
for (p, r) in [(2, 2), (3, 4), (4, 6), (3, 3), (4, 5), (5, 7)]:
    out.append("target r=%d: Q_%d full set in a generator? %s; Q_%d? %s" % (r, p, into(frozenset(range(1 << p)), gens(p, r)), p - 1, into(frozenset(range(1 << (p - 1))), gens(p - 1, r))))
cases = [("even", 2, 3), ("even", 2, 4), ("even", 2, 5), ("even", 3, 4), ("even", 3, 5), ("even", 3, 6), ("even", 4, 5), ("even", 4, 6),
         ("odd", 3, 4), ("odd", 3, 5), ("odd", 3, 6), ("odd", 4, 5), ("odd", 4, 6), ("odd", 5, 6)]
for kind, p, n in cases:
    r = 2 * p - 2 if kind == "even" else 2 * p - 3
    modes = [kind, "clamp"] + (["even_rand"] if kind == "even" else [])
    for m in modes:
        nf, b1, b3 = run(p, n, r, m, rng)
        out.append("%-5s p=%d n=%d r=%d mode=%-9s faces=%3d (expect %3d) simplicial-fails=%d other-face-fails=%d" % (kind, p, n, r, m, nf, 2 ** (n - p) * len(list(itertools.combinations(range(n), p))), b1, b3))
text = "\n".join(out) + "\n"
if len(sys.argv) > 1: open(sys.argv[1], "w", encoding="utf8").write(text)
else: print(text, end="")
