# Independent face classes in Čech complexes of hypercube graphs

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

Let N(G,r) be the Čech complex of a graph G at scale r. For the hypercube graph Q_n and the scales r = 2p-2 with p >= 2 and r = 2p-3 with p >= 3, at which Q_p is the smallest cube whose complex is not a simplex, the 2^(n-p) C(n,p) coordinate p-faces of Q_n carry independent split copies of the reduced homology of N(Q_p,r), with any coefficients. The retraction onto a face is the coordinate clamp of Adams, Shukla and Singh (Čech complexes of hypercube graphs, Homology Homotopy Appl. 27 (2025) 83-105, doi:10.4310/HHA.2025.v27.n1.a6, arXiv:2212.05871), corrected by one bit so that it sends every other face into a simplex. Consequently the homology of N(Q_n,2p-2) has rank at least 2^(n-p) C(n,p) in degree 2^p-2, and H_4(N(Q_n,3)) has rank at least 3·2^(n-3) C(n,3), which equals the computed values 24 and 120 for n = 4, 5. With the free summand of rank at least 3 in degree 3·2^(p-2)-2 of N(Q_p,2p-3) for p >= 3 (Collapsibility numbers and high-degree homology of Čech complexes of hypercube graphs, doi:10.5281/zenodo.24569831, Theorem 1.2), the same degree of N(Q_n,2p-3) has a free summand of rank at least 3·2^(n-p) C(n,p) for n >= p.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/check_face_retractions.py`: for every coordinate p-face of Q_n in fourteen cases (p,n) with n <= 6, checks that the corrected retraction is simplicial on the maximal generators of Adams-Shukla-Singh Lemma 2.2, restricts to the projection on its own face, and sends every other face into one simplex; the plain clamp is run as a control and fails the last test. At even scales it also runs `even_rand`, which flips a randomly chosen coordinate of each odd word instead of the first one (the proof allows any single coordinate); the first lines of the output check that N(Q_p,r) is not a simplex while N(Q_(p-1),r) is. Run `python3 check_face_retractions.py`; standard library only, under a minute.
- `check_output.txt`: the output of that run.
- `SHA256SUMS`: sha256 of each file above.
