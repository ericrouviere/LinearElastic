# Methods

*Describes the model as implemented in `src/` and the simulation protocol of `exps/011` (physical parameters as in
`exps/010`). Experiments 001–010 differ in the selection weight (a sigmoid, with $\alpha$ set by the mean rather than the
spread of the random-sequence binding energies) and started every replicate in environment 1.*

## A. Elastic network model

A protein is a two-dimensional elastic network of $N = 115$ beads. The reference positions
$\mathbf{R}^0 \in \mathbb{R}^{2N}$ form a triangular lattice of 11 columns, alternating between 10 and 11 beads, with
each coordinate displaced by Gaussian noise of standard deviation $0.15$ (in units of the lattice spacing). Beads
that are nearest neighbors on the undisplaced lattice are connected by $M = 304$ springs. The geometry is drawn once
and is the same for all sequences.

A sequence $s = (s_1, \dots, s_N)$ assigns each bead one of $q = 5$ types. The spring $l = (i, j)$ has stiffness
$k_l = K_{s_i s_j}$ and rest-length mismatch $d_l = D_{s_i s_j}$, read from two symmetric $q \times q$ tables
drawn once. The entries of $\ln K$ are uniformly distributed and rescaled to span $[\ln k_{\min}, \ln k_{\max}]$
($k_{\min} = 10^{-2}$, $k_{\max} = 10$), and $D_{ab} \sim \mathcal{N}(0, \sigma_D^2)$ with $\sigma_D = 0.1$.

The elastic energy is linearized about $\mathbf{R}^0$. For a displacement $\mathbf{u} \in \mathbb{R}^{2N}$, the
elongation of spring $l$ is $\mathbf{t}_l \cdot \mathbf{u}$, where $\mathbf{t}_l$ holds the unit vector
$\hat{\mathbf{n}}_{ij}$ (pointing from $j$ to $i$ in $\mathbf{R}^0$) at bead $i$ and $-\hat{\mathbf{n}}_{ij}$ at
bead $j$. Then

$$
E_s(\mathbf{u}) = \sum_{l=1}^{M} \tfrac{1}{2} k_l \left(\mathbf{t}_l \cdot \mathbf{u} - d_l\right)^2
= \tfrac{1}{2}\mathbf{u}^\top H \mathbf{u} - \mathbf{h}^\top \mathbf{u} + c,
\tag{1}
$$

$$
H = T^\top \Gamma T, \qquad \mathbf{h} = T^\top \Gamma \mathbf{d}, \qquad c = \tfrac{1}{2}\mathbf{d}^\top \Gamma \mathbf{d},
\tag{2}
$$

where the rows of $T \in \mathbb{R}^{M \times 2N}$ are the $\mathbf{t}_l$, $\Gamma = \mathrm{diag}(k_l)$ and
$\mathbf{d} = (d_l)$. Nonzero mismatches frustrate the network and generate the prestress $\mathbf{h}$. The Hessian
$H$ has three zero modes (two translations and a rotation), which are orthogonal to $\mathbf{h}$. We write
$G = H^{+}$ for its pseudoinverse and $\mathrm{Det}'$ for the product of its nonzero eigenvalues.

The network is in thermal equilibrium at inverse temperature $\beta = 1$. Integrating over the $2N - 3$ nonzero
modes gives the free energy

$$
F(s) = c - \tfrac{1}{2}\mathbf{h}^\top G\, \mathbf{h} + \frac{1}{2\beta}\ln \mathrm{Det}'\!\left(\frac{\beta H}{2\pi}\right),
\tag{3}
$$

the energy of the ground state $\mathbf{u}^* = G\mathbf{h}$ plus the free energy of thermal fluctuations.

## B. Ligand binding

A ligand state $\lambda \in \{0, 1, 2\}$, with $\lambda = 0$ the unbound (solvent) state, changes the springs of a
binding site $\mathcal{B}$: $k_l \to k_l + \kappa_l$ and $d_l \to d_l + \delta_l^{\lambda}$ for $l \in \mathcal{B}$.
A spring of $\mathcal{B}$ that is absent from the network has $k_l = d_l = 0$. The free energy $F_\lambda(s)$ of
state $\lambda$ follows from Eq. (3) with the modified $H$, $\mathbf{h}$ and $c$. The binding free energy of
ligand $\lambda$ is

$$
\Delta F_\lambda(s) = F_\lambda(s) - F_0(s), \qquad \lambda = 1, 2,
\tag{4}
$$

and $\Delta F_\lambda < 0$ means the ligand binds.

The binding site is a single spring, with unit elongation mode $\mathbf{b}$, between two interior beads that are
next-nearest neighbors. That pair is not connected in the network, so the ligand spring does not depend on the
sequence. Every state inserts the same stiffness $\kappa = 1$, and the two ligands impose opposite strains:
$\delta^0 = 0$, $\delta^{1} = +\delta$ and $\delta^{2} = -\delta$, with $\delta = 0.2$. Ligand 1 favors opening the
site and ligand 2 favors closing it. Because all states share the Hessian $H' = H + \kappa\,\mathbf{b}\mathbf{b}^\top$,
the fluctuation term cancels in Eq. (4) and

$$
\Delta F_{1,2}(s) = \tfrac{1}{2}\kappa\left(1 - \kappa\chi\right)\delta^2 \mp \kappa\,\delta\, e_0,
\qquad e_0 = \mathbf{b}^\top G' \mathbf{h}, \quad \chi = \mathbf{b}^\top G' \mathbf{b}, \quad G' = H'^{+}.
\tag{5}
$$

The sequence therefore enters binding only through the pre-strain of the site in the unbound state, $e_0$, and its
compliance $\chi$. Since $\kappa\chi < 1$, $\Delta F_1 + \Delta F_2 > 0$: no sequence binds both ligands, and
specializing on one ligand weakens binding to the other.

## C. Sequence, mutations and phenotype

The phenotype under selection is the pair of binding energies $(\Delta F_1, \Delta F_2)$. The fitness in environment
$k \in \{1, 2\}$ is the binding score

$$
f_k(s) = -\Delta F_k(s).
\tag{6}
$$

A point mutation $s_i \to a$ changes the stiffness and mismatch of the (at most six) springs of bead $i$, a low-rank
change of $H$. We compute its effect exactly from the parent's $G$ and $\ln \mathrm{Det}' H$ with the Woodbury
identity and the matrix determinant lemma, which avoids a new diagonalization. Ligand binding is computed the same
way. The effect of a mutation on binding is

$$
\Delta\Delta F_\lambda(i, a) = \Delta F_\lambda\!\left(s^{(i \to a)}\right) - \Delta F_\lambda(s),
\tag{7}
$$

and the evolvability of $s$ is the root-mean-square effect of its $(q-1)N$ single mutants on the binding phenotype,
$\mathcal{E}(s) = \big[\frac{1}{(q-1)N}\sum_{i,\, a \neq s_i} \lVert \Delta\Delta F(i, a) \rVert^2\big]^{1/2}$.

## D. Evolutionary dynamics

A population of $P = 100$ sequences, initially independent and uniformly random, evolves by Wright–Fisher dynamics
in an environment $k(t) \in \{1, 2\}$ that alternates every $\tau$ generations. Generation $t$ consists of two steps:

1. **Selection.** The next population is drawn with replacement from the current one, choosing sequence $s$ with
   probability proportional to
   $$
   w(s) = \exp\!\left[\alpha\, f_{k(t)}(s)\right].
   \tag{8}
   $$
2. **Mutation.** Each site of each sequence mutates with probability $\mu$ to one of the other $q - 1$ types,
   chosen uniformly.

The environment switches, $k \to 3 - k$, after every generation $t$ that is a multiple of $\tau$, so it is a square
wave of period $2\tau$; the initial environment is drawn uniformly for each replicate. The selection strength is
$\alpha = a/\sigma_0$ with $a = 1$, where $\sigma_0$ is the spread of the binding energies of $10^3$ random
sequences: the standard deviation of the distance of the points $(\Delta F_1, \Delta F_2)$ from their centroid.

Each run lasts $N_g = 2048$ generations with $\mu = 3 \times 10^{-3}$ (about $0.35$ mutations per sequence per
generation) and $\tau = 2^0, 2^1, \dots, 2^8$. For each $\tau$ we evolve 100 independent replicate populations on
the same network geometry and tables. From each, we sample 50 sequences: one random individual at each of 50 random
generations in the second half of the run, within a window that spans a whole number of environmental periods.

## Parameters

| Symbol | Meaning | Value |
|---|---|---|
| $N$, $M$ | beads, springs | 115, 304 |
| $q$ | bead types | 5 |
| $k_{\min}$, $k_{\max}$ | range of spring stiffnesses | $10^{-2}$, $10$ |
| $\sigma_D$ | spread of rest-length mismatches | 0.1 |
| $\beta$ | inverse temperature | 1 |
| $\kappa$, $\delta$ | ligand spring stiffness, imposed strain | 1, 0.2 |
| $P$ | population size | 100 |
| $\mu$ | mutation rate per site per generation | $3 \times 10^{-3}$ |
| $a$ | normalized selection strength, $\alpha = a/\sigma_0$ | 1 |
| $\tau$ | environmental timescale (generations) | $1, 2, 4, \dots, 256$ |
| $N_g$ | generations per run | 2048 |
| — | replicates per $\tau$, sequences sampled per replicate | 100, 50 |
