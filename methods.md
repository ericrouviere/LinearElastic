# Methods

*Describes the model as implemented in `src/` and the simulation protocol of `exps/011` (physical parameters as in `exps/010`). Experiments 001–010 differ in the selection weight (a sigmoid, with the binding-energy scale set by the mean rather than the spread of the random-sequence binding energies) and started every replicate in environment 1.*

## Sequence

The model represents a protein as an elastic network of $N = 115$ nodes, standing for amino acid residues, at fixed reference positions $\mathbf{R}^0 \in \mathbb{R}^{2N}$ in two dimensions. The nodes form a triangular lattice of 11 columns that alternate between 10 and 11 nodes, and each coordinate is displaced by Gaussian noise of standard deviation 0.15 (in units of the lattice spacing) to break the symmetry of the lattice. This gives the model a fixed structure and thus does not allow for unfolding. Each node carries one of $q = 5$ amino acid types, $s_i \in \{1, \dots, q\}$, and the full sequence is $\mathbf{s} = (s_1, \dots, s_N)$.

## Interactions

Nodes $i$ and $j$ that are nearest neighbors on the undisplaced lattice are connected by a spring; we write $\langle ij \rangle$ for the $M = 304$ edges of the network. The spring between $i$ and $j$ has a stiffness $k_{ij}(s_i, s_j)$, a rest-length mismatch $d_{ij}(s_i, s_j)$ and a contact energy $\varepsilon_{ij}(s_i, s_j)$ that depend on the types of the two amino acids, read from three symmetric $q \times q$ tables. The values of $\ln k$ are drawn uniformly and rescaled to span $[\ln k_{\min}, \ln k_{\max}]$, with $k_{\min} = 10^{-2}$ and $k_{\max} = 10$, and the mismatches and contact energies are drawn from normal distributions with mean zero and standard deviations $\sigma_D = 0.1$ and $\sigma_\varepsilon$, respectively ($\sigma_\varepsilon = 0$ in the simulations reported here). All three tables are fixed during evolution. A mismatch $d_{ij} > 0$ means that the spring prefers to be longer than in $\mathbf{R}^0$, so nonzero mismatches frustrate the network and create internal stress. The contact energy is a direct interaction between the two amino acids that does not depend on the conformation; it contributes to the stability of a sequence but not to its mechanics.

## Energy landscape

We treat the network in linear elasticity. For a displacement $\mathbf{u} \in \mathbb{R}^{2N}$ of the nodes from $\mathbf{R}^0$, the extension of spring $\langle ij \rangle$ is $\mathbf{t}_{ij} \cdot \mathbf{u}$, where $\mathbf{t}_{ij} \in \mathbb{R}^{2N}$ holds the unit vector $\hat{\mathbf{n}}_{ij}$ pointing from $j$ to $i$ at node $i$, $-\hat{\mathbf{n}}_{ij}$ at node $j$, and zeros elsewhere. Ligand binding is modeled by an additional spring of stiffness $\kappa$ between two active-site nodes $b_1$ and $b_2$, next-nearest neighbors in the interior of the network that are not otherwise connected, whose rest-length offset $\lambda$ is set by the ligand. The energy of the model is

$$U(\mathbf{u}, \mathbf{s}, \lambda) = \sum_{\langle ij \rangle} \left[ \tfrac{1}{2} k_{ij}(s_i, s_j) \left(\mathbf{t}_{ij} \cdot \mathbf{u} - d_{ij}(s_i, s_j)\right)^2 + \varepsilon_{ij}(s_i, s_j) \right] + \tfrac{1}{2} \kappa \left(\mathbf{b} \cdot \mathbf{u} - \lambda\right)^2, \tag{1}$$

with $\mathbf{b} = \mathbf{t}_{b_1 b_2}$. The last term implements ligand binding. Different ligands are represented by different values of $\lambda$: $\lambda_1$ and $\lambda_2$ characterize the effects of two different ligands, and $\lambda_0$ represents the effect of the solvent interacting with the active site. To impose a binding trade-off, we choose ligand values that are opposite relative to the solvent condition, $\lambda_0 = 0$, $\lambda_1 = \delta$ and $\lambda_2 = -\delta$, with $\delta = 0.2$ and $\kappa = 1$: ligand 1 favors opening the active site and ligand 2 favors closing it. Expanding Eq. (1) gives the quadratic energy landscape

$$U(\mathbf{u}, \mathbf{s}, \lambda) = \tfrac{1}{2} \mathbf{u}^\top H(\mathbf{s})\, \mathbf{u} - \mathbf{h}(\mathbf{s}, \lambda)^\top \mathbf{u} + c(\mathbf{s}, \lambda), \tag{2}$$

$$H = \sum_{\langle ij \rangle} k_{ij} \mathbf{t}_{ij} \mathbf{t}_{ij}^\top + \kappa \mathbf{b} \mathbf{b}^\top, \qquad \mathbf{h} = \sum_{\langle ij \rangle} k_{ij} d_{ij} \mathbf{t}_{ij} + \kappa \lambda \mathbf{b}, \qquad c = \sum_{\langle ij \rangle} \left( \tfrac{1}{2} k_{ij} d_{ij}^2 + \varepsilon_{ij} \right) + \tfrac{1}{2} \kappa \lambda^2. \tag{3}$$

The Hessian $H$ depends on the sequence but not on the ligand. It has three zero modes, two translations and a rotation, which are orthogonal to $\mathbf{h}$. We write $C = H^{+}$ for its pseudoinverse, the compliance matrix, and $\mathrm{Det}'$ for the product of its nonzero eigenvalues.

## Free energy

At inverse temperature $\beta$, integrating the Boltzmann weight over the $2N - 3$ nonzero modes gives the free energy of the model as a function of sequence and ligand,

$$F(\mathbf{s}, \lambda) = -\frac{1}{\beta} \ln \int d\mathbf{u}\, e^{-\beta U(\mathbf{u}, \mathbf{s}, \lambda)} = c - \tfrac{1}{2} \mathbf{h}^\top C\, \mathbf{h} + \frac{1}{2\beta} \ln \mathrm{Det}'\!\left(\frac{\beta H}{2\pi}\right). \tag{4}$$

The first two terms are the energy of the ground state $\mathbf{u}^* = C \mathbf{h}$, and the last term is the free energy of thermal fluctuations about it. We use $\beta = 1$.

The network can also unfold. Folding is a hard stability threshold $F_{\mathrm{u}}$ on the free energy of the network in solvent: a sequence folds if $F(\mathbf{s}, \lambda_0) \le F_{\mathrm{u}}$ and is unfolded otherwise,

$$\phi(\mathbf{s}) = \begin{cases} 1 & \text{if } F(\mathbf{s}, \lambda_0) \le F_{\mathrm{u}}, \\ 0 & \text{otherwise}. \end{cases} \tag{5}$$

$F_{\mathrm{u}}$ is a parameter of the model; $F_{\mathrm{u}} = \infty$ means that no sequence unfolds.

## Binding free energies

The binding free energy is the difference in free energy between ligand-bound and solvent-bound conditions,

$$\Delta F_k(\mathbf{s}) = F(\mathbf{s}, \lambda_k) - F(\mathbf{s}, \lambda_0), \qquad k = 1, 2. \tag{6}$$

Because the ligand changes $\mathbf{h}$ and $c$ but not $H$ or the sequence, the fluctuation term and the contact energies cancel, and Eq. (4) gives

$$\Delta F_k(\mathbf{s}) = \tfrac{1}{2} \kappa \left(1 - \kappa \chi\right) \lambda_k^2 - \kappa \lambda_k e_0, \qquad e_0 = \mathbf{b}^\top C\, \mathbf{h}(\mathbf{s}, \lambda_0), \quad \chi = \mathbf{b}^\top C\, \mathbf{b}. \tag{7}$$

The binding energies $\Delta F_1$ and $\Delta F_2$ depend only on the sequence, yielding two sequence-to-function maps that can be used to evolve a population. The sequence enters them through two properties of the active site in the solvent-bound state: its pre-strain $e_0$, the extension of the active site in the ground state, and its compliance $\chi$. A pre-strain toward opening ($e_0 > 0$) favors ligand 1 at the expense of ligand 2. Since $\kappa \chi < 1$, $\Delta F_1 + \Delta F_2 = \kappa (1 - \kappa \chi)\, \delta^2 > 0$, so no sequence binds both ligands and tighter binding of one comes at the expense of the other. Sequences that strike an optimal balance between the two, by maximizing the weaker of the two affinities $G = \min(-\Delta F_1, -\Delta F_2)$, define generalists.

## Selection

Selection rewards tighter binding, and folding, with increased probability of reproduction. With $\Delta F^n_k$ denoting the binding energy to ligand $k \in \{1, 2\}$ of the $n$th individual. To implement a selection for stability, a non functional unfolded state is defined to have a free energy $F_{\mathrm{unfolded}}$ . If the free energy in the unbound state, $F(\mathbf{s}, \lambda_0)$ is greater than $F_{\mathrm{unfolded}}$ then the system is taken to be unfolded and is does not reproduce. This is captured by variable  $\phi^n =  \Theta(F_{\mathrm{unfolded}} -  F(\mathbf{s}, \lambda_0)$), where $\Theta$ is the Heaviside function. The population at generation $t + 1$ is obtained by sampling $P$ individuals at generation $t$ according to a multinomial distribution with probabilities

$$p^n_t = \frac{\phi^n\, e^{-a \Delta F^n_{k(t)} / \Delta F_0}}{\sum_{n'=1}^{P} \phi^{n'}\, e^{-a \Delta F^{n'}_{k(t)} / \Delta F_0}}, \qquad k(t) = 1 + \left[\left(\left\lfloor \frac{t - 1}{\tau} \right\rfloor + k_0 - 1\right) \bmod 2\right]. \tag{8}$$

A sequence that does not fold therefore has zero probability of reproduction, whatever its binding energy. The switching selective pressure is implemented by $k(t)$, which alternates between $k = 1$ and $k = 2$ every $\tau$ generations, starting from an environment $k_0$ drawn uniformly for each replicate. Under static selection, $k(t) = k$ is held fixed. $\Delta F_0$ is the scale of the binding energies. Unlike in a model with bounded energies, it has no simple closed form here, so we take it as the spread of the binding energies of $10^3$ uniformly random sequences: the standard deviation of the distance of the points $(\Delta F_1, \Delta F_2)$ from their centroid. The parameter $a > 0$ defines the selective pressure, and we use $a = 1$.

## Evolutionary dynamics

We simulated the evolution of a population of fixed size $P = 100$ using Wright–Fisher dynamics, with selection switching between acting on $\Delta F_1$ and $\Delta F_2$ every $\tau$ generations. Populations were initialized with independent, uniformly random sequences. Each generation consists of selection, Eq. (8), followed by mutation: each amino acid mutates with probability $\mu$ to one of the other $q - 1$ types, chosen uniformly. We used $\mu = 3 \times 10^{-3}$ (about 0.35 mutations per sequence per generation) and $\tau = 1, 2, 4, \dots, 256$, and ran each trajectory for $N_g = 2048$ generations. For each $\tau$, 100 replicate populations are evolved on the same network and interaction tables, and 50 sequences are sampled from each: one random individual at each of 50 random generations in the last half of the trajectory, within a window that spans a whole number of environmental periods.

## Parameters

| Symbol | Meaning | Value |
|---|---|---|
| $N$, $M$ | nodes, springs | 115, 304 |
| $q$ | amino acid types | 5 |
| $k_{\min}$, $k_{\max}$ | range of spring stiffnesses | $10^{-2}$, $10$ |
| $\sigma_D$ | spread of rest-length mismatches | 0.1 |
| $\sigma_\varepsilon$ | spread of contact energies | 0 |
| $\kappa$ | ligand spring stiffness | 1 |
| $\lambda_0, \lambda_1, \lambda_2$ | solvent and ligand rest-length offsets | $0, 0.2, -0.2$ |
| $\beta$ | inverse temperature | 1 |
| $F_{\mathrm{u}}$ | stability threshold on $F(\mathbf{s}, \lambda_0)$ | $\infty$ (no unfolding) |
| $P$ | population size | 100 |
| $\mu$ | mutation rate per site per generation | $3 \times 10^{-3}$ |
| $a$ | selective pressure | 1 |
| $\tau$ | environmental timescale (generations) | $1, 2, 4, \dots, 256$ |
| $N_g$ | generations per trajectory | 2048 |
| — | replicates per $\tau$, sequences sampled per replicate | 100, 50 |
