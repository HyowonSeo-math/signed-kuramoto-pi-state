# Manuscript figures and experiment guide

This guide follows **Section IV: Numerical Simulations** of *Frequency synchronization of signed Kuramoto networks via a π-state transformation* (Choi and Seo). Manuscript figure numbers are distinguished from the numerical prefixes used in the output filenames. Figures in `figures/` are crops from the supplied manuscript PDF; they are included only for comparison and browsing.

## Section IV-A — Certified heterogeneous example

**Fig. 1(a) and Fig. 1(b); Table I.** The network has $N_p=7$ conformists, $N_m=3$ contrarians, and $K=10$. The weights are $a_{ij}=0.9+0.2\sin(ij)$ for conformist sources and $a_{ij}=0.2+0.2\sin(ij)$ for contrarian sources, with $a_{ii}=0$. Natural frequencies and transformed initial phases are sampled independently from $U(-0.5,0.5)$ and $U(-\pi/3,\pi/3)$.

The full-box certificate uses $C=2.2$, $C_1=0.75$, $\Gamma(C)\approx1.441886$, and $\Lambda(C_1)\approx2.857261$. **Fig. 1(a)** compares the lifted transformed phase diameter $D_\theta(t)$ with the scalar comparison solution $z(t)$ and the certified radii. **Fig. 1(b)** compares $D(t)$ with its exponential envelope after the observed entry time.

- `01_phase_entry.*`: Fig. 1(a).
- `02_frequency_envelope.*`: Fig. 1(b).
- `00_phase_frequency_onecolumn.*`: auxiliary, vertically stacked version of both panels.
- `baseline_certificate.csv`, `baseline_sample.mat`: numerical certificate and baseline inputs.

## Section IV-B — Behavior outside the certified region

**Fig. 2(a) and Fig. 2(b).** With the same signed baseline network, the gain and natural-frequency mismatch are varied over

$$
K\in\{0.5,1,2,3,5,10,20\},\qquad
\Delta\omega\in\{0.25,0.5,1,2,4\}.
$$

There are **20 paired trials per grid point** (Section IV-B). The numerical frequency-synchronization fraction is displayed in **Fig. 2(a)**; the sufficient-condition certification fraction is displayed in **Fig. 2(b)**. Files: `08_map_frequency.*`, `09_map_certificate.*`, and `map_trials.csv`.

**Fig. 3(a) and Fig. 3(b).** The physical initial phases are sampled independently from $U(0,2\pi)$ for **300 trials**, with $K=10$ and the baseline weights and natural frequencies. **Fig. 3(a)** displays a selected realization with increasing lifted phase diameter, together with the critical radius $C_{\mathrm{net}}\approx2.754891$ for which $\Gamma(C_{\mathrm{net}})=0$. **Fig. 3(b)** shows the corresponding frequency diameter and the $10^{-6}$ threshold. The circular-width condition for $\pi$-type grouping is checked separately from the lifted phase diameter. Files: `11_unrestricted_phase_diameter.*`, `12_unrestricted_frequency_diameter.*`, `initial_trials.csv`, `initial_trials_summary.csv`, `unrestricted_initial_illustration.mat`.

The code selects an illustrative run using a deterministic rule based on lifted-diameter growth and finite-horizon diagnostics. **This need not be the identical realization displayed in the manuscript's Fig. 3.** No unpublished exact sample is claimed.

**Fig. 4(a) and Fig. 4(b).** The natural frequencies, initial phases, and $K=10$ are retained, but source weights are changed to $0.3+0.2\sin(ij)$ for conformists and $0.6+0.2\sin(ij)$ for contrarians. The modified network violates the pairwise sufficient condition; the paper reports two positive eigenvalues of the linearization at transformed consensus. **Fig. 4(a)** compares transformed phase diameters; **Fig. 4(b)** compares frequency diameters. Files: `13_modified_network_phase.*`, `14_modified_network_frequency.*`, `modified_network_summary.csv`.

## Section IV-C — Quantitative convergence and network structure

**Table II — Coupling gain.** $K\in\{0.5,1,2,3,5,10,20,40\}$, with **100 paired trials per gain**. The table lists numerical synchronization counts, certified counts, maximum final phase diameters among synchronized trials, uniform guaranteed rates where available, and minimum terminal spectral rates. Outputs include `gain_trials.csv` and `gain_uniform_certificate.csv`. The optional plots `03_gain_dispersion.*` and `04_gain_rates.*` visualize gain dependence; they are **not manuscript Fig. 3 or Fig. 4**.

**Sparse-source construction.** The Section III design example uses a directed two-source architecture with $2(N-1)$ nonzero edges and a minority of conformists. The code also explores gain and source-weight threshold factors. Outputs: `sparse_size.csv`, `design_trials.csv`, `05_sparse_minority.*`. These plots are auxiliary to the discussion; there is no matching numbered figure in the manuscript.

**Table III — Topology and density.** ER, WS (rewiring probability 0.1), and BA undirected support graphs with generally asymmetric directed interaction weights; 70% conformists. The coupling gain is chosen to satisfy $K\langle k\rangle/N=10$, where $\langle k\rangle$ is the realized mean degree. There are **30 trials per topology, size, and density**. The corresponding outputs are `topology_trials.csv`, `topology_summary.csv`, and auxiliary figure `06_topology.*`.

**Table IV — Fixed expected degree.** ER support graphs with expected degree $10$ are evaluated at $N\in\{20,50,100,200,400\}$ for **30 paired trials per size**, comparing $K=10$ with $K=10N/\langle k\rangle$. Outputs: `fixed_degree_trials.csv`, `fixed_degree_summary.csv`, and auxiliary figure `07_fixed_degree.*`.

## Execution and generated files

```matlab
% All manuscript experiments and Tables I–IV
[results, outdir] = run_signed_kuramoto('full', [], struct('Visible','off'));

% Figs. 1–4 and Table I
[results, outdir] = run_signed_kuramoto('ab', [], struct('Visible','off'));

% Tables II–IV
[results, outdir] = run_signed_kuramoto('c', [], struct('Visible','off'));

% Only Fig. 2 and Fig. 3 experiments
[results, outdir] = run_signed_kuramoto('additional', [], struct('Visible','off'));
```

The `quick` mode uses reduced sample sizes and a smaller parameter grid; it is for preliminary runs, **not** replication of the manuscript's trial counts. Numerical figures are exported individually as `.fig`, `.pdf`, `.eps`, and `.png`, not necessarily in the side-by-side layout used in the paper.

### Other diagnostic figures

The MATLAB file also generates `00_pi_snapshots_onecolumn`, `fig1_snapshot_t0`, `fig1_snapshot_t10`, and `10_phase_omega20_blue_only`. These are additional visualization or stress-test outputs and **are not the manuscript's Figs. 1–4**.

To regenerate saved figures without new ODE runs:

```matlab
[results, outdir] = run_signed_kuramoto('replot', 'EXISTING_OUTPUT_DIRECTORY', struct('Visible','off'));
```
