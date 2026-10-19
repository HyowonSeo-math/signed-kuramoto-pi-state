Frequency Synchronization of Signed Kuramoto Networks via a π-State Transformation
MATLAB code for the numerical experiments in “Frequency synchronization of signed Kuramoto networks via a π-state transformation” by Sun-Ho Choi and Hyowon Seo.
The simulations concern directed, weighted conformist–contrarian Kuramoto networks. They test a sufficient condition for phase cohesion and exponential frequency synchronization, examine its conservativeness, and compare coupling gain, source weights, graph topology, and network size. The experiment descriptions and figure numbers below follow Section IV of the manuscript.
Numerical results in the paper
Manuscript	Experiment	Simulation output
Fig. 1(a)	Certified transformed phase-diameter bound and scalar comparison equation	`01_phase_entry.*`
Fig. 1(b)	Frequency-diameter decay and theoretical exponential envelope	`02_frequency_envelope.*`
Fig. 2(a)	Frequency synchronization fraction over the $(K,\Delta\omega)$ grid	`08_map_frequency.*`
Fig. 2(b)	Fraction certified by Theorem II.1 on the same grid	`09_map_certificate.*`
Fig. 3(a)	Lifted phase diameter with unrestricted initial phases	`11_unrestricted_phase_diameter.*`
Fig. 3(b)	Frequency diameter and numerical threshold for the same trajectory	`12_unrestricted_frequency_diameter.*`
Fig. 4(a)	Baseline and modified-network phase diameters	`13_modified_network_phase.*`
Fig. 4(b)	Baseline and modified-network frequency diameters	`14_modified_network_frequency.*`
Table I	Full-box deterministic certificate	`baseline_certificate.csv`, `configuration.mat`
Table II	Coupling-gain experiments (100 paired trials per gain)	`gain_trials.csv`, `gain_uniform_certificate.csv`
Table III	ER, WS, and BA graph ensembles (30 trials per setting)	`topology_trials.csv`, `topology_summary.csv`
Table IV	Fixed-expected-degree ER networks (30 paired trials per size)	`fixed_degree_trials.csv`, `fixed_degree_summary.csv`
For figure captions, exact parameters, auxiliary plots, and the relationship between the paper and code, see FIGURES.md. The plot previews below were cropped from the manuscript PDF; they are not output files from a fresh execution of this repository.
Manuscript figure previews
Fig. 1 — Certified phase and frequency estimates
![Fig. 1: certified example](figures/manuscript_fig1.png)
Fig. 2 — Gain–frequency mismatch grid
![Fig. 2: synchronization and certification fractions](figures/manuscript_fig2.png)
Fig. 3 — Unrestricted initial phases
![Fig. 3: lifted phase and frequency diameters](figures/manuscript_fig3.png)
Fig. 4 — Modified source weights
![Fig. 4: baseline versus modified network](figures/manuscript_fig4.png)
Requirements
MATLAB with `ode45`, `tiledlayout`, `exportgraphics`, and `savefig`.
No additional MATLAB toolboxes are required by the source code.
The code uses graphics export options supported by recent MATLAB versions; compatibility with older releases and GNU Octave has not been checked.
Running the experiments
Place `run_signed_kuramoto.m` in the MATLAB working directory or add its folder to the MATLAB path.
```matlab
% Algebraic checks without ODE integration
[results, outdir] = run_signed_kuramoto('selftest');

% Reduced trial counts for a preliminary check
[results, outdir] = run_signed_kuramoto('quick', [], struct('Visible','off'));

% All experiments reported in manuscript Section IV
[results, outdir] = run_signed_kuramoto('full', [], struct('Visible','off'));
```
The `full` configuration uses the 20 trials per grid point for Fig. 2, 300 unrestricted-initial-phase trials for Fig. 3, 100 trials per gain for Table II, and 30 trials per graph setting for Tables III–IV, as specified in the manuscript. The initial sampling and network parameters are set within the code.
To run only the experiments from individual sections:
```matlab
% Section IV-A and IV-B (Figs. 1–4 and Table I)
[results, outdir] = run_signed_kuramoto('ab', [], struct('Visible','off'));

% Section IV-C (Tables II–IV and auxiliary plots)
[results, outdir] = run_signed_kuramoto('c', [], struct('Visible','off'));

% Fig. 2 and Fig. 3 alone (without baseline and Section IV-C experiments)
[results, outdir] = run_signed_kuramoto('additional', [], struct('Visible','off'));
```
Other options: `gain`, `sparse`, `topology`, `degree`, `modified`, and `core`. The `replot` mode regenerates figures from saved numerical results, where applicable:
```matlab
[results, outdir] = run_signed_kuramoto('replot', 'PATH_TO_PREVIOUS_RESULTS', struct('Visible','off'));
```
Override numerical settings with the third argument, for example:
```matlab
opts = struct('Visible','off','Seed',20261010,'TMax',150);
[results, outdir] = run_signed_kuramoto('ab', 'results_AB', opts);
```
New runs require an output folder without an existing `results.mat`. A `full` run can take substantial time and disk space because individual trial inputs are saved.
Model and numerical criteria
The original model is
$$
\dot\theta_i=\omega_i+\frac{K}{N}s_i\sum_{j=1}^{N}a_{ij}\sin(\theta_j-\theta_i),
\qquad s_i\in{-1,+1}.
$$
The interaction weight $a_{ij}\geq 0$ represents input from oscillator $j$ to oscillator $i$. The code integrates the transformed dynamics in the rotating frame
$$
\widetilde\theta_i(t)=\theta_i(t)-\frac{1-s_i}{2}\pi-\bar\omega t,
\qquad \bar\omega=\frac{1}{N}\sum_i\omega_i,
$$
with signed source matrix $B=A\operatorname{diag}(s)$. The transformed dynamics remain non-cooperative: conformist sources contribute attraction and contrarian sources contribute repulsion. The phase diameter uses unwrapped real-valued lifts of the transformed phases; circular-width diagnostics are stored separately.
The solver uses `ode45` with relative tolerance $10^{-9}$ and absolute tolerance $10^{-11}$. Each trial runs for at least 20 time units and stops at the first completed integration block satisfying
$$
\max_{t\in[t_f-5,t_f]}D(t)\leq 10^{-6},
$$
or when the simulation reaches $t_f=100$. Certification by Theorem II.1 is evaluated independently of this finite-horizon numerical criterion. Trials not meeting the numerical tolerance by the stopping time are not, on that basis alone, proved to be non-synchronizing.
Output files and reproducibility
The default output folder is `results_signed_kuramoto_<mode>_<timestamp>/`. Output files depend on the selected mode:
`configuration.mat`: run parameters, MATLAB version, timestamp, and seed.
`results.mat`: experiment records, summaries, and selected trajectories.
`trials/*.mat`: adjacency weights, oscillator signs, initial phases, natural frequencies, and individual results.
`*_trials.csv`, `*_summary.csv`: trial-level records and aggregate statistics.
`.fig`, `.pdf`, `.eps`, `.png`: generated plots (MATLAB editable, vector, and 600-dpi raster formats).
The default seed is `20261010`; trial seeds are assigned deterministically and input arrays are saved for each run. The baseline example uses $N_p=7$, $N_m=3$, $K=10$, $C=2.2$, and $C_1=0.75$. For the gain study, the code selects certificate radii algorithmically, consistent with the method described in Section IV-C of the paper.
Figure output layouts may differ from the two-column composite arrangement in the manuscript: several paper panels are exported as individual plots. Manuscript tables were generated from earlier simulation runs; the provided code has not been executed or numerically compared against every printed table entry in this release. Reproduction of the manuscript's exact stochastic summaries should be checked in MATLAB before making a stronger claim.
Citation and reuse
See CITATION.cff for software citation metadata. The article's journal and DOI may be added after publication. Code use and redistribution conditions depend on the license chosen by the copyright holders; a `LICENSE.template` file is included as a reminder to select a license before granting additional rights.
