ufunction [results, outdir] = run_signed_kuramoto(mode, outdir, user)
% Frequency synchronization of signed Kuramoto networks via a pi-state transformation.
% A(i,j) is the directed weight from node j to node i.
% Phase diameters are computed in unwrapped transformed coordinates.

if nargin < 1 || isempty(mode), mode = 'full'; end
if nargin < 2 || isempty(outdir)
    if strcmpi(mode,'replot')
        error('Specify the output folder of a previous run for replot mode.');
    end
    home = getenv('USERPROFILE');
    if isempty(home), home = getenv('HOME'); end
    if isempty(home), home = pwd; end
    outdir = fullfile(home,'signed_kuramoto_results', ...
        ['results_signed_kuramoto_' lower(mode) '_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if nargin < 3 || isempty(user), user = struct(); end

if strcmpi(mode,'replot')
    dataFile = fullfile(outdir,'results.mat');
    assert(exist(dataFile,'file')==2, 'results.mat was not found in the specified folder.');
    loaded = load(dataFile,'results');
    results = loaded.results;
    cfg = results.configuration;
    permitted = {'Visible','FigureWidth','FigureHeight'};
    names = fieldnames(user);
    for k = 1:numel(names)
        assert(ismember(names{k},permitted), 'Unknown replot option: %s', names{k});
        cfg.(names{k}) = user.(names{k});
    end
    if isfield(results,'baselineTrajectory')
        qb = results.baselineTrajectory;
        cb = results.baselineCertificate;
        plot_baseline(qb,cb,cfg,outdir);
        plot_pi_snapshots(qb,results.baselineInputs.s,cfg,outdir);
        if isfield(results,'frequency20')
            plot_omega20_blue_only(qb,results.frequency20,cb,cfg,outdir);
        end
        if isfield(results,'modifiedNetworkTrajectory')
            plot_modified_network(qb,results.modifiedNetworkTrajectory,cfg,outdir);
        end
    end
    if isfield(results,'gain'), plot_gain(results.gain,cfg,outdir); end
    if isfield(results,'sparseSize'), plot_sparse(results.sparseSize,cfg,outdir); end
    if isfield(results,'topologySummary'), plot_topology(results.topologySummary,cfg,outdir); end
    if isfield(results,'degreeSummary'), plot_degree(results.degreeSummary,cfg,outdir); end
    if isfield(results,'map'), plot_map(results.map,cfg,outdir); end
    if isfield(results,'initialIllustration')
        initialData = load(fullfile(outdir,'unrestricted_initial_illustration.mat'), ...
            'illustrative','A','s');
        [pInitial,mInitial] = pair_margins(initialData.A,initialData.s);
        plot_unrestricted_initial(initialData.illustrative.trajectory,pInitial,mInitial,cfg,outdir);
    end
    fprintf('Figures regenerated in %s\n',outdir);
    return
end

cfg = defaults(mode);
names = fieldnames(user);
for k = 1:numel(names)
    assert(isfield(cfg,names{k}), 'Unknown option: %s', names{k});
    cfg.(names{k}) = user.(names{k});
end
assert(cfg.TMax >= cfg.MinTime && cfg.MinTime > cfg.TailWindow);
assert(cfg.OutputStep > 0 && cfg.MaxStep > 0);
assert(cfg.RelTol > 0 && cfg.AbsTol > 0 && cfg.FrequencyTol > 0);
assert(~exist(fullfile(outdir,'results.mat'),'file'), ...
    'The folder already contains results.mat. Specify a new output folder.');
prepare_output_directory(outdir);
cfg.MATLABVersion = version;
cfg.Mode = mode;
cfg.Created = datestr(now,30);
save(fullfile(outdir,'configuration.mat'),'cfg');
tests = selftest();
results = struct('configuration',cfg,'selftest',tests);
if strcmpi(mode,'selftest')
    save(fullfile(outdir,'results.mat'),'results');
    fprintf('Algebraic checks passed.\n');
    return
end

runAB = any(strcmpi(mode,{'full','quick','ab','modified','core'}));
runGain = any(strcmpi(mode,{'full','quick','c','gain','core'}));
runSparse = any(strcmpi(mode,{'full','quick','c','sparse','core'}));
runTopology = any(strcmpi(mode,{'full','quick','c','topology'}));
runDegree = any(strcmpi(mode,{'full','quick','c','degree'}));
runAdditional = any(strcmpi(mode,{'full','quick','ab','additional'}));
runConvergence = any(strcmpi(mode,{'full','quick','ab','additional','core'}));

N = 10;
s = [ones(7,1); -ones(3,1)];
A = baseline_matrix(N,s);
[P,M] = pair_margins(A,s);
rng(cfg.Seed,'twister');
omega = rand(N,1)-0.5;
x0 = (2*rand(N,1)-1)*pi/3;
baselineSeed = cfg.Seed;

if runAB
    fprintf('Sections A-B: baseline and perturbed-network simulations\n');
initialTheta = x0+(1-s)*pi/2;
save(fullfile(outdir,'baseline_sample.mat'), ...
    'baselineSeed','omega','x0','initialTheta','s','A');

[P,M] = pair_margins(A,s);
c = certificate(P,M,N,10,2*pi/3,1,2.2,0.75);
assert(c.FrequencyCertified);
q = simulate(A,s,omega,x0,10,c,cfg);
results.baseline = record(q,c,A,s,omega,x0,10,cfg, ...
    'baseline','deterministic',1,1,baselineSeed,NaN);
results.baselineInputs = struct('Seed',baselineSeed,'omega',omega, ...
    'thetaTilde0',x0,'theta0',initialTheta,'s',s);
results.baselineTrajectory = q;
save_trial(outdir,results.baseline,A,s,omega,x0,c,q);
plot_baseline(q,c,cfg,outdir);
plot_pi_snapshots(q,s,cfg,outdir);

% Only the natural frequencies change in this comparison.
omega20 = 20*omega;
c20 = certificate(P,M,N,10,span(x0),span(omega20),2.2,0.75);
q20 = simulate(A,s,omega20,x0,10,c20,cfg);
results.frequency20 = q20;
save(fullfile(outdir,'frequency20_baseline_sample.mat'), ...
    'omega','omega20','x0','A','s','q20','c20','baselineSeed');
fprintf('20x-frequency experiment: Delta_omega=%.6f, tail frequency pass=%d\n', ...
    q20.DeltaOmega,q20.frequencyTailPass);
% The reference bounds in this plot apply to the baseline system only.
plot_omega20_blue_only(q,q20,c,cfg,outdir);

Amodified = modified_source_weight_matrix(N,s);
[Pmodified,Mmodified] = pair_margins(Amodified,s);
cModified = certificate(Pmodified,Mmodified,N,10, ...
    2*pi/3,1,2.2,0.75);
assert(~cModified.PhaseCertified && ~cModified.FrequencyCertified, ...
    'The modified weights unexpectedly satisfy the baseline certificate.');
qModified = simulate(Amodified,s,omega,x0,10,cModified,cfg);

Bmodified = Amodified*diag(s);
Jmodified = Bmodified-diag(sum(Bmodified,2));
modifiedEigenvalues = sort(real(eig(Jmodified)),'descend');
positiveEigenvalues = modifiedEigenvalues(modifiedEigenvalues>1e-9);
assert(numel(positiveEigenvalues)>=2, ...
    'Expected two positive eigenvalues of the modified Jacobian.');
modifiedPairMargin = min(Pmodified-Mmodified);
modifiedAveragedMargin = averaged_network_margin(Amodified,s);
assert(abs(modifiedPairMargin+1.1771475243483405)<1e-8, ...
    'Unexpected pairwise margin for the prescribed modified weights.');
assert(abs(modifiedAveragedMargin+4.604323563895003)<1e-8, ...
    'Unexpected averaged margin for the prescribed modified weights.');

results.modifiedNetwork = record(qModified,cModified,Amodified,s,omega,x0, ...
    10,cfg,'modified_network','modified_source_weights',1,1,baselineSeed,NaN);
results.modifiedNetworkTrajectory = qModified;
results.modifiedNetworkEigenvalues = modifiedEigenvalues;
save_trial(outdir,results.modifiedNetwork,Amodified,s,omega,x0,cModified,qModified);
save(fullfile(outdir,'modified_network_comparison.mat'), ...
    'A','Amodified','s','omega','x0','q','qModified','c','cModified', ...
    'baselineSeed','modifiedEigenvalues','modifiedPairMargin', ...
    'modifiedAveragedMargin');
modifiedSummary = table(modifiedPairMargin,modifiedAveragedMargin,cModified.Gamma, ...
    positiveEigenvalues(1),positiveEigenvalues(2), ...
    q.frequencyTailPass,qModified.frequencyTailPass, ...
    q.Dfrequency(end),qModified.Dfrequency(end),qModified.t(end), ...
    'VariableNames',{'ModifiedPairMargin','ModifiedAveragedMargin', ...
    'ModifiedGammaAtC22', ...
    'PositiveEigenvalue1_KoverN','PositiveEigenvalue2_KoverN', ...
    'BaselineFrequencyTailPass','ModifiedFrequencyTailPass', ...
    'BaselineFinalDfrequency','ModifiedFinalDfrequency','ModifiedTend'});
writetable(modifiedSummary,fullfile(outdir,'modified_network_summary.csv'));
fprintf(['Modified weights: min(P-M)=%.6f, averaged margin=%.6f, ', ...
    'Gamma(2.2)=%.6f, ', ...
    'positive eigenvalues (units K/N)=(%.6f, %.6f)\n'], ...
    modifiedPairMargin,modifiedAveragedMargin,cModified.Gamma, ...
    positiveEigenvalues(1),positiveEigenvalues(2));
fprintf(['Paired same-omega simulation: baseline D(T)=%.6g, ', ...
    'modified D(T)=%.6g at modified T=%.2f; ', ...
    'modified frequency-tail pass=%d\n'], ...
    q.Dfrequency(end),qModified.Dfrequency(end), ...
    qModified.t(end),qModified.frequencyTailPass);
plot_modified_network(q,qModified,cfg,outdir);
results.baselineCertificate = c;
results.frequency20Certificate = c20;
end

if runGain, run_gain_experiment(); end
if runSparse, run_sparse_experiment(); end
if runTopology, run_topology_experiment(); end
if runDegree, run_degree_experiment(); end

if runAdditional
    fprintf('Additional: gain-mismatch map and unrestricted initial phases\n');
    rows = struct([]);
    for rep = 1:cfg.MapTrials
        seed = cfg.Seed+400000+rep; rng(seed,'twister');
        baseOm = rand(N,1); baseOm = (baseOm-mean(baseOm))/span(baseOm);
        xx = (2*rand(N,1)-1)*pi/3;
        for jd = 1:numel(cfg.Mismatches)
            om = cfg.Mismatches(jd)*baseOm;
            for ik = 1:numel(cfg.MapGains)
                K = cfg.MapGains(ik); id = (jd-1)*numel(cfg.MapGains)+ik;
                c = certificate(P,M,N,K,span(xx),span(om),[],[]);
                q = simulate(A,s,om,xx,K,c,cfg);
                r = record(q,c,A,s,om,xx,K,cfg,'map','deterministic', ...
                    id,rep,seed,NaN);
                rows = append(rows,r); save_trial(outdir,r,A,s,om,xx,c,[]);
            end
        end
    end
    results.map = rows; write_rows(outdir,'map_trials',rows);
    plot_map(rows,cfg,outdir);
    rows = struct([]);
    illustrative = struct([]); bestIllustrationScore = -1;
    growthTolerance = 0.05;
    for rep = 1:cfg.InitialTrials
        seed = cfg.Seed+500000+rep; rng(seed,'twister');
        theta0 = 2*pi*rand(N,1); xx = theta0-(1-s)*pi/2;
        c = certificate(P,M,N,10,span(xx),span(omega),[],[]);
        q = simulate(A,s,omega,xx,10,c,cfg);
        r = record(q,c,A,s,omega,xx,10,cfg,'initial','deterministic',1,rep,seed,NaN);
        r.LiftedDiameterIncrease = max(q.Dtheta)-q.Dtheta(1);
        r.LiftedDiameterGrowth = r.LiftedDiameterIncrease>growthTolerance;
        rows = append(rows,r); save_trial(outdir,r,A,s,omega,xx,c,[]);
        if r.LiftedDiameterGrowth && r.PiPatternTailPass
            score=3;
        elseif r.LiftedDiameterGrowth && r.FrequencyTailPass
            score=2;
        elseif r.LiftedDiameterGrowth
            score=1;
        else
            score=0;
        end
        if score>bestIllustrationScore
            illustrative = struct('Trial',rep,'Seed',seed, ...
                'theta0',theta0,'thetaTilde0',xx,'trajectory',q, ...
                'Record',r,'Score',score);
            bestIllustrationScore=score;
        end
    end
    results.initial = rows; write_rows(outdir,'initial_trials',rows);
    nTrials=numel(rows); nCertified=sum([rows.Certified]);
    nFrequency=sum([rows.FrequencyTailPass]);
    nPiPattern=sum([rows.PiPatternTailPass]);
    nGrowth=sum([rows.LiftedDiameterGrowth]);
    nUncertFrequency=sum(~[rows.Certified] & [rows.FrequencyTailPass]);
    nUncertPiPattern=sum(~[rows.Certified] & [rows.PiPatternTailPass]);
    initialSummary = table(nTrials,nCertified,nFrequency,nPiPattern, ...
        nGrowth,nUncertFrequency,nUncertPiPattern, ...
        nFrequency/nTrials,nPiPattern/nTrials, ...
        'VariableNames',{'Trials','Certified','FrequencyTailPass', ...
        'PiPatternTailPass','LiftedDiameterGrowth', ...
        'UncertifiedFrequencyTailPass','UncertifiedPiPatternTailPass', ...
        'FrequencyTailFraction','PiPatternTailFraction'});
    results.initialSummary = initialSummary;
    writetable(initialSummary,fullfile(outdir,'initial_trials_summary.csv'));
    fprintf('Unrestricted initial phases: frequency criterion %d/%d, ', ...
        nFrequency,nTrials);
    fprintf('pi-pattern criterion %d/%d, lifted diameter growth %d/%d\n', ...
        nPiPattern,nTrials,nGrowth,nTrials);
    if nFrequency<nTrials || nPiPattern<nTrials
        fprintf('Some trials did not meet the finite-horizon criteria.\n');
    end
    save(fullfile(outdir,'unrestricted_initial_illustration.mat'), ...
        'illustrative','omega','A','s','growthTolerance');
    writetable(struct2table(illustrative.Record), ...
        fullfile(outdir,'unrestricted_initial_illustration.csv'));
    results.initialIllustration = illustrative.Record;
    plot_unrestricted_initial(illustrative.trajectory,P,M,cfg,outdir);
    fprintf('Unrestricted illustration: trial %d (seed %d), ', ...
        illustrative.Trial,illustrative.Seed);
    fprintf('lifted diameter growth %.4f\n', ...
        illustrative.Record.LiftedDiameterIncrease);
    if illustrative.Score==0
        fprintf('No trial had lifted-diameter growth above %.3f.\n',growthTolerance);
    end
end

if runConvergence
fprintf('Numerical convergence checks\n');
results.convergence = convergence_check(cfg);
write_rows(outdir,'numerical_convergence',results.convergence);
end
if isfield(results,'baseline')
    write_rows(outdir,'baseline_certificate',results.baseline);
end
save(fullfile(outdir,'results.mat'),'results');
fprintf('Completed. Results in %s\n',outdir);

    function run_gain_experiment()
    fprintf('Section C: paired gain sweep and three different rate estimates\n');
    rows = struct([]);
    for rep = 1:cfg.GainTrials
        seed = cfg.Seed + 1000 + rep; rng(seed,'twister');
        om = rand(N,1)-0.5; xx = (2*rand(N,1)-1)*pi/3;
        for ik = 1:numel(cfg.Gains)
            K = cfg.Gains(ik);
            c = certificate(P,M,N,K,span(xx),span(om),[],[]);
            q = simulate(A,s,om,xx,K,c,cfg);
            r = record(q,c,A,s,om,xx,K,cfg,'gain','deterministic',ik,rep,seed,NaN);
            rows = append(rows,r);
            save_trial(outdir,r,A,s,om,xx,c,[]);
        end
    end
    results.gain = rows;
    write_rows(outdir,'gain_trials',rows);
    plot_gain(rows,cfg,outdir);
    save(fullfile(outdir,'results.mat'),'results');
    end

    function run_sparse_experiment()
    fprintf('Section C: sparse minority-conformist network and design rules\n');
    rows = struct([]); rTarget = 2.5; C = 2.2; C1 = 0.75; ep = 0.2;
    for iz = 1:numel(cfg.SparseSizes)
        n = cfg.SparseSizes(iz); np = max(1,round(0.1*n));
        ss = -ones(n,1); ss(1:np) = 1;
        seed = cfg.Seed + 2000 + iz; rng(seed,'twister');
        om = linspace(-0.5,0.5,n)'; om = om(randperm(n));
        xx = linspace(-pi/3,pi/3,n)'; xx = xx(randperm(n));
        wRequired = max((ep+1/(2*sin(C1/2)))/cos(C/2), ...
            (ep+rTarget)/cos(C1));
        w = 1.05*wRequired; AA = broadcast_matrix(n,1,n,w,ep);
        [p,m] = pair_margins(AA,ss);
        assert(max(abs(p-w))<1e-11 && max(abs(m-ep))<1e-11);
        for ig = 1:2
            if ig == 1, K = n; else, K = 10; end
            c = certificate(p,m,n,K,span(xx),span(om),C,C1);
            q = simulate(AA,ss,om,xx,K,c,cfg);
            r = record(q,c,AA,ss,om,xx,K,cfg,'sparse_size', ...
                'two_sources',2*iz-2+ig,ig,seed,NaN);
            r.TargetRate = rTarget; r.SourceWeight = w;
            rows = append(rows,r); save_trial(outdir,r,AA,ss,om,xx,c,[]);
        end
    end
    results.sparseSize = rows; write_rows(outdir,'sparse_size',rows);
    n = 100; ss = -ones(n,1); ss(1:10) = 1;
    om = linspace(-0.5,0.5,n)'; xx = linspace(-pi/3,pi/3,n)';
    rows = struct([]); factors = [0.5 0.9 1.05];
    for kind = 1:2
        for j = 1:numel(factors)
            if kind == 1
                w = 0.9;
                KRequired = n*max(1/(2*(w*cos(C/2)-ep)*sin(C1/2)), ...
                    rTarget/(w*cos(C1)-ep));
                K = factors(j)*KRequired; threshold = KRequired;
            else
                K = n;
                wRequired = max((ep+n/(2*K*sin(C1/2)))/cos(C/2), ...
                    (ep+n*rTarget/K)/cos(C1));
                w = factors(j)*wRequired; threshold = wRequired;
            end
            AA = broadcast_matrix(n,1,n,w,ep); [p,m] = pair_margins(AA,ss);
            c = certificate(p,m,n,K,span(xx),span(om),C,C1);
            q = simulate(AA,ss,om,xx,K,c,cfg);
            r = record(q,c,AA,ss,om,xx,K,cfg,'design','two_sources', ...
                3*(kind-1)+j,kind,cfg.Seed,NaN);
            r.TargetRate = rTarget; r.SourceWeight = w;
            r.DesignFactor = factors(j); r.DesignThreshold = threshold;
            rows = append(rows,r); save_trial(outdir,r,AA,ss,om,xx,c,[]);
        end
    end
    results.design = rows; write_rows(outdir,'design_trials',rows);
    plot_sparse(results.sparseSize,cfg,outdir);
    end

    function run_topology_experiment()
        fprintf('Section C: random graph ensembles, with retained failures\n');
        rows = struct([]); id = 0; families = {'ER','WS','BA'};
        for jf = 1:3
            if jf == 3, densities = cfg.BADensities; else, densities = cfg.Densities; end
            for iz = 1:numel(cfg.Sizes)
                n = cfg.Sizes(iz);
                for jd = 1:numel(densities)
                    density = densities(jd); id = id+1;
                    for rep = 1:cfg.GraphTrials
                        seed = cfg.Seed+100000+id*1000+rep; rng(seed,'twister');
                        H = make_graph(families{jf},n,density,cfg.WSRewire);
                        [AA,ss,om,xx] = random_inputs(H);
                        degree = nnz(H)/n;
                        if degree > 0, K = 10*n/degree; else, K = 0; end
                        [p,m] = pair_margins(AA,ss);
                        c = certificate(p,m,n,K,span(xx),span(om),[],[]);
                        q = simulate(AA,ss,om,xx,K,c,cfg);
                        r = record(q,c,AA,ss,om,xx,K,cfg, ...
                            'topology',families{jf},id,rep,seed,density);
                        rows = append(rows,r); save_trial(outdir,r,AA,ss,om,xx,c,[]);
                    end
                    fprintf('  %s, N=%d, target density %.2f completed\n',families{jf},n,density);
                end
            end
        end
        results.topology = rows;
        write_rows(outdir,'topology_trials',rows);
        results.topologySummary = ensemble_summary(rows);
        write_rows(outdir,'topology_summary',results.topologySummary);
        plot_topology(results.topologySummary,cfg,outdir);
        save(fullfile(outdir,'results.mat'),'results');
    end

    function run_degree_experiment()
        fprintf('Section C: fixed expected degree, paired gain normalization\n');
        rows = struct([]);
        for iz = 1:numel(cfg.DegreeSizes)
            n = cfg.DegreeSizes(iz); prob = min(1,10/(n-1));
            for rep = 1:cfg.GraphTrials
                seed = cfg.Seed+300000+iz*1000+rep; rng(seed,'twister');
                H = make_graph('ER',n,prob,cfg.WSRewire);
                [AA,ss,om,xx] = random_inputs(H); degree = nnz(H)/n;
                [p,m] = pair_margins(AA,ss);
                if degree > 0, normalizedK = 10*n/degree; else, normalizedK = 0; end
                for ig = 1:2
                    if ig == 1, K = 10; else, K = normalizedK; end
                    c = certificate(p,m,n,K,span(xx),span(om),[],[]);
                    q = simulate(AA,ss,om,xx,K,c,cfg);
                    r = record(q,c,AA,ss,om,xx,K,cfg,'fixed_degree','ER', ...
                        2*iz-2+ig,rep,seed,prob);
                    r.GainNormalized = (ig==2);
                    rows = append(rows,r); save_trial(outdir,r,AA,ss,om,xx,c,[]);
                end
            end
        end
        results.fixedDegree = rows;
        write_rows(outdir,'fixed_degree_trials',rows);
        results.degreeSummary = ensemble_summary(rows);
        write_rows(outdir,'fixed_degree_summary',results.degreeSummary);
        plot_degree(results.degreeSummary,cfg,outdir);
    end
end

function cfg = defaults(mode)
assert(any(strcmpi(mode,{'quick','full','ab','c','gain','sparse','topology','degree', ...
    'additional','core','modified','selftest'})), 'Unknown mode.');
cfg = struct('Seed',20261010,'TMax',100,'MinTime',20,'TailWindow',5, ...
    'BlockTime',10,'OutputStep',0.05,'MaxStep',0.1,'RelTol',1e-9, ...
    'AbsTol',1e-11,'FrequencyTol',1e-6,'EnvelopeTol',2e-7, ...
    'PiRadius',0.75,'Visible','on','WSRewire',0.1, ...
    'FigureWidth',3.5,'FigureHeight',3.0,'PanelHeight',2.4,'StackedHeight',5.8, ...
    'GainTrials',6, ...
    'GraphTrials',5,'MapTrials',3,'InitialTrials',10, ...
    'Gains',[0.5 1 2 5 10 20],'Sizes',[20 50], ...
    'Densities',[0.3 0.7 0.9],'BADensities',[0.1 0.3 0.45], ...
    'DegreeSizes',[20 100],'SparseSizes',[10 50 100], ...
    'MapGains',[1 3 10],'Mismatches',[0.5 1 4]);
if any(strcmpi(mode,{'full','c','gain','sparse','topology','degree','additional','ab'}))
    cfg.GainTrials = 100; cfg.GraphTrials = 30;
    cfg.MapTrials = 20; cfg.InitialTrials = 300;
    cfg.Gains = [0.5 1 2 3 5 10 20 40];
    cfg.Sizes = [20 50 100 200]; cfg.Densities = [0.1 0.3 0.5 0.7 0.9];
    cfg.DegreeSizes = [20 50 100 200 400]; cfg.SparseSizes = [10 50 100 200];
    cfg.MapGains = [0.5 1 2 3 5 10 20]; cfg.Mismatches = [0.25 0.5 1 2 4];
elseif strcmpi(mode,'core')
    cfg.GainTrials = 20;
end
end

function A = baseline_matrix(N,s)
z = (1:N)'*(1:N);
A = repmat(0.2+0.7*(s(:)'>0),N,1)+0.2*sin(z);
A(1:N+1:end) = 0;
end

function A = modified_source_weight_matrix(N,s)
indices = (1:N)'*(1:N);
sourceWeights = 0.6-0.3*(s(:)'>0);
A = repmat(sourceWeights,N,1)+0.2*sin(indices);
A(1:N+1:end) = 0;
end

function value = averaged_network_margin(A,s)
rowPositive = sum(A(:,s>0),2);
rowNegative = sum(A(:,s<0),2);
maxRowDistance = 0;
for ii=1:size(A,1)
    for jj=ii+1:size(A,1)
        maxRowDistance = max(maxRowDistance, ...
            sum(abs(A(ii,:)-A(jj,:))));
    end
end
value = min(rowPositive)-max(rowNegative)-maxRowDistance;
end

function A = broadcast_matrix(N,h,q,w,ep)
A = sparse(N,N); A(:,h) = w; A(:,q) = ep;
A(h,h) = 0; A(q,q) = 0;
assert(nnz(A)==2*(N-1));
end

function [P,M] = pair_margins(A,s)
% Pairwise positive and negative contributions include both direct links.
N = numel(s); A = full(A); A(1:N+1:end) = 0;
P = zeros(N); M = zeros(N);
for k = 1:N
    v = A(:,k);
    if s(k)>0
        v(k) = max(v); P = P+bsxfun(@min,v,v');
    else
        M = M+bsxfun(@max,v,v');
    end
end
idx = triu(true(N),1); P = P(idx); M = M(idx);
end

function c = certificate(P,M,N,K,d0,delta,C,C1)
c = struct('C',NaN,'C1',NaN,'Gamma',NaN,'Lambda',NaN, ...
    'Margin0',min(P-M),'Cinf',NaN,'CinfImproved',NaN, ...
    'Rate',NaN,'LimitingRate',NaN,'EntryBound',NaN,'Capacity',NaN, ...
    'PhaseCertified',false,'FrequencyCertified',false,'GainCertifiable',false, ...
    'Delta',delta,'InitialDiameter',d0,'ZeroPairs',mean(P<=1e-14));
gamma = @(z) min(P*cos(z/2)-M);
lambda = @(z) min(P*cos(z)-M);
c.GainCertifiable = d0 < pi && gamma(d0)>0;
if isempty(C)
    if ~c.GainCertifiable, return; end
    ratio = max(M./P); critical = acos(ratio);
    lo = d0+1e-8; hi = min(pi,2*critical)-1e-8;
    if lo>=hi, return; end
    capacityShape = @(z) 2*gamma(z)*sin(min(z,critical)/2);
    opt = fminbnd(@(z)-capacityShape(z),lo,hi,optimset('TolX',1e-10));
    choices = [lo opt hi];
    if critical>lo && critical<hi, choices(end+1) = critical; end
    values = arrayfun(capacityShape,choices);
    [~,j] = max(values); C = choices(j);
end
c.C = C; c.Gamma = gamma(C);
c.Capacity = (2*K/N)*c.Gamma*sin(C/2);
c.PhaseCertified = C>0 && C<pi && d0<C && c.Gamma>0 && delta<c.Capacity;
if ~c.PhaseCertified, return; end
c.Cinf = 2*asin(N*delta/(2*K*c.Gamma));
c.LimitingRate = (K/N)*lambda(c.Cinf);
if isempty(C1)
    critical = acos(max(M./P));
    if c.Cinf>=min(C,critical), return; end
    C1 = min(C,(c.Cinf+critical)/2);
end
c.C1 = C1; c.Lambda = lambda(C1);
c.FrequencyCertified = C1>c.Cinf && C1<=C && C1<pi/2 && c.Lambda>0;
if c.FrequencyCertified
    c.Rate = K*c.Lambda/N;
    denominator = (2*K/N)*c.Gamma*sin(C1/2)-delta;
    c.EntryBound = max(d0-C1,0)/denominator;
    c.CinfImproved = 2*asin(N*delta/(2*K*gamma(C1)));
end
end

function f = rhs_transformed(x,B,omegaTilde,kappa)
f = omegaTilde+kappa*(cos(x).*(B*sin(x))-sin(x).*(B*cos(x)));
end

function f = rhs_original(theta,A,s,omega,K)
f = omega+(K/numel(s))*s.*(cos(theta).*(A*sin(theta))-sin(theta).*(A*cos(theta)));
end

function q = simulate(A,s,omega,x0,K,c,cfg)
% Integrate the transformed dynamics and evaluate finite-horizon diagnostics.
N = numel(s); B = sparse(A)*spdiags(s,0,N,N); reference = mean(omega);
fun = @(t,x) rhs_transformed(x,B,omega-reference,K/N);
opts = odeset('RelTol',cfg.RelTol,'AbsTol',cfg.AbsTol,'MaxStep',cfg.MaxStep);
entry = NaN; xEntry = [];
if c.FrequencyCertified
    if span(x0)<=c.C1, entry = 0; xEntry = x0; end
    opts = odeset(opts,'Events',@(t,x) entry_event(t,x,c.C1));
end
t = 0; X = x0'; start = 0; current = x0;
while start < cfg.TMax-1e-10
    finish = min(cfg.TMax,start+cfg.BlockTime);
    grid = linspace(start,finish,max(2,ceil((finish-start)/cfg.OutputStep)+1));
    if c.FrequencyCertified
        [tb,xb,te,xe] = ode45(fun,grid,current,opts);
        if isnan(entry) && ~isempty(te), entry = te(1); xEntry = xe(1,:)'; end
    else
        [tb,xb] = ode45(fun,grid,current,opts);
    end
    assert(tb(end)>=finish-1e-8 && all(isfinite(xb(:))), ...
        'ODE integration did not finish; inspect this trial.');
    t = [t; tb(2:end)]; X = [X; xb(2:end,:)];
    start = tb(end); current = xb(end,:)';
    tail = t>=start-cfg.TailWindow;
    Vtail = velocities(X(tail,:),B,omega-reference,K/N);
    dvTail = max(Vtail,[],2)-min(Vtail,[],2);
    if start>=cfg.MinTime && max(dvTail)<=cfg.FrequencyTol, break; end
end
V = velocities(X,B,omega-reference,K/N);
dtheta = max(X,[],2)-min(X,[],2); dv = max(V,[],2)-min(V,[],2);
tail = t>=t(end)-cfg.TailWindow;
freqOK = max(dv(tail))<=cfg.FrequencyTol;
arcTail = zeros(sum(tail),1); ids = find(tail);
for j = 1:numel(ids), arcTail(j) = circular_width(X(ids(j),:)); end
rp = abs(mean(exp(1i*X(:,s>0)),2)); rm = abs(mean(exp(1i*X(:,s<0)),2));
gap = abs(angle(mean(exp(1i*X(:,s>0)),2).*conj(mean(exp(1i*X(:,s<0)),2))));
gap(rp<1e-8 | rm<1e-8) = NaN;
envelope = NaN(size(t)); violation = NaN;
if c.FrequencyCertified && ~isnan(entry)
    vEntry = fun(entry,xEntry); v0 = span(vEntry); after = t>=entry;
    envelope(after) = v0*exp(-c.Rate*(t(after)-entry));
    violation = max(dv(after)-envelope(after));
end
fit = struct('rate',NaN,'r2',NaN,'first',NaN,'last',NaN);
gapSpectral = NaN;
if freqOK
    anchor = 0; if ~isnan(entry), anchor = entry; end
    fit = fit_rate(t,dv,anchor,cfg);
    gapSpectral = spectral_gap(A,s,X(end,:)',K);
end
q = struct('t',t,'thetaTilde',X,'Dtheta',dtheta,'Dfrequency',dv, ...
    'Rp',rp,'Rm',rm,'PiGapError',gap,'entry',entry, ...
    'envelope',envelope,'envelopeViolation',violation,'frequencyTailPass',freqOK, ...
    'piPatternTailPass',freqOK && max(arcTail)<=cfg.PiRadius, ...
    'arcTailMax',max(arcTail),'arcFinal',arcTail(end),'fit',fit, ...
    'spectralGap',gapSpectral,'omegaReference',reference, ...
    'DeltaOmega',span(omega));
if c.FrequencyCertified
    assert(max(dtheta)<c.C+cfg.EnvelopeTol,'Certified phase bound violated.');
    if cfg.TMax>c.EntryBound+cfg.OutputStep
        assert(isfinite(entry),'Certified entry not observed; inspect integration.');
    end
    if isfinite(entry)
        assert(entry<=c.EntryBound+1e-5*(1+c.EntryBound),'Entry-time bound violated.');
        assert(violation<=cfg.EnvelopeTol,'Frequency envelope violated.');
    end
end
end

function V = velocities(X,B,omegaTilde,kappa)
S = sin(X); C = cos(X);
V = bsxfun(@plus,kappa*(C.*(B*S')'-S.*(B*C')'),omegaTilde');
end

function [value,isterminal,direction] = entry_event(~,x,C1)
value = span(x)-C1; isterminal = 0; direction = -1;
end

function z = span(x)
z = max(x)-min(x);
end

function w = circular_width(x)
z = sort(mod(x(:),2*pi));
w = 2*pi-max(diff([z;z(1)+2*pi]));
end

function g = spectral_gap(A,s,x,K)
% Remove the single rotational mode using phase differences.
N = numel(s); J = full(A).*repmat(s',N,1).*cos(bsxfun(@minus,x',x));
J(1:N+1:end) = 0; J = (K/N)*(J-diag(sum(J,2)));
Jq = J(2:end,2:end)-repmat(J(1,2:end),N-1,1);
g = -max(real(eig(Jq)));
end

function f = fit_rate(t,dv,anchor,cfg)
% Fit the last qualifying decay interval; no rate maximization is used.
f = struct('rate',NaN,'r2',NaN,'first',NaN,'last',NaN);
ia = find(t>=anchor,1); scale = dv(ia);
eligible = t>=anchor & dv<=0.1*scale & dv>=max(1e-8,100*cfg.AbsTol);
edges = diff([false; eligible; false]); starts = find(edges==1); ends = find(edges==-1)-1;
good = find(ends-starts+1>=8,1,'last'); if isempty(good), return; end
ix = starts(good):ends(good); y = log(dv(ix)); p = polyfit(t(ix),y,1);
res = y-polyval(p,t(ix)); ss = sum((y-mean(y)).^2);
if ss<=eps, return; end
f.r2 = 1-sum(res.^2)/ss; f.first = t(ix(1)); f.last = t(ix(end));
if f.r2>=0.98 && p(1)<0, f.rate = -p(1); end
end

function r = record(q,c,A,s,omega,x0,K,cfg,experiment,topology,id,trial,seed,density)
N = numel(s); net = A*s;
r = struct('Experiment',experiment,'Topology',topology,'CaseID',id, ...
    'Trial',trial,'Seed',seed,'N',N,'Np',sum(s>0), ...
    'TargetDensity',density,'Density',nnz(A)/(N*(N-1)), ...
    'MeanDegree',nnz(A)/N,'Edges',nnz(A),'K',K,'DeltaOmega',span(omega), ...
    'Dtheta0',span(x0),'C',c.C,'C1',c.C1,'Gamma',c.Gamma, ...
    'Lambda',c.Lambda,'Margin0',c.Margin0,'PositiveMargin',c.Margin0>0, ...
    'GainCertifiable',c.GainCertifiable,'Certified',c.FrequencyCertified, ...
    'PhaseCertified',c.PhaseCertified,'Cinf',c.Cinf, ...
    'CinfImproved',c.CinfImproved,'RateGuaranteed',c.Rate, ...
    'RateLimitSupremum',c.LimitingRate,'Capacity',c.Capacity, ...
    'EntryBound',c.EntryBound,'EntryObserved',q.entry, ...
    'RateFitted',q.fit.rate,'FitR2',q.fit.r2,'FitStart',q.fit.first, ...
    'FitEnd',q.fit.last,'SpectralGapFinal',q.spectralGap, ...
    'Tend',q.t(end),'FrequencyTailPass',q.frequencyTailPass, ...
    'PiPatternTailPass',q.piPatternTailPass,'MaxDtheta',max(q.Dtheta), ...
    'FinalDtheta',q.Dtheta(end),'FinalDfrequency',q.Dfrequency(end), ...
    'ArcFinal',q.arcFinal,'ArcTailMax',q.arcTailMax, ...
    'PiGapErrorFinal',q.PiGapError(end),'MinRp',min(q.Rp),'MinRm',min(q.Rm), ...
    'CohesiveThroughPi',max(q.Dtheta)<pi,'ZeroPairFraction',c.ZeroPairs, ...
    'HasNegativeNetInput',any(net<0),'EnvelopeViolation',q.envelopeViolation, ...
    'TailWindow',cfg.TailWindow,'GainNormalized',false, ...
    'SourceWeight',NaN,'TargetRate',NaN,'DesignFactor',NaN,'DesignThreshold',NaN);
end

function rows = append(rows,r)
if isempty(rows), rows = r; else, rows(end+1) = r; end
end

function save_trial(outdir,r,A,s,omega,thetaTilde0,certificateInfo,trajectory)
filename = sprintf('%s_%s_c%03d_r%04d.mat',r.Experiment,r.Topology,r.CaseID,r.Trial);
A = sparse(A); omegaReference = mean(omega);
theta0 = thetaTilde0+(1-s)*pi/2;
save(fullfile(outdir,'trials',filename),'A','s','omega','omegaReference', ...
    'thetaTilde0','theta0','r','certificateInfo','trajectory');
end

function write_rows(outdir,name,rows)
if ~isempty(rows), writetable(struct2table(rows),fullfile(outdir,[name '.csv'])); end
end

function H = make_graph(family,N,density,rewire)
switch family
    case 'ER'
        H = triu(rand(N)<density,1); H = H|H';
    case 'WS'
        halfDegree = min(floor((N-1)/2),max(1,round(density*(N-1)/2)));
        H = false(N);
        for i = 1:N
            for d = 1:halfDegree
                j = mod(i-1+d,N)+1; H(i,j) = true; H(j,i) = true;
            end
        end
        for i = 1:N
            for d = 1:halfDegree
                j = mod(i-1+d,N)+1;
                if rand < rewire
                    candidates = find(~H(i,:)); candidates(candidates==i) = [];
                    if ~isempty(candidates)
                        new = candidates(randi(numel(candidates)));
                        H(i,j) = false; H(j,i) = false;
                        H(i,new) = true; H(new,i) = true;
                    end
                end
            end
        end
    case 'BA'
        ms = 1:N-1; achieved = (2*ms*N-ms.*(ms+1))/(N*(N-1));
        [~,m] = min(abs(achieved-density));
        H = false(N); H(1:m+1,1:m+1) = true; H(1:N+1:end) = false;
        for v = m+2:N
            weights = sum(H(1:v-1,1:v-1),2);
            for j = 1:m
                cumulative = cumsum(weights); u = rand*cumulative(end);
                k = find(cumulative>u,1); H(v,k) = true; H(k,v) = true;
                weights(k) = 0;
            end
        end
    otherwise
        error('Unknown graph family.');
end
assert(isequal(H,H') && ~any(diag(H)));
end

function [A,s,omega,x0] = random_inputs(H)
N = size(H,1); s = -ones(N,1); order = randperm(N); s(order(1:round(.7*N))) = 1;
weights = repmat(.2+.7*(s'>0),N,1)+.2*(2*rand(N)-1);
A = sparse(double(H).*weights); omega = rand(N,1)-.5;
x0 = (2*rand(N,1)-1)*pi/3;
end

function out = ensemble_summary(rows)
out = struct([]); ids = unique([rows.CaseID]);
for id = ids
    z = rows([rows.CaseID]==id); n = numel(z);
    [lo,hi] = wilson(sum([z.FrequencyTailPass]),n);
    [clo,chi] = wilson(sum([z.Certified]),n);
    r = struct('Topology',z(1).Topology,'CaseID',id,'N',z(1).N, ...
        'TargetDensity',z(1).TargetDensity,'MeanDensity',mean([z.Density]), ...
        'GainNormalized',z(1).GainNormalized,'Trials',n, ...
        'PositiveMarginFraction',mean([z.PositiveMargin]), ...
        'GainCertifiableFraction',mean([z.GainCertifiable]), ...
        'CertifiedFraction',mean([z.Certified]),'CertCILow',clo,'CertCIHigh',chi, ...
        'FrequencyFraction',mean([z.FrequencyTailPass]), ...
        'FrequencyCILow',lo,'FrequencyCIHigh',hi, ...
        'PiPatternFraction',mean([z.PiPatternTailPass]), ...
        'Unresolved',sum(~[z.FrequencyTailPass]), ...
        'ZeroPairFraction',mean([z.ZeroPairFraction]), ...
        'NegativeNetFraction',mean([z.HasNegativeNetInput]), ...
        'NegativeNetTrials',sum([z.HasNegativeNetInput]), ...
        'NegativeNetFrequencyPass',sum([z.HasNegativeNetInput] & [z.FrequencyTailPass]), ...
        'OtherTrials',sum(~[z.HasNegativeNetInput]), ...
        'OtherFrequencyPass',sum(~[z.HasNegativeNetInput] & [z.FrequencyTailPass]), ...
        'MedianSpectralGap',finite_median([z.SpectralGapFinal]));
    out = append(out,r);
end
end

function [lo,hi] = wilson(k,n)
z = 1.95996398454005; p = k/n; denominator = 1+z^2/n;
center = (p+z^2/(2*n))/denominator;
width = z*sqrt(p*(1-p)/n+z^2/(4*n^2))/denominator;
lo = max(0,center-width); hi = min(1,center+width);
end

function m = finite_median(x)
x = x(isfinite(x)); if isempty(x), m = NaN; else, m = median(x); end
end

function prepare_output_directory(outdir)
if exist(outdir,'dir')~=7
    [ok,msg] = mkdir(outdir);
    assert(ok,'Cannot create output directory %s: %s',outdir,msg);
end
trialdir = fullfile(outdir,'trials');
if exist(trialdir,'dir')~=7
    [ok,msg] = mkdir(trialdir);
    assert(ok,'Cannot create trial directory %s: %s',trialdir,msg);
end
testfile = fullfile(outdir,'.write_test');
fid = fopen(testfile,'w');
assert(fid~=-1,'Cannot write to output directory: %s',outdir);
fclose(fid);
delete(testfile);
end

function h = newfigure(cfg,layoutKind)
if nargin<2, layoutKind = 'single'; end
if strcmp(layoutKind,'two_panel')
    figHeight = cfg.PanelHeight;
elseif strcmp(layoutKind,'stacked')
    figHeight = cfg.StackedHeight;
else
    figHeight = cfg.FigureHeight;
end
h = figure('Units','inches','Position',[1 1 cfg.FigureWidth figHeight], ...
    'Color','w','Visible',cfg.Visible,'Renderer','painters');
if strcmp(layoutKind,'single')
    axes('Parent',h,'FontName','Times New Roman','FontSize',9, ...
        'LineWidth',0.65,'TickDir','out');
    hold on; box on; grid off;
end
end

function lg = manuscript_legend(items,location)
lg = legend(items,'Interpreter','latex','Location',location, ...
    'AutoUpdate','off','Box','on');
set(lg,'FontSize',8,'Color','white');
end

function export_figure(h,outdir,name)
assert(exist(outdir,'dir')==7, ...
    'Output directory is missing: %s. Run from an extracted, permanent folder.',outdir);
set(h,'Color','w'); drawnow;
set(h,'Units','inches'); pos=get(h,'Position'); width=pos(3);
base = fullfile(outdir,name);
savefig(h,[base '.fig']);
exportgraphics(h,[base '.pdf'],'ContentType','vector', ...
    'BackgroundColor','white','Units','inches','Width',width,'Padding','tight');
exportgraphics(h,[base '.eps'],'ContentType','vector', ...
    'BackgroundColor','white','Units','inches','Width',width,'Padding','tight');
exportgraphics(h,[base '.png'],'Resolution',600, ...
    'BackgroundColor','white','Units','inches','Width',width,'Padding','tight');
end

function [items, ax] = plot_phase_panel(ax,q,c,bound,drawLegend)
blue=[0 0.4470 0.6980]; orange=[0.8350 0.3690 0];
gray=[0.42 0.42 0.42]; red=[0.6350 0.0780 0.1840];
green=[0.25 0.50 0.20]; purple=[0.4940 0.1840 0.5560];
hold(ax,'on'); box(ax,'on'); grid(ax,'off');
p1=plot(ax,q.t,q.Dtheta,'-','Color',blue,'LineWidth',1.35, ...
    'DisplayName','$D_{\theta}(t)$');
p2=plot(ax,q.t,bound,'--','Color',orange,'LineWidth',1.25, ...
    'DisplayName','$z(t)$');
p3=yline(ax,c.C,':','Color',gray,'LineWidth',0.9,'DisplayName','$C$');
p4=yline(ax,c.C1,'-.','Color',red,'LineWidth',1.0,'DisplayName','$C_1$');
p5=yline(ax,c.Cinf,'--','Color',green,'LineWidth',0.95, ...
    'DisplayName','$C_{\infty}$');
items=[p1 p2 p3 p4 p5];
phaseMax=min(2,q.t(end));
if isfinite(q.entry) && q.entry <= phaseMax
    p6=xline(ax,q.entry,':','Color',purple,'LineWidth',1.15, ...
        'DisplayName',['$t_{\mathrm{entry}}=',num2str(q.entry,'%.4f'),'$']);
    items=[items p6];
end
xlim(ax,[0 phaseMax]);
ylim(ax,[0 max(2.4,c.C+0.22)]);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$D_{\theta}(t)$','Interpreter','latex');
if drawLegend
    lg=legend(ax,items,'Interpreter','latex','Location','northeast', ...
        'AutoUpdate','off','Box','on');
    lg.NumColumns=2;
    set(lg,'FontSize',8,'Color','white');
    lg.Units = 'normalized';
    oldUnits = ax.Units;
    ax.Units = 'normalized';
    axPosition = ax.Position;
    ax.Units = oldUnits;
    legendPosition = lg.Position;
    legendPosition(2) = legendPosition(2) - 0.15*axPosition(4);
    lg.Position = legendPosition;
end
end

function [items,ax] = plot_frequency_panel(ax,q,cfg,drawLegend)
blue=[0 0.4470 0.6980]; orange=[0.8350 0.3690 0];
purple=[0.4940 0.1840 0.5560];
threshold=max(1e-10,10*cfg.AbsTol);
observed=q.Dfrequency; observed(observed<threshold)=NaN;
upper=q.envelope; upper(upper<threshold)=NaN;
hold(ax,'on'); box(ax,'on'); grid(ax,'off');
set(ax,'YScale','log');
p1=semilogy(ax,q.t,observed,'-','Color',blue,'LineWidth',1.35, ...
    'DisplayName','$D(t)$');
p2=semilogy(ax,q.t,upper,'--','Color',orange,'LineWidth',1.25, ...
    'DisplayName','Theoretical bound');
items=[p1 p2];
if isfinite(q.entry) && q.entry <= min(10,q.t(end))
    p3=xline(ax,q.entry,':','Color',purple,'LineWidth',1.15, ...
        'DisplayName','$t_{\mathrm{entry}}$');
    items=[items p3];
end
xlim(ax,[0 min(3,q.t(end))]);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$D(t)$','Interpreter','latex');
if drawLegend
    lg=legend(ax,items,'Interpreter','latex','Location','northeast', ...
        'AutoUpdate','off','Box','on');
    set(lg,'FontSize',8,'Color','white');
end
end

function plot_baseline(q,c,cfg,outdir)
opts=odeset('RelTol',cfg.RelTol,'AbsTol',cfg.AbsTol,'MaxStep',cfg.MaxStep);
a=c.Delta/sin(c.Cinf/2);
[~,bound]=ode45(@(t,z)q.DeltaOmega-a*sin(z/2), ...
    q.t,q.Dtheta(1),opts);
if isfinite(q.entry)
    fprintf('Observed baseline entry time = %.12f (seed %d)\n', ...
        q.entry,cfg.Seed);
end

h=newfigure(cfg); ax=gca;
plot_phase_panel(ax,q,c,bound,true);
text(ax,0.02,0.98,'(a)','Units','normalized', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'VerticalAlignment','top','HandleVisibility','off');
export_figure(h,outdir,'01_phase_entry');
h=newfigure(cfg); ax=gca;
plot_frequency_panel(ax,q,cfg,true);
text(ax,0.02,0.98,'(b)','Units','normalized', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'VerticalAlignment','top','HandleVisibility','off');
export_figure(h,outdir,'02_frequency_envelope');

h=newfigure(cfg,'stacked');
tl=tiledlayout(h,2,1,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl,1);
set(ax1,'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.65,'TickDir','out');
plot_phase_panel(ax1,q,c,bound,true);
text(ax1,0.02,0.98,'(a)','Units','normalized', ...
    'FontSize',9,'FontWeight','bold','VerticalAlignment','top');
ax2=nexttile(tl,2);
set(ax2,'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.65,'TickDir','out');
plot_frequency_panel(ax2,q,cfg,true);
text(ax2,0.02,0.98,'(b)','Units','normalized', ...
    'FontSize',9,'FontWeight','bold','VerticalAlignment','top');
export_figure(h,outdir,'00_phase_frequency_onecolumn');
end

function plot_omega20_blue_only(qBaseline,q20,c,cfg,outdir)
opts=odeset('RelTol',cfg.RelTol,'AbsTol',cfg.AbsTol,'MaxStep',cfg.MaxStep);
a=c.Delta/sin(c.Cinf/2);
[~,baselineComparison]=ode45( ...
    @(t,z)qBaseline.DeltaOmega-a*sin(z/2), ...
    qBaseline.t,qBaseline.Dtheta(1),opts);
h=newfigure(cfg); ax=gca;
[items,ax]=plot_phase_panel(ax,qBaseline,c,baselineComparison,true);
set(items(1),'XData',q20.t,'YData',q20.Dtheta, ...
    'DisplayName','$D_{\theta}^{(20\omega)}(t)$');
lg=findobj(h,'Type','Legend');
if ~isempty(lg)
    names=lg.String;
    names{1}='$D_{\theta}^{(20\omega)}(t)$';
    lg.String=names;
end
visible=q20.Dtheta(q20.t <= min(2,q20.t(end)));
if ~isempty(visible) && max(visible) > ax.YLim(2)
    ylim(ax,[0 1.08*max(visible)]);
end
text(ax,0.02,0.98,'(a)','Units','normalized', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'VerticalAlignment','top','HandleVisibility','off');
export_figure(h,outdir,'10_phase_omega20_blue_only');
end

function plot_pi_snapshots(q,s,cfg,outdir)
times = [0 10];
assert(q.t(1)<=times(1) && q.t(end)>=times(end), ...
    'Baseline trajectory must cover both t=0 and t=10.');

h = newfigure(cfg,'two_panel');
set(h,'Position',[1 1 cfg.FigureWidth 2.05]);
tiles = tiledlayout(h,1,2,'Padding','compact','TileSpacing','compact');
for k = 1:numel(times)
    ax = nexttile(tiles,k);
    draw_pi_snapshot(ax,q,s,times(k),true,k);
end
export_figure(h,outdir,'00_pi_snapshots_onecolumn');

for k = 1:numel(times)
    h = figure('Units','inches','Position',[1 1 1.68 1.88], ...
        'Color','w','Visible',cfg.Visible,'Renderer','painters');
    ax = axes('Parent',h,'Position',[0.07 0.09 0.86 0.78]);
    draw_pi_snapshot(ax,q,s,times(k),true,k);
    export_figure(h,outdir,sprintf('fig1_snapshot_t%d',times(k)));
end
end

function draw_pi_snapshot(ax,q,s,tk,withPanelLetter,k)
idx = find(abs(q.t-tk)<=1e-10,1,'first');
if isempty(idx)
    x = interp1(q.t,q.thetaTilde,tk,'pchip').';
else
    x = q.thetaTilde(idx,:).';
end
theta = x + (1-s(:))*pi/2 + q.omegaReference*tk;

red  = [0.82 0.11 0.12];
blue = [0.07 0.29 0.72];
phi = linspace(0,2*pi,361);
hold(ax,'on');
plot(ax,cos(phi),sin(phi),'k--','LineWidth',0.9, ...
    'HandleVisibility','off');
scatter(ax,cos(theta(s>0)),sin(theta(s>0)),38,red,'filled', ...
    'MarkerEdgeColor',red,'LineWidth',0.35, ...
    'HandleVisibility','off');
scatter(ax,cos(theta(s<0)),sin(theta(s<0)),38,blue,'filled', ...
    'MarkerEdgeColor',blue,'LineWidth',0.35, ...
    'HandleVisibility','off');
plot(ax,0,0,'ko','MarkerSize',3.5,'MarkerFaceColor','k', ...
    'HandleVisibility','off');
axis(ax,'equal');
xlim(ax,[-1.20 1.20]);
ylim(ax,[-1.20 1.20]);
axis(ax,'off');
if withPanelLetter
    letters = {'(a)','(b)'};
    heading = [letters{k} ' $t=' num2str(tk) '$'];
else
    heading = ['$t=' num2str(tk) '$'];
end
title(ax,heading,'Interpreter','latex', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','normal');
end

function plot_unrestricted_initial(q,P,M,cfg,outdir)
blue=[0 0.4470 0.6980]; gray=[0.42 0.42 0.42];
Cnet=NaN;
if all(P>0)
    ratio=max(M./P);
    if ratio>=0 && ratio<1
        Cnet=2*acos(ratio);
    end
end

h=newfigure(cfg); ax=gca;
hold(ax,'on'); box(ax,'on'); grid(ax,'off');
p1=plot(ax,q.t,q.Dtheta,'-','Color',blue,'LineWidth',1.35, ...
    'DisplayName','$D_{\theta}(t)$');
items=p1;
if isfinite(Cnet)
    p2=yline(ax,Cnet,'--','Color',gray,'LineWidth',1.1, ...
        'DisplayName','$n {C}_{\mathrm{net}}$');
    items=[p1 p2];
end
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$D_{\theta}(t)$','Interpreter','latex');
xlim(ax,[0 q.t(end)]);
lg=legend(ax,items,'Interpreter','latex','Location','northeast', ...
    'AutoUpdate','off','Box','on');
set(lg,'FontSize',8,'Color','white');
text(ax,0.02,0.15,'(a)','Units','normalized', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'VerticalAlignment','top','HandleVisibility','off');
export_figure(h,outdir,'11_unrestricted_phase_diameter');

h=newfigure(cfg); ax=gca;
hold(ax,'on'); box(ax,'on'); grid(ax,'off');
threshold=max(1e-10,10*cfg.AbsTol);
observed=q.Dfrequency;
observed(observed<threshold)=NaN;
p1=semilogy(ax,q.t,observed,'-','Color',blue,'LineWidth',1.35, ...
    'DisplayName','$D(t)$');
p2=yline(ax,cfg.FrequencyTol,'--','Color',gray,'LineWidth',1.1, ...
    'DisplayName','Numerical tolerance');
set(ax,'YScale','log');
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$D(t)$','Interpreter','latex');
xlim(ax,[0 q.t(end)]);
lg=legend(ax,[p1 p2],'Interpreter','latex','Location','northeast', ...
    'AutoUpdate','off','Box','on');
set(lg,'FontSize',8,'Color','white');
text(ax,0.02,0.15,'(b)','Units','normalized', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'VerticalAlignment','top','HandleVisibility','off');
export_figure(h,outdir,'12_unrestricted_frequency_diameter');
end

function plot_modified_network(qBaseline,qModified,cfg,outdir)
blue = [0 0.4470 0.6980];
orange = [0.8350 0.3690 0];
tmax = min([10,qBaseline.t(end),qModified.t(end)]);

h = newfigure(cfg); ax = gca;
p1 = plot(ax,qBaseline.t,qBaseline.Dtheta,'-', ...
    'Color',blue,'LineWidth',1.35,'DisplayName','Baseline network');
hold(ax,'on');
p2 = plot(ax,qModified.t,qModified.Dtheta,'--', ...
    'Color',orange,'LineWidth',1.25,'DisplayName','Modified network');
box(ax,'on'); grid(ax,'off'); xlim(ax,[0 tmax]);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$D_{\theta}(t)$','Interpreter','latex');
lg = legend(ax,[p1 p2],'Interpreter','latex','Location','northwest', ...
    'AutoUpdate','off','Box','on');
set(lg,'FontSize',8,'Color','white');
text(ax,0.98,0.98,'(a)','Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'HandleVisibility','off');
export_figure(h,outdir,'13_modified_network_phase');

h = newfigure(cfg); ax = gca;
d0=qBaseline.Dfrequency;
d1=qModified.Dfrequency;
p1=plot(ax,qBaseline.t,d0,'-', ...
    'Color',blue,'LineWidth',1.35,'DisplayName','Baseline network');
hold(ax,'on');
p2=plot(ax,qModified.t,d1,'--', ...
    'Color',orange,'LineWidth',1.25,'DisplayName','Modified network');
box(ax,'on'); grid(ax,'off'); xlim(ax,[0 tmax]);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$D(t)$','Interpreter','latex');
lg=legend(ax,[p1 p2],'Interpreter','latex','Location','north', ...
    'AutoUpdate','off','Box','on');
set(lg,'FontSize',8,'Color','white');
text(ax,0.98,0.98,'(b)','Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'HandleVisibility','off');
export_figure(h,outdir,'14_modified_network_frequency');
end

function plot_gain(rows,cfg,outdir)
ks = unique([rows.K]); mat = NaN(numel(ks),5);
s = [ones(7,1);-ones(3,1)]; [P,M] = pair_margins(baseline_matrix(10,s),s);
uniformRows = struct([]);
for j = 1:numel(ks)
    z = rows([rows.K]==ks(j)); good = z([z.FrequencyTailPass]);
    cb = certificate(P,M,10,ks(j),2*pi/3,1,[],[]);
    upper = cb.Cinf;
    if cb.FrequencyCertified, upper = cb.CinfImproved; end
    cb.K = ks(j); uniformRows = append(uniformRows,cb);
    mat(j,:) = [finite_median([good.FinalDtheta]),upper, ...
        cb.Rate,finite_median([good.RateFitted]), ...
        finite_median([good.SpectralGapFinal])];
end
write_rows(outdir,'gain_uniform_certificate',uniformRows);
blue=[0 0.4470 0.6980]; orange=[0.8350 0.3690 0];
purple=[0.4940 0.1840 0.5560];
h = newfigure(cfg);
p1 = loglog(ks,mat(:,1),'-','Color',blue,'LineWidth',1.35, ...
    'Marker','none', ...
    'DisplayName','$\mathrm{Median}\ D_{\theta}(t_{\mathrm{end}})$');
set(gca,'NextPlot','add');
p2 = loglog(ks,mat(:,2),'--','Color',orange,'LineWidth',1.25, ...
    'Marker','none','DisplayName','Theoretical asymptotic bound');
set(gca,'XScale','log','YScale','log');
xlabel('$K$','Interpreter','latex'); ylabel('$D_{\theta}$','Interpreter','latex');
lg=manuscript_legend([p1 p2],'southwest'); lg.NumColumns=1;
export_figure(h,outdir,'03_gain_dispersion');
h = newfigure(cfg);
p1 = plot(ks,mat(:,3),'--','Color',orange,'LineWidth',1.25, ...
    'Marker','none','DisplayName','Guaranteed rate');
set(gca,'NextPlot','add');
p2 = plot(ks,mat(:,4),'-','Color',blue,'LineWidth',1.35, ...
    'Marker','none','DisplayName','Fitted decay rate');
p3 = plot(ks,mat(:,5),'-.','Color',purple,'LineWidth',1.15, ...
    'Marker','none','DisplayName','Spectral rate estimate');
xlabel('$K$','Interpreter','latex'); ylabel('Rate','Interpreter','latex');
lg=manuscript_legend([p1 p2 p3],'northwest'); lg.NumColumns=1;
export_figure(h,outdir,'04_gain_rates');
end

function plot_sparse(rows,cfg,outdir)
h = newfigure(cfg); ax = gca;
colors = [0 0.4470 0.6980; 0.8350 0.3690 0];
markers = {'o','s'};
items = gobjects(1,4);
for ig = 1:2
    z = rows([rows.Trial]==ig);
    [x,order] = sort([z.N]); z = z(order);
    if ig==1, name = 'N'; else, name = '10'; end
    items(2*ig-1) = plot(ax,x,[z.RateGuaranteed],[markers{ig} '--'], ...
        'Color',colors(ig,:),'LineWidth',1.15,'MarkerSize',5, ...
        'MarkerFaceColor','w', ...
        'DisplayName',['$K=' name ',\ \mathrm{guaranteed}$']);
    items(2*ig) = plot(ax,x,[z.SpectralGapFinal],[markers{ig} '-'], ...
        'Color',colors(ig,:),'LineWidth',1.40,'MarkerSize',5, ...
        'MarkerFaceColor','w', ...
        'DisplayName',['$K=' name ',\ \mathrm{spectral}$']);
end
set(ax,'TickLabelInterpreter','latex');
xlabel(ax,'$N$','Interpreter','latex');
ylabel(ax,'Rate','Interpreter','latex');
lg = manuscript_legend(items,'best');
lg.NumColumns = 2;
export_figure(h,outdir,'05_sparse_minority');
end

function plot_topology(rows,cfg,outdir)
h = newfigure(cfg); ax = gca;
families = {'ER','WS','BA'};
colors = [0 0.4470 0.6980; 0.8350 0.3690 0; 0 0.6200 0.4510];
markers = {'o','s','^'};
items = gobjects(1,6);
n = max([rows.N]);
for j = 1:3
    z = rows(strcmp({rows.Topology},families{j}) & [rows.N]==n);
    [x,order] = sort([z.MeanDensity]); z = z(order);
    p = [z.FrequencyFraction];
    pc = [z.CertifiedFraction];
    items(2*j-1) = errorbar(ax,x,p,p-[z.FrequencyCILow], ...
        [z.FrequencyCIHigh]-p,['-' markers{j}], ...
        'Color',colors(j,:),'LineWidth',1.35,'MarkerSize',5, ...
        'CapSize',3,'MarkerFaceColor','w', ...
        'DisplayName',['$\mathrm{' families{j} '}$: observed']);
    items(2*j) = errorbar(ax,x,pc,pc-[z.CertCILow], ...
        [z.CertCIHigh]-pc,['--' markers{j}], ...
        'Color',colors(j,:),'LineWidth',1.10,'MarkerSize',5, ...
        'CapSize',3,'MarkerFaceColor','w', ...
        'DisplayName',['$\mathrm{' families{j} '}$: certified']);
end
set(ax,'TickLabelInterpreter','latex');
xlabel(ax,'Mean density','Interpreter','latex');
ylabel(ax,'Fraction','Interpreter','latex');
ylim(ax,[-0.03 1.05]);
lg = manuscript_legend(items,'best');
lg.NumColumns = 2;
export_figure(h,outdir,'06_topology');
end

function plot_degree(rows,cfg,outdir)
h = newfigure(cfg); ax = gca;
items = gobjects(1,2);
for ig = 0:1
    z = rows([rows.GainNormalized]==logical(ig));
    [x,order] = sort([z.N]); z = z(order);
    p = [z.FrequencyFraction];
    if ig==0
        name = '$K=10$'; fmt = 'o-'; col = [0 0.4470 0.6980];
    else
        name = '$K\langle k\rangle/N=10$';
        fmt = 's--'; col = [0.8350 0.3690 0];
    end
    items(ig+1) = errorbar(ax,x,p,p-[z.FrequencyCILow], ...
        [z.FrequencyCIHigh]-p,fmt, ...
        'Color',col,'LineWidth',1.35-(ig*0.20), ...
        'MarkerSize',5,'CapSize',4,'MarkerFaceColor','w', ...
        'DisplayName',name);
end
set(ax,'TickLabelInterpreter','latex');
xlabel(ax,'$N$','Interpreter','latex');
ylabel(ax,'Fraction','Interpreter','latex');
ylim(ax,[-0.03 1.05]);
lg = manuscript_legend(items,'best');
lg.NumColumns = 1;
export_figure(h,outdir,'07_fixed_degree');
end

function plot_map(rows,cfg,outdir)
a = NaN(numel(cfg.Mismatches),numel(cfg.MapGains)); b = a;
for j = 1:size(a,1)
    for k = 1:size(a,2)
        id = (j-1)*size(a,2)+k; z = rows([rows.CaseID]==id);
        a(j,k) = mean([z.FrequencyTailPass]);
        b(j,k) = mean([z.Certified]);
    end
end
cmap = parula(256);
for k = 1:2
    if k==1
        values=a; name='08_map_frequency';
    else
        values=b; name='09_map_certificate';
    end
    h = newfigure(cfg); ax = gca;
    nx = size(values,2); ny = size(values,1);
    colormap(ax,cmap); clim(ax,[0 1]);
    for j = 1:ny
        for l = 1:nx
            v = values(j,l);
            if isfinite(v)
                index = 1+round(255*min(1,max(0,v)));
                face = cmap(index,:);
                label = sprintf('%.2f',v);
                luminance = face*[0.2126;0.7152;0.0722];
                if luminance < 0.48
                    txt = [1 1 1];
                else
                    txt = [0 0 0];
                end
            else
                face = [0.90 0.90 0.90]; label = '--';
                txt = [0.25 0.25 0.25];
            end
            patch(ax,[l-.5 l+.5 l+.5 l-.5], ...
                [j-.5 j-.5 j+.5 j+.5],face, ...
                'EdgeColor',[1 1 1],'LineWidth',0.45, ...
                'HandleVisibility','off');
            text(ax,l,j,label,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','Color',txt, ...
                'FontName','Times New Roman','FontSize',7.5, ...
                'Interpreter','latex','HandleVisibility','off');
        end
    end
    set(ax,'YDir','normal','XTick',1:nx, ...
        'XTickLabel',cellstr(num2str(cfg.MapGains(:))), ...
        'YTick',1:ny, ...
        'YTickLabel',cellstr(num2str(cfg.Mismatches(:))), ...
        'TickLabelInterpreter','latex');
    xlim(ax,[0.5 nx+0.5]); ylim(ax,[0.5 ny+0.5]);
    xlabel(ax,'$K$','Interpreter','latex');
    ylabel(ax,'$\Delta\omega$','Interpreter','latex');
    cb = colorbar(ax,'eastoutside');
    set(cb,'FontName','Times New Roman','FontSize',8, ...
        'TickLabelInterpreter','latex','Ticks',0:0.25:1);
    cb.Label.String = 'Fraction';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 8;
    export_figure(h,outdir,name);
end
end

function rows = convergence_check(cfg)
rows = struct([]);
N = 10; s = [ones(7,1);-ones(3,1)]; A = baseline_matrix(N,s);
x0 = linspace(-pi/3,pi/3,N)'; omega = linspace(-.5,.5,N)';
for j = 1:3
    Ks = [10 2 10]; scales = [1 1 20]; K = Ks(j); om = scales(j)*omega;
    B = sparse(A)*spdiags(s,0,N,N); ref = mean(om);
    fun = @(t,x) rhs_transformed(x,B,om-ref,K/N);
    t = (0:.01:10)';
    opts = odeset('RelTol',1e-11,'AbsTol',1e-13,'MaxStep',.01);
    [~,X] = ode45(fun,t,x0,opts);
    theta0 = x0+(1-s)*pi/2;
    [~,Y] = ode45(@(tt,y)rhs_original(y,A,s,om,K),t,theta0,opts);
    Y = bsxfun(@minus,Y,((1-s)*pi/2)');
    XR = bsxfun(@minus,X,X(:,1)); YR = bsxfun(@minus,Y,Y(:,1));
    Z1 = rk4(fun,x0,.01,10); Z2 = rk4(fun,x0,.005,10); Z2 = Z2(1:2:end,:);
    Z1 = bsxfun(@minus,Z1,Z1(:,1)); Z2 = bsxfun(@minus,Z2,Z2(:,1));
    r = struct('Case',j,'K',K,'FrequencyScale',scales(j), ...
        'OriginalVsTransformed',max(abs(XR(:)-YR(:))), ...
        'RK4_h_vs_h2',max(abs(Z1(:)-Z2(:))), ...
        'RK4_h_vs_ode45',max(abs(Z1(:)-XR(:))), ...
        'RK4_h2_vs_ode45',max(abs(Z2(:)-XR(:))));
    assert(r.OriginalVsTransformed<2e-6,'Transformation check failed.');
    rows = append(rows,r);
end
assert(cfg.RelTol>0);
end

function X = rk4(fun,x,h,T)
steps = round(T/h); X = zeros(steps+1,numel(x)); X(1,:) = x';
for k = 1:steps
    t = (k-1)*h; a = fun(t,x); b = fun(t+h/2,x+h*a/2);
    c = fun(t+h/2,x+h*b/2); d = fun(t+h,x+h*c);
    x = x+h*(a+2*b+2*c+d)/6; X(k+1,:) = x';
end
end

function tests = selftest()
N = 10; s = [ones(7,1);-ones(3,1)]; A = ones(N)-eye(N);
[p,m] = pair_margins(A,s);
assert(max(abs(p-7))<1e-12 && max(abs(m-3))<1e-12);
assert(abs(spectral_gap(A,s,zeros(N,1),10)-4)<1e-10);
A = baseline_matrix(N,s); [p,m] = pair_margins(A,s);
c = certificate(p,m,N,10,2*pi/3,1,2.2,.75);
assert(abs(c.Gamma-1.4418857633120108)<1e-10);
assert(abs(c.Lambda-2.8572612829996995)<1e-10);
assert(abs(c.EntryBound-23.901933686984094)<1e-7);
assert(abs(c.Cinf-.7082464128595691)<1e-10);
pd = []; md = [];
for j = 2:N
    for i = 1:j-1
        z = setdiff(1:N,[i j]); pos = z(s(z)>0); neg = z(s(z)<0);
        pd(end+1,1) = sum(min(A(i,pos),A(j,pos)))+A(i,j)*(s(j)>0)+A(j,i)*(s(i)>0);
        md(end+1,1) = sum(max(A(i,neg),A(j,neg)))+A(i,j)*(s(j)<0)+A(j,i)*(s(i)<0);
    end
end
assert(max(abs(p-pd))<1e-12 && max(abs(m-md))<1e-12);
x = linspace(-1,1,N)'; om = linspace(-.5,.5,N)'; theta = x+(1-s)*pi/2;
f = rhs_original(theta,A,s,om,10);
g = rhs_transformed(x,A*diag(s),om,1);
assert(max(abs(f-g))<1e-12);
assert(abs(circular_width([0 2*pi .1])-.1)<1e-12);
AA = broadcast_matrix(N,1,N,.9,.2); [p,m] = pair_margins(AA,s);
assert(max(abs(p-.9))<1e-12 && max(abs(m-.2))<1e-12);
tests = struct('Passed',true,'Gamma22',c.Gamma,'Lambda075',c.Lambda, ...
    'Cinf',c.Cinf,'EntryBound',c.EntryBound,'RHSdifference',max(abs(f-g)));
end
