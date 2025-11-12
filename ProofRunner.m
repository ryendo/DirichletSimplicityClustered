classdef ProofRunner < handle
% A runner for Algorithm 1 & 2.
%
% This class can:
%   • Initialize INTLAB & repository paths (my_intlab_config_alone / my_intlab_config)
%   • Run Algorithm 1
%       - Single direction, direction-interval, or full sweep
%       - Clear "what is being computed" message
%       - Estimated time/step count and live percent-complete
%       - Resume from existing CSV (skips completed bins)
%   • Prepare Algorithm 2 (runs prep.sh)
%   • Run Algorithm 2 (runs main_algo2.sh), optionally launch+monitor
%       - Live progress by reading prep/algo2_list_j.csv
%       - ETA estimation based on observed completion rate
%   • Point/box checks of λ₂ upper / λ₃ lower bounds near equilateral
%
%
% QUICK START
%   s = ProofRunner;
%   s.setupAll();                      % INTLAB + prep.sh
%   s.runAlgo1All();                   % full Ω_up sweep with ETA/progress
%   s.summarizeAlgo1CSV();             % check sup(mu1) < inf(mu2)
%
%   s.runAlgo1Direction(I_pi/12);        % single direction near equilateral
%   s.runAlgo1Interval([I_pi/20,I_pi/20+intval('1e-5')]); % small interval
%
%   s.boundsAtPoint(intval('1')/2,sqrt(intval('3'))/2-intval('1e-3'));       % eigenvalue bounds at (s,t)
%   s.boundsOnBox(intval('1')/2,intval('1')/2,sqrt(intval('3'))/2-2*intval('1e-3'),sqrt(intval('3'))/2-intval('1e-3'));     %[s_inf,s_sup,t_inf,t_sup]
%
%   s.prepareAlgo2();                  % runs ./prep.sh
%   s.runAlgo2('detached',true);       % launch main_algo2.sh (Unix)
%   s.monitorAlgo2();                  % live progress from CSV
%
% OUTPUT (Algorithm 1)
%   - Prints: what is computed, planned steps, warm-up timing
%   - During run: k/Total, % complete, ETA
%
% DEPENDENCIES (already in your repo)
%   - INTLAB; verifyeig; calc_eigen_bounds_any_order_for_quotients
%   - calc_eigen_bounds_any_order (for point/box checks)
%   - prep.sh, main_algo2.sh, prep/algo2_list_j.csv (Algorithm 2)
%
% Author: Ryoki Endo and Xuefeng Liu

properties
    % ===== Algorithm 1 global parameters =====
    omega_N (1,1) double = 1000      % bins over [0, π/3]
    mesh_N  (1,1) double = 32        % FEM mesh parameter
    mesh_LG (1,1) double = 8         % FEM mesh parameter for L-G method
    ord     (1,1) double = 5         % FEM polynomial order
    ep      (1,1) string = '1e-5'    % t ∈ [0, ep]

    % Output CSV (if empty → auto timestamp in results/)
    outFile (1,:) char = ''

    % Behavior
    verbose (1,1) logical = true
    resume  (1,1) logical = true

    % Bounds callback for point/box evidence
    boundsFcn function_handle = @ProofRunner.default_boundsFcn

    % ===== Algorithm 2 paths =====
    prepScript   (1,:) char = 'prep.sh'
    mainScript2  (1,:) char = 'main_algo2.sh'
    listCsvPath  (1,:) char = 'prep/algo2_list_j.csv'
    lockPath     (1,:) char = 'prep/algo2_list_j.lock'

    % Monitor settings
    pollSeconds (1,1) double = 5
end

methods
    function self = ProofRunner(varargin)
        % Name-value pairs
        if mod(nargin,2) ~= 0
            error('Use name-value pairs, e.g., ProofRunner(''omega_N'',1000).');
        end
        for k = 1:2:nargin
            name = char(varargin{k});
            if ~isprop(self,name), error('Unknown option: %s',name); end
            self.(name) = varargin{k+1};
        end
        if ~exist('results','dir'), mkdir('results'); end
        if isempty(self.outFile)
            ts = datestr(now,'yyyy-mm-dd_HH-MM-SS');
            self.outFile = fullfile('results',sprintf('quotients_%s.csv',ts));
        end
    end

    %==================== Setup (INTLAB + prep.sh) ====================%
    function setupAll(self)
        % Full one-shot setup: INTLAB + prep.sh
        self.setupIntlab();
        self.prepareAlgo2();
    end

    function setupIntlab(self)
        % Initialize INTLAB and repo paths by calling your existing config.
        try
            if exist('my_intlab_config_alone','file')==2
                my_intlab_config_alone;
            elseif exist('my_intlab_config','file')==2
                my_intlab_config;
            else
                warning(['No my_intlab_config_alone/my_intlab_config found. ', ...
                         'Please ensure INTLAB and repo paths are on MATLAB path.']);
            end
        catch ME
            warning('INTLAB setup call failed: %s', ME.message);
        end
        format compact short infsup
    end

    function prepareAlgo2(self)
        % Run ./prep.sh as-is to prepare Algorithm 2 task list.
        if ~isunix
            error('Algorithm 2 preparation requires a Unix-like shell.');
        end
        self.assertFile(self.prepScript,'prep.sh not found.');
        cmd = sprintf('bash %s', self.prepScript);
        if self.verbose, fprintf('Running: %s\n', cmd); end
        st = system(cmd);
        if st ~= 0
            error('prep.sh failed (exit code %d). Check its console output.', st);
        end
        if ~isfile(self.listCsvPath)
            error('Expected %s not found after prep.', self.listCsvPath);
        end
        if self.verbose
            fprintf('Preparation complete. Task list: %s\n', self.listCsvPath);
        end
    end

    %==================== Algorithm 1 (Ω_up) ====================%
    function runAlgo1All(self),  self.runAlgo1Interval([self.I_intval(0), self.I_intval('pi')/3], self.omega_N); end
    function runAlgo1Direction(self, delta), self.runAlgo1Interval([delta, delta], 1); end

    function runAlgo1Interval(self, interval, K)
        % Detailed teacher-facing preface
        self.setupIntlab();
        self.printAlgo1Preface(interval, K);

        % Matrices & basis (exactly as your main_algo1.m does internally)
        [mat,u1,u2,u3] = self.prepareFEM();

        % Partition interval
        if nargin<3 || isempty(K)
            K = self.inferBinCount(interval, self.omega_N);
        end
        edges   = linspace(interval(1), interval(2), K+1)
        idxList = 1:K;

        % Resume mask
        doneMask = false(size(idxList));
        if self.resume && isfile(self.outFile)
            have = self.readResultsIndices(self.outFile);
            doneMask(ismember(idxList, have)) = true;
            if self.verbose
                fprintf('Resume: %d/%d bins already present → skipped.\n', sum(doneMask), K);
            end
        end
        todo = idxList(~doneMask);
        total = numel(todo);
        if total == 0
            if self.verbose
                fprintf('Nothing to do; interval fully present in %s\n', self.outFile);
            end
            return;
        end

        % Ensure header
        self.ensureCsvHeader(self.outFile);

        % Warm-up for ETA
        step_t = NaN;
        if self.verbose && total>=1
            d0 = [edges(todo(1)), edges(todo(1)+1)];
            t0 = tic;
            try
                [~,~] = self.quotient_calc(self.interval_t(), self.interval_delta(d0), mat, u1, u2, u3);
                step_t = toc(t0);
                if isfinite(step_t) && step_t>0
                    fprintf('Planned steps: %d | warm-up per-step ≈ %.2fs | ETA ≈ %s\n', ...
                        total, step_t, self.formatDuration(step_t*total));
                end
            catch ME
                warning('Warm-up failed (%s). ETA will be rough.', ME.message);
            end
        end

        % Main loop with live percent
        t_start = tic;
        if self.verbose
            fprintf('--- Running Algorithm 1 on %d bin(s) ---\n', total);
        end
        for kk = 1:total
            i = todo(kk);
            d = [edges(i), edges(i+1)];
            try
                [~, interval_mu] = self.quotient_calc(self.interval_t(), self.interval_delta(d), mat, u1, u2, u3);
                row = [i, self.callI_inf(interval_mu(1)), self.callI_sup(interval_mu(1)), ...
                          self.callI_inf(interval_mu(2)), self.callI_sup(interval_mu(2))];
                writematrix(row, self.outFile, 'WriteMode','append');
            catch ME
                warning('Bin %d failed: %s', i, ME.message);
                row = [i, NaN, NaN, NaN, NaN];
                writematrix(row, self.outFile, 'WriteMode','append');
            end

            if self.verbose
                pct = 100*kk/total;
                elapsed = toc(t_start);
                eta_sec = (elapsed/max(kk,1))*(total-kk);
                fprintf('\r%s %d/%d (%.1f%%) ETA %s', ...
                        self.progressBar(kk,total,32), kk,total,pct, self.formatDuration(eta_sec));
                if kk==total, fprintf('\n'); end
            end
        end
        if self.verbose
            fprintf('Done. Results → %s\n', self.outFile);
        end
    end

    function summarizeAlgo1CSV(self, file)
        if nargin<2 || isempty(file), file = self.outFile; end
        if ~isfile(file), error('CSV not found: %s', file); end

        T = [];
        try
            Tt = readtable(file);
            if all(ismember({'idx','inf_mu1','sup_mu1','inf_mu2','sup_mu2'}, Tt.Properties.VariableNames))
                T = [Tt.idx, Tt.inf_mu1, Tt.sup_mu1, Tt.inf_mu2, Tt.sup_mu2];
            end
        catch, end
        if isempty(T)
            T = readmatrix(file);
            if any(isnan(T(1,:))), T = T(2:end,:); end
        end
        if size(T,2) < 5, error('Unexpected CSV format: %s', file); end

        sup_mu1 = T(:,3);
        inf_mu2 = T(:,4); % Column order: idx, inf1, sup1, inf2, sup2
        ok = sup_mu1 < inf_mu2;

        fprintf('File: %s\n', file);
        fprintf('Bins: %d | OK: %d | NG: %d\n', numel(ok), sum(ok), sum(~ok));
        if all(ok)
            fprintf('>>> PROOF SUCCESSFUL: sup(mu1) < inf(mu2) holds for all bins.\n');
        else
            bad = find(~ok);
            fprintf('>>> PROOF FAILED in bins: ');
            fprintf('%d ', bad);
            fprintf('\n');
        end
    end

    %==================== Algorithm 2 (Ω_down) ====================%
    function runAlgo2(self, varargin)
        % Run main_algo2.sh as-is. Options:
        %   'detached' (logical, default=false): if true, launch with nohup & (Unix)
        %   'logFile'  (char): optional shell redirection target for detached mode
        p = inputParser;
        addParameter(p,'detached',false,@islogical);
        addParameter(p,'logFile','Each_Process/log/main_algo2_launcher.log',@(s)ischar(s)||isstring(s));
        parse(p, varargin{:});
        opt = p.Results;

        if ~isunix
            error('Algorithm 2 requires a Unix-like shell. Please use bash.');
        end
        self.assertFile(self.mainScript2, 'main_algo2.sh not found.');
        if opt.detached
            % Launch in background; user can monitor via monitorAlgo2()
            logFile = char(opt.logFile);
            [logDir,~] = fileparts(logFile);
            if ~isempty(logDir) && ~exist(logDir,'dir'), mkdir(logDir); end
            cmd = sprintf('nohup bash %s > %s 2>&1 &', self.mainScript2, logFile);
            if self.verbose, fprintf('Launching (detached): %s\n', cmd); end
            st = system(cmd);
            if st ~= 0
                error('Failed to launch main_algo2.sh (exit %d).', st);
            end
            if self.verbose
                fprintf('main_algo2.sh started. Tail logs or call monitorAlgo2() for progress.\n');
            end
            else
            cmd = sprintf('bash %s', self.mainScript2);
            if self.verbose, fprintf('Running: %s\n', cmd); end
            st = system(cmd);
            if st ~= 0
                error('main_algo2.sh failed (exit %d). Check Each_Process/log/*.log', st);
            end
            if self.verbose, fprintf('Algorithm 2 completed.\n'); end
        end
    end

    function monitorAlgo2(self, varargin)
        % Live progress for Algorithm 2 by reading prep/algo2_list_j.csv
        % Options:
        %   'csv' (char)         : override CSV path
        %   'pollSeconds' (double): override polling interval
        p = inputParser;
        addParameter(p,'csv', self.listCsvPath, @(s)ischar(s)||isstring(s));
        addParameter(p,'pollSeconds', self.pollSeconds, @isnumeric);
        parse(p, varargin{:});
        csvPath = char(p.Results.csv);
        poll    = p.Results.pollSeconds;

        self.assertFile(csvPath, 'CSV list not found. Run prepareAlgo2() first.');
        fprintf('--- Algorithm 2 monitor ---\n');
        fprintf('File: %s (poll=%.1fs)\n', csvPath, poll);

        % Heuristic ETA using completion rate
        donePrev = 0; t0 = tic; tLast = tic;
        while true
            [n0,n2,n1,tot] = self.readAlgo2States(csvPath);
            done = n1;
            pct  = 100 * done / max(tot,1);

            % rate & eta
            dt  = toc(tLast); tLast = tic;
            dN  = max(done - donePrev, 0);
            rate = dN / max(dt, eps);        % jobs/sec
            rem  = max(tot - done, 0);
            eta  = rem / max(rate, eps);     % sec
            donePrev = done;

            fprintf('\rTotal=%d | pending=%d | running=%d | done=%d | %.1f%% | ETA %s', ...
                tot, n0, n2, n1, pct, self.formatDuration(eta));

            if done >= tot
                fprintf('\nAll tasks finished. Elapsed %s.\n', self.formatDuration(toc(t0)));
                break;
            end
            pause(poll);
        end
    end

    %==================== Partial evidence near equilateral ====================%
    function [lam2, lam3] = boundsAtPoint(self, a, b, varargin)
        p = inputParser;
        addParameter(p,'mesh_N', self.mesh_N);
        addParameter(p,'N_LG',   self.mesh_LG);
        addParameter(p,'N_rho',   self.mesh_N);
        addParameter(p,'ord',    self.ord);
        parse(p, varargin{:}); q = p.Results;
        isLG = 1;

        tri = [self.callI_intval(0),self.callI_intval(0); self.callI_intval(1),self.callI_intval(0); a, b];
        lams = self.boundsFcn(tri, q.N_LG,q.N_rho,q.ord,isLG);
        lam2 = lams(1);
        lam3 = lams(2);

        if self.verbose
            fprintf('Point (a=[%.6g,%.6g], b=[%.6g,%.6g]): \n λ2 ∈ [%.12g,%.12g]  λ3  ∈ [%.12g,%.12g],\n  gap ≥ %.12g\n', self.callI_inf(a),self.callI_sup(a), self.callI_inf(b),self.callI_sup(b), self.callI_inf(lam2), self.callI_sup(lam2), self.callI_inf(lam3), self.callI_sup(lam3), self.callI_inf(lam3 - lam2));
        end
    end

    function [up2, lo3] = boundsOnBox(self, a_lo, a_hi, b_lo, b_hi, varargin)
        p = inputParser;
        addParameter(p,'mesh_N', self.mesh_N);
        addParameter(p,'N_LG',   self.mesh_LG);
        addParameter(p,'N_rho',   self.mesh_N);
        addParameter(p,'ord',    self.ord);
        parse(p, varargin{:}); q = p.Results;
        isLG = 1;

        lams1 = self.boundsFcn([self.callI_intval('0'),self.callI_intval('0');self.callI_intval('1'),self.callI_intval('0');a_lo, b_lo], q.N_LG,q.N_rho,q.ord,isLG);
        lams2 = self.boundsFcn([self.callI_intval('0'),self.callI_intval('0');self.callI_intval('1'),self.callI_intval('0');a_hi, b_hi], q.N_LG,q.N_rho,q.ord,isLG);
        up2 = self.callI_sup(lams1(1));
        lo3 = self.callI_inf(lams2(2));
        
        if self.verbose
            fprintf('Box [%.6g,%.6g]×[%.6g,%.6g]:  sup λ2 ≤ %.12g,  inf λ3 ≥ %.12g,  min gap ≥ %.12g\n', ...
                self.callI_inf(a_lo), self.callI_sup(a_hi), self.callI_inf(b_lo), self.callI_sup(b_hi), up2, lo3, self.callI_inf(self.callI_intval(lo3) - self.callI_intval(up2)));
        end
    end

    function explain(self)
        % One-shot explanation for advisors/users
        fprintf([ ...
            "ALGORITHM 1 (Ω_up):\n" ...
            "- Discretize δ∈[0,π/3] into bins. For each bin (an angular interval),\n" ...
            "  we compute rigorous interval bounds (μ₁, μ₂) for the difference quotients\n" ...
            "  of λ₂ and λ₃ using verified FEM and interval arithmetic.\n" ...
            "- Output CSV columns: idx, inf_mu1, sup_mu1, inf_mu2, sup_mu2.\n" ...
            "- PROOF condition (per bin): sup(μ₁) < inf(μ₂).\n" ...
            "- During run: we display steps, percent complete, and ETA.\n\n" ...
            "ALGORITHM 2 (Ω_down):\n" ...
            "- The domain is partitioned into subregions. Each task j computes\n" ...
            "  verified bounds for λ₂ and λ₃ on its subregion.\n" ...
            "- We run your unmodified prep.sh and main_algo2.sh, and monitor progress\n" ...
            "  from prep/algo2_list_j.csv (pending/running/done) with ETA.\n\n" ...
            "PARTIAL EVIDENCE:\n" ...
            "- boundsAtPoint(a,b): compute λ₂ upper and λ₃ lower bounds at the vertex\n" ...
            "  (0.5+a, sqrt(3)/2+b). boundsOnBox(...) checks corners of a given box.\n" ...
        ]);
    end
end

%==================== Private: Algorithm 1 core (unchanged math) ====================%
methods (Access = private)
    function printAlgo1Preface(self, interval, K)
        fprintf('--- Algorithm 1: rigorous difference-quotient bounds (Ω_up) ---\n');
        fprintf('What is computed:\n');
        fprintf('  - Verified interval bounds for the difference quotients of λ₂ and λ₃\n');
        fprintf('    over a direction interval (near equilateral triangle).\n');
        fprintf('Parameters:\n');
        fprintf('  mesh_N=%d, ord=%d, ep=%g, omega_N=%d\n', self.mesh_N, self.ord, self.ep, self.omega_N);
        fprintf('Direction interval: [%.9f, %.9f] rad\n', self.callI_inf(interval), self.callI_sup(interval));
        if nargin>=3 && ~isempty(K)
            fprintf('Planned bins (K): %d \n', K);
        else
            fprintf('Planned bins (K): inferred from interval length and omega_N\n');
        end
        fprintf('Output CSV: %s\n', self.outFile);
    end

    function [mat,u1,u2,u3] = prepareFEM(self)
        tri_I = [ self.callI_intval('0'), self.callI_intval('0'); ...
                  self.callI_intval('1'), self.callI_intval('0'); ...
                  1/self.callI_intval(2), sqrt(self.callI_intval(3))/2 ];
        [lamhs, inner_uhat, Agrad, AL2, Axx, Axy, Ayy, ~] = ...
            calc_eigen_bounds_any_order_for_quotients(tri_I, self.mesh_N, self.ord); %#ok<ASGLU>
        mat = struct('xx',Axx,'xy',Axy,'yy',Ayy,'grad',Agrad,'l2',AL2);

        u1 = inner_uhat(:,1);
        u2 = inner_uhat(:,2);
        u3 = inner_uhat(:,3);
        u1 = u1 / sqrt(self.L2(mat,u1,u1));
        u2 = u2 / sqrt(self.L2(mat,u2,u2));
        u3 = u3 - u2 * self.L2(mat,u3,u2);
        u3 = u3 / sqrt(self.L2(mat,u3,u3));
    end

    function tI = interval_t(self)
        epI = self.callI_intval(num2str(self.ep));
        tI  = self.callI_hull(0, epI);
    end

    function deltaI = interval_delta(self, d_pair)
        deltaI = self.callI_hull(d_pair(1), d_pair(2));
    end

    function [approx_mu, interval_mu] = quotient_calc(self, t, delta, mat, u1, u2, u3)
        d = self.d_a(delta); a = d(1); b = d(2);
        x = self.callI_intval('0.5'); tx = x + t*a;
        y = sqrt(self.callI_intval(3)) / 2; ty = y + t*b;
        Pt = [ t*a/(ty^2), -a*y/(ty^2); -a*y/(ty^2), -b*(y+ty)/(ty^2) ];

        S    = [1, (tx-x)/y; 0, ty/y];
        Sinv = [1, (x-tx)/ty; 0, y/ty];
        pu_  = self.callI_intval(self.callI_sup(norm(S*S',2)^2));
        pl   = self.callI_intval(self.callI_inf(1 / norm(S*S',2)^2));
        % pu := self.callI_intval(self.callI_sup(norm(Sinv*Sinv',2)^2)); %#ok<NASU,DEF> (kept for parity)

        lam = [ self.callI_intval('16')/self.callI_intval('3')*self.callI_intval('pi')^2, ...
                self.callI_intval('112')/self.callI_intval('9')*self.callI_intval('pi')^2, ...
                self.callI_intval('112')/self.callI_intval('9')*self.callI_intval('pi')^2, ...
                self.callI_intval('64')/self.callI_intval('3')*self.callI_intval('pi')^2 ];
        lam1=lam(1); lam2=lam(2); lam3=lam(3); lam4=lam(4);

        rho = lam4 * pl;

        lam1h = self.callI_intval(self.callI_sup(self.grad(mat,u1,u1)/self.L2(mat,u1,u1)));
        lam3h = self.callI_intval(self.callI_sup(self.grad(mat,u3,u3)/self.L2(mat,u3,u3)));
        lam1_tilde = (lam1) * pu_;
        lam3_tilde = (lam3) * pu_;

        matF=[self.L2(mat,u1,u2); self.L2(mat,u1,u3)];
        matG=self.L2(mat,u1,u1);
        matH=[ self.L2(mat,u2,u2), self.L2(mat,u2,u3); self.L2(mat,u3,u2), self.L2(mat,u3,u3)];
        etaF=norm(matF*matF',2); etaG=norm(eye(1)-matG,2); etaH=norm(eye(2)-matH,2);
        epb_hat12 = etaF / ((1-etaG)*(1-etaH));

        dbE1E1h = sqrt((lam1h - lam1) / (lam4 - lam1h));
        thb2    = (lam4 - lam2) * (epb_hat12 + dbE1E1h)^2;
        dbE2E2h = sqrt((lam3h - lam2 + thb2) / (lam4 - lam2));
        dbbarE2E2h = sqrt(2 - 2*sqrt(1 - dbE2E2h^2));

        dbE1E1t = sqrt((lam1_tilde - lam1) / (rho - lam1));
        thb2t   = (rho - lam2) * dbE1E1t^2;
        dbE2E2t = sqrt((lam3_tilde - lam2 + thb2t) / (rho - lam2));
        dbbarE2E2t = sqrt(2 - 2*sqrt(1 - dbE2E2t^2));

        eta_hat   = sqrt(lam3 * dbbarE2E2h^2 + lam3h - lam2);
        eta_tilde = sqrt(lam3 * dbbarE2E2t^2 + lam3_tilde - lam2);

        Err_F = self.callI_intval(self.callI_sup(norm(Pt,2) * (2*sqrt(lam3h)*eta_hat + sqrt(lam3)*eta_tilde)));
        Err_b = self.callI_intval(self.callI_sup(2*dbbarE2E2h + dbbarE2E2t));

        M=[ self.F(mat,Pt,u3,u3), self.F(mat,Pt,u3,u2); self.F(mat,Pt,u2,u3), self.F(mat,Pt,u2,u2)];
        N=[ self.L2(mat,u2,u2),   self.L2(mat,u2,u3);   self.L2(mat,u3,u2),   self.L2(mat,u3,u3)];

        M_ = M + self.callI_hull(-Err_F, Err_F);
        N_ = N + self.callI_hull(-Err_b, Err_b);

        [V,D]   = eig(self.callI_mid(M),  self.callI_mid(N));
        [mu2,~] = verifyeig(M,  D(1,1), V(:,1), N);
        [mu3,~] = verifyeig(M,  D(2,2), V(:,2), N);
        approx_mu = [mu2, mu3];
        [~,ind]   = sort(self.callI_mid(approx_mu));
        approx_mu = approx_mu(ind);

        [V,D]     = eig(self.callI_mid(M_), self.callI_mid(N_));
        [mu2_,~]  = verifyeig(M_, D(1,1), V(:,1), N_);
        [mu3_,~]  = verifyeig(M_, D(2,2), V(:,2), N_);
        interval_mu = [mu2_, mu3_];
        [~,ind]     = sort(self.callI_mid(interval_mu));
        interval_mu = interval_mu(ind);
    end

    %---------- Interval wrappers ----------%
    function y = callI_intval(~, x), y = I_intval(x); end
    function y = callI_hull(~, a,b), y = I_hull(a,b); end
    function y = callI_mid(~, x),     y = I_mid(x);   end
    function y = callI_inf(~, x),     y = I_inf(x);   end
    function y = callI_sup(~, x),     y = I_sup(x);   end

    %---------- Bilinear & direction ----------%
    function f = L2(~,mat,u,v),   f = u' * mat.l2   * v; end
    function f = grad(~,mat,u,v), f = u' * mat.grad * v; end
    function f = F(~,mat,Pt,u,v)
        uxvx = u' * mat.xx * v; uyvy = u' * mat.yy * v;
        uxvy = u' * mat.xy * v; uyvx = v' * mat.xy * u;
        f = Pt(1,1)*uxvx + Pt(1,2)*uyvx + Pt(2,1)*uxvy + Pt(2,2)*uyvy;
    end
    function d = d_a(~, delta), d = [sin(delta), -cos(delta)]; end

    %---------- CSV utils ----------%
    function ensureCsvHeader(~, file)
        if ~isfile(file)
            fid = fopen(file,'w'); if fid<0, error('Failed to create %s',file); end
            fprintf(fid,'idx,inf_mu1,sup_mu1,inf_mu2,sup_mu2\n');
            fclose(fid);
        end
    end

    function have = readResultsIndices(~, file)
        have = [];
        try
            Tt = readtable(file);
            if any(strcmp('idx', Tt.Properties.VariableNames))
                have = unique(Tt.idx(~isnan(Tt.idx))); have = have(:).'; return;
            end
        catch, end
        try
            M = readmatrix(file);
            if any(isnan(M(1,:))), M = M(2:end,:); end
            if ~isempty(M) && size(M,2)>=1
                have = unique(M(:,1)); have = have(~isnan(have)); have = have(:).';
            end
        catch, end
    end

    %---------- Helpers ----------%
    function K = inferBinCount(~, interval, omegaN)
        fullLen = self.callI_intval('pi')/3;
        frac = max(0, min(1, (interval(2)-interval(1))/fullLen));
        K = max(1, round(omegaN * frac));
    end

    function s = progressBar(~, k,N,width)
        if nargin<4, width = 32; end
        filled = round(width * k / max(N,1));
        s = ['[', repmat('=',1,filled), repmat('.',1,width-filled), ']'];
    end

    function s = formatDuration(~, sec)
        if ~isfinite(sec) || sec<0, s = 'n/a'; return; end
        mins = floor(sec/60); hrs = floor(mins/60);
        s = sprintf('%02dh:%02dm:%02ds', hrs, rem(mins,60), round(rem(sec,60)));
    end

    function assertFile(~, pth, msg)
        if ~isfile(pth), error('%s (missing: %s)', msg, pth); end
    end

    function [n0,n2,n1,tot] = readAlgo2States(~, csvPath)
        % States: "0"=pending, "2"=running (locked), "1"=done
        M = readmatrix(csvPath, 'OutputType','double');
        if isempty(M), n0=0;n2=0;n1=0;tot=0; return; end
        % Expected format: j,state per row (2 columns)
        if size(M,2) < 2
            % If commas absent or header present, try to coerce
            M = M(~any(isnan(M),2),:);
            if size(M,2) < 2, n0=0;n2=0;n1=0;tot=0; return; end
        end
        st = M(:,2);
        n0 = sum(st==0); n2 = sum(st==2); n1 = sum(st==1);
        tot = numel(st);
    end
end

%==================== Static defaults ====================%
methods (Static)
    function lams = default_boundsFcn(tri_intval,N_LG,N_rho,ord,isLG)
        if exist('calc_eigen_bounds_any_order','file')==2
            lams = calc_eigen_bounds_any_order(tri_intval,N_LG,N_rho,ord,isLG);            
        else
            error(['default_boundsFcn: calc_eigen_bounds_any_order not found.\n' ...
                   'Provide a custom function via the ''boundsFcn'' property.']);
        end
    end
end

end
