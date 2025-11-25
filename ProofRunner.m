classdef ProofRunner < handle
% A runner for Algorithm 1 & 2.
%
% This class manages the execution flow but delegates the core rigorous
% computation of Algorithm 1 to the external function 'verification_step_1.m'.
%
% QUICK START
%   s = ProofRunner;
%   s.setupAll();                      % INTLAB + prep.sh
%   s.runAlgo1All();                   % full Omega_up sweep
%   s.summarizeAlgo1CSV();             % check sup(mu1) < inf(mu2)
%   s.summarizeAlgo2Results();         % check sup(lam2) < inf(lam3) for Algo 2
%
% Author: Ryoki Endo and Xuefeng Liu

properties
    % ===== Algorithm 1 global parameters =====
    omega_N (1,1) double = 1000      % bins over [0, pi/3]
    mesh_N  (1,1) double = 32        % FEM mesh parameter
    mesh_LG (1,1) double = 8         % FEM mesh parameter for L-G method
    ord     (1,1) double = 5         % FEM polynomial order
    ep      (1,:) char = '1e-5'      % t in [0, ep]

    % Output CSV (if empty -> auto timestamp in results/)
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
        % Name-value pairs constructor
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
        self.setupIntlab();
        self.prepareAlgo2();
    end

    function setupIntlab(self)
        try
            if exist('my_intlab_config_alone','file')==2
                my_intlab_config_alone;
            elseif exist('my_intlab_config','file')==2
                my_intlab_config;
            else
                warning('No my_intlab_config found.');
            end
        catch ME
            warning('INTLAB setup call failed: %s', ME.message);
        end
        format compact short infsup
    end

    function prepareAlgo2(self)
        if ~isunix
            error('Algorithm 2 preparation requires a Unix-like shell.');
        end
        self.assertFile(self.prepScript,'prep.sh not found.');
        cmd = sprintf('bash %s', self.prepScript);
        if self.verbose, fprintf('Running: %s\n', cmd); end
        st = system(cmd);
        if st ~= 0
            error('prep.sh failed (exit code %d).', st);
        end
    end

    %==================== Algorithm 1 (Omega_up) ====================%
    function runAlgo1All(self)
        % Run the full sweep from 0 to pi/3
        full_interval = [I_intval(0), I_intval('pi')/3];
        self.runAlgo1Interval(full_interval, self.omega_N);
    end

    function runAlgo1Direction(self, delta)
        % Run for a single direction
        self.runAlgo1Interval([delta, delta], 1);
    end

    function runAlgo1Interval(self, interval, K)
        % Main loop for Algorithm 1
        self.setupIntlab();
        
        if nargin<3 || isempty(K)
            K = self.inferBinCount(interval, self.omega_N);
        end

        self.printAlgo1Preface(interval, K);

        % 1. Prepare FEM matrices and Basis (Pre-computation)
        [mat, u1, u2, u3] = self.prepareFEM();

        % 2. Prepare the perturbation interval 't'
        %    Corresponds to t = I_hull(0, ep) in the script.
        t_interval = I_hull(0, I_intval(self.ep));

        % 3. Partition the angular interval
        edges = linspace(interval(1), interval(2), K+1);
        idxList = 1:K;

        % 4. Resume capability
        doneMask = false(size(idxList));
        if self.resume && isfile(self.outFile)
            have = self.readResultsIndices(self.outFile);
            doneMask(ismember(idxList, have)) = true;
            if self.verbose
                fprintf('Resume: %d/%d bins already present -> skipped.\n', sum(doneMask), K);
            end
        end
        todo = idxList(~doneMask);
        total = numel(todo);

        if total == 0
            if self.verbose, fprintf('Nothing to do.\n'); end
            return;
        end

        self.ensureCsvHeader(self.outFile);

        % Warm-up for ETA
        if self.verbose && total >= 1
            d0 = I_hull(edges(todo(1)), edges(todo(1)+1));
            t0 = tic;
            try
                % Call the external function
                verification_step_1(t_interval, d0, mat, u1, u2, u3);
                step_t = toc(t0);
                fprintf('Planned steps: %d | warm-up per-step ~ %.2fs | ETA ~ %s\n', ...
                    total, step_t, self.formatDuration(step_t*total));
            catch ME
                warning('Warm-up failed (%s).', ME.message);
            end
        end

        % 5. Main Execution Loop
        t_start = tic;
        if self.verbose
            fprintf('--- Running Algorithm 1 on %d bin(s) ---\n', total);
        end

        for kk = 1:total
            i = todo(kk);
            
            % Define the interval 'delta' for the current sector
            delta = I_hull(edges(i), edges(i+1));
            
            try
                % --- CALL EXTERNAL FUNCTION ---
                interval_mu = verification_step_1(t_interval, delta, mat, u1, u2, u3);
                
                % Store results: index and rigorous bounds
                row = [i, I_inf(interval_mu(1)), I_sup(interval_mu(1)), ...
                          I_inf(interval_mu(2)), I_sup(interval_mu(2))];
                
                writematrix(row, self.outFile, 'WriteMode','append');
                
                % detailed print per step
                fprintf('delta=[%.17g,%.17g]: gap >= %.17g\n', ...
                   I_inf(edges(i)), I_sup(edges(i+1)), I_inf(interval_mu(2) - interval_mu(1)));
                
            catch ME
                warning('Bin %d failed: %s', i, ME.message);
                writematrix([i, NaN, NaN, NaN, NaN], self.outFile, 'WriteMode','append');
            end

            if self.verbose
                pct = 100*kk/total;
                elapsed = toc(t_start);
                eta_sec = (elapsed/max(kk,1))*(total-kk);
                fprintf('\r%s %d/%d (%.1f%%) ETA %s', ...
                        self.progressBar(kk,total,32), kk, total, pct, self.formatDuration(eta_sec));
                if kk==total, fprintf('\n'); end
            end
        end
        
        if self.verbose
            fprintf('Done. Results -> %s\n', self.outFile);
        end
    end

    function summarizeAlgo1CSV(self, file)
        if nargin<2 || isempty(file), file = self.outFile; end
        if ~isfile(file), error('CSV not found: %s', file); end

        T = readmatrix(file);
        % Skip header if it exists and causes NaN in first row
        if size(T,1) > 0 && any(isnan(T(1,:))), T = T(2:end,:); end
        
        if size(T,2) < 5
            error('Unexpected CSV format. Expected at least 5 cols.');
        end

        % Column order: idx, inf_mu1, sup_mu1, inf_mu2, sup_mu2
        sup_mu1 = T(:,3);
        inf_mu2 = T(:,4);
        
        ok = sup_mu1 < inf_mu2;

        fprintf('File: %s\n', file);
        fprintf('Bins: %d | OK: %d | NG: %d\n', numel(ok), sum(ok), sum(~ok));
        if all(ok)
            fprintf('>>> PROOF SUCCESSFUL: sup(mu1) < inf(mu2) for all checked bins.\n');
        else
            fprintf('>>> PROOF FAILED in bins: %s\n', num2str(find(~ok)'));
        end
    end

    %==================== Algorithm 2 (Omega_down) ====================%
    function runAlgo2(self, varargin)
        p = inputParser;
        addParameter(p,'detached',false,@islogical);
        addParameter(p,'logFile','Each_Process/log/main_algo2_launcher.log');
        parse(p, varargin{:});
        opt = p.Results;

        if ~isunix
            error('Algorithm 2 requires Unix shell.');
        end
        self.assertFile(self.mainScript2, 'main_algo2.sh not found.');
        
        if opt.detached
            logFile = char(opt.logFile);
            [logDir,~] = fileparts(logFile);
            if ~isempty(logDir) && ~exist(logDir,'dir'), mkdir(logDir); end
            cmd = sprintf('nohup bash %s > %s 2>&1 &', self.mainScript2, logFile);
            if self.verbose, fprintf('Launching (detached): %s\n', cmd); end
            st = system(cmd);
            if st ~= 0, error('Launch failed (exit %d).', st); end
        else
            cmd = sprintf('bash %s', self.mainScript2);
            if self.verbose, fprintf('Running: %s\n', cmd); end
            st = system(cmd);
            if st ~= 0, error('Failed (exit %d).', st); end
        end
    end

    function monitorAlgo2(self, varargin)
        p = inputParser;
        addParameter(p,'csv', self.listCsvPath);
        addParameter(p,'pollSeconds', self.pollSeconds);
        parse(p, varargin{:});
        
        csvPath = char(p.Results.csv);
        poll    = p.Results.pollSeconds;
        self.assertFile(csvPath, 'CSV list not found.');

        fprintf('--- Algorithm 2 monitor (%s) ---\n', csvPath);
        donePrev = 0; tLast = tic; t0 = tic;

        while true
            [n0,n2,n1,tot] = self.readAlgo2States(csvPath);
            done = n1;
            
            dt = toc(tLast); tLast = tic;
            rate = max(done - donePrev, 0) / max(dt, eps);
            rem  = max(tot - done, 0);
            eta  = rem / max(rate, eps);
            donePrev = done;

            fprintf('\rTot=%d | Pend=%d | Run=%d | Done=%d (%.1f%%) | ETA %s', ...
                tot, n0, n2, n1, 100*done/max(tot,1), self.formatDuration(eta));

            if done >= tot
                fprintf('\nAll tasks finished.\n');
                break;
            end
            pause(poll);
        end
    end

    function summarizeAlgo2Results(self, resultsDir)
        % summarizeAlgo2Results reads all CSV files from Algorithm 2 
        % and checks the verified inclusion property.
        
        if nargin < 2, resultsDir = 'results'; end
        
        % Look for files matching the pattern output by Algorithm2.m
        % Pattern: step2_bounds_<min_j>_<max_j>.csv
        files = dir(fullfile(resultsDir, 'step2_bounds_*.csv'));
        
        if isempty(files)
            fprintf('No Algorithm 2 results found in %s\n', resultsDir);
            return;
        end

        total_checked = 0;
        failed_count = 0;
        
        fprintf('--- Verifying Algorithm 2 Results in %s ---\n', resultsDir);
        
        for k = 1:length(files)
            fname = fullfile(files(k).folder, files(k).name);
            try
                % Read matrix, skipping the header line
                % CSV Format: i, j, inf_lam2, sup_lam2, inf_lam3, sup_lam3
                data = readmatrix(fname, 'NumHeaderLines', 1);
                
                if isempty(data)
                    continue; 
                end
                
                % Extract relevant columns
                % col 4: sup(lambda_2), col 5: inf(lambda_3)
                sup_lam2 = data(:, 4);
                inf_lam3 = data(:, 5);
                
                % Check separation: sup(lam2) < inf(lam3)
                fails = sup_lam2 >= inf_lam3;
                n_fail = sum(fails);
                
                total_checked = total_checked + size(data, 1);
                failed_count = failed_count + n_fail;
                
                if n_fail > 0
                    fprintf('FAIL: %s (%d violations found)\n', files(k).name, n_fail);
                    % Optional: Print specific failing indices
                    % idx = find(fails);
                    % disp(data(idx, 1:2)); 
                end
                
            catch ME
                warning('Could not read or parse %s: %s', fname, ME.message);
            end
        end
        
        fprintf('--------------------------------------------------\n');
        fprintf('Total subregions checked: %d\n', total_checked);
        
        if total_checked > 0
            if failed_count == 0
                fprintf('>>> PROOF SUCCESSFUL (Algo 2): sup(lam2) < inf(lam3) verified for all loaded data.\n');
            else
                fprintf('>>> PROOF FAILED (Algo 2): %d violations found.\n', failed_count);
            end
        else
            fprintf('>>> WARNING: No valid data entries processed.\n');
        end
    end

    %==================== Partial Evidence ====================%
    function [lam2, lam3] = boundsAtPoint(self, a, b, varargin)
        p = inputParser;
        addParameter(p,'N_LG', self.mesh_LG);
        addParameter(p,'N_rho', self.mesh_N);
        addParameter(p,'ord', self.ord);
        parse(p, varargin{:}); q = p.Results;

        tri = [I_intval(0),I_intval(0); I_intval(1),I_intval(0); a, b];
        lams = self.boundsFcn(tri, q.N_LG, q.N_rho, q.ord, 1);
        lam2 = lams(1); lam3 = lams(2);
        
        if self.verbose
            fprintf('Point Bounds:\n lam2 in [%.17g, %.17g]\n lam3 in [%.17g, %.17g]\n', ...
                I_inf(lam2), I_sup(lam2), I_inf(lam3), I_sup(lam3));
        end
    end
    
    function [up2, lo3] = boundsOnBox(self, r_lo, r_hi, h_lo, h_hi, varargin)
        % BOUNDSONBOX Computes rigorous eigenvalue bounds over a box in (r,h) space.
        % The function converts (r,h) back to (x,y) to construct the triangle vertices.
        %
        % Inputs:
        %   r_lo, r_hi : Interval bounds for r = sqrt(x^2+y^2)
        %   h_lo, h_hi : Interval bounds for h = y

        p = inputParser;
        addParameter(p,'N_LG', self.mesh_LG);
        addParameter(p,'N_rho', self.mesh_N);
        addParameter(p,'ord', self.ord);
        parse(p, varargin{:}); q = p.Results;
        
        % 1. Construct the "Smallest" Triangle (for Upper Bound of lambda2)
        %    Smallest r and smallest h yield the smallest domain -> Largest Eigenvalues.
        %    Coordinate conversion: x = sqrt(r^2 - h^2), y = h
        x_small = sqrt(r_lo^2 - h_lo^2);
        y_small = h_lo;
        
        tri_small = [I_intval(0), I_intval(0); ...
                    I_intval(1), I_intval(0); ...
                    x_small,     y_small];

        % 2. Construct the "Largest" Triangle (for Lower Bound of lambda3)
        %    Largest r and largest h yield the largest domain -> Smallest Eigenvalues.
        x_large = sqrt(r_hi^2 - h_hi^2);
        y_large = h_hi;
        
        tri_large = [I_intval(0), I_intval(0); ...
                    I_intval(1), I_intval(0); ...
                    x_large,     y_large];
        
        % Compute eigenvalues
        % Calculate lambda2 on the smallest domain (Potential Maximum of lambda2)
        lams_sup = self.boundsFcn(tri_small, q.N_LG, q.N_rho, q.ord, 1);
        
        % Calculate lambda3 on the largest domain (Potential Minimum of lambda3)
        lams_inf = self.boundsFcn(tri_large, q.N_LG, q.N_rho, q.ord, 1);
        
        % Extract rigorous bounds
        % lams_sup(2) is the interval for lambda2 on the small domain. We take its sup.
        % lams_inf(3) is the interval for lambda3 on the large domain. We take its inf.
        up2 = I_sup(lams_sup(1)); 
        lo3 = I_inf(lams_inf(2));
        
        if self.verbose
            fprintf('Box Bounds [(r,h)]: sup(lam2) <= %.17g, inf(lam3) >= %.17g\n', up2, lo3);
        end
    end

    function explain(self)
        fprintf('ProofRunner manages the verification of eigenvalue gaps.\n');
        fprintf('Algorithm 1: Uses verification_step_1.m to verify near-equilateral regions.\n');
        fprintf('Algorithm 2: Uses prep.sh and main_algo2.sh for global regions.\n');
    end
end

%==================== Private Methods ====================%
methods (Access = private)
    function printAlgo1Preface(self, interval, K)
        fprintf('--- Algorithm 1: Rigorous difference-quotients ---\n');
        fprintf('Parameters: mesh_N=%d, ord=%d, ep=%s, omega_N=%d\n', ...
            self.mesh_N, self.ord, self.ep, self.omega_N);
        fprintf('Bins: %d\nOutput: %s\n', K, self.outFile);
    end

    function [mat,u1,u2,u3] = prepareFEM(self)
        % Pre-calculation on the exact equilateral triangle
        tri_I = [ I_intval('0'), I_intval('0'); ...
                  I_intval('1'), I_intval('0'); ...
                  1/I_intval('2'), sqrt(I_intval('3'))/2 ];
              
        [~, inner_uhat, Agrad, AL2, Axx, Axy, Ayy, ~] = ...
            calc_eigen_bounds_any_order_for_quotients(tri_I, self.mesh_LG, self.ord);
        
        mat = struct('xx',Axx,'xy',Axy,'yy',Ayy,'grad',Agrad,'l2',AL2);
        u1 = inner_uhat(:,1);
        u2 = inner_uhat(:,2);
        u3 = inner_uhat(:,3);
    end

    function ensureCsvHeader(~, file)
        if ~isfile(file)
            fid = fopen(file,'w');
            fprintf(fid,'idx,inf_mu1,sup_mu1,inf_mu2,sup_mu2\n');
            fclose(fid);
        end
    end

    function have = readResultsIndices(~, file)
        have = [];
        if ~isfile(file), return; end
        try
            M = readmatrix(file);
            % Handle possible header row causing NaNs
            if ~isempty(M) && any(isnan(M(1,:))), M = M(2:end,:); end
            if ~isempty(M), have = unique(M(:,1))'; end
        catch
            % If readmatrix fails (e.g. empty file), ignore
        end
    end

    function K = inferBinCount(~, interval, omegaN)
        fullLen = I_intval('pi')/3;
        % Using I_mid to get a standard double for calculation
        len = interval(2) - interval(1);
        frac = I_mid(len / fullLen);
        K = max(1, round(omegaN * frac));
    end

    function s = progressBar(~, k, N, width)
        pct = k / max(N,1);
        filled = round(width * pct);
        s = ['[', repmat('=',1,filled), repmat('.',1,width-filled), ']'];
    end

    function s = formatDuration(~, sec)
        if ~isfinite(sec) || sec<0, s = 'n/a'; return; end
        m = floor(sec/60); sec = sec - m*60;
        h = floor(m/60); m = m - h*60;
        s = sprintf('%02dh:%02dm:%02ds', h, m, round(sec));
    end

    function assertFile(~, pth, msg)
        if ~isfile(pth), error('%s (missing: %s)', msg, pth); end
    end

    function [n0,n2,n1,tot] = readAlgo2States(~, csvPath)
        n0=0; n2=0; n1=0; tot=0;
        try
            M = readmatrix(csvPath);
            if size(M,2) >= 2
                M = M(~any(isnan(M),2),:);
                st = M(:,2);
                n0=sum(st==0); n2=sum(st==2); n1=sum(st==1); tot=numel(st);
            end
        catch, end
    end
end

methods (Static)
    function lams = default_boundsFcn(tri, NL, NR, ord, lg)
        lams = calc_eigen_bounds_any_order(tri, NL, NR, ord, lg);
    end
end

end