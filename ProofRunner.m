classdef ProofRunner < handle
% A runner for Algorithm 1 & 2.
%
% This class manages the execution flow for proving the spectral gap.
% Algorithm 1: verification_step_1.m (Local/Perturbation near equilateral)
% Algorithm 2: verification_step_2.m (Global/Domain Monotonicity)
%
% QUICK START
%   s = ProofRunner;
%   s.setupAll();                      % INTLAB setup
%   s.runAlgo1All();                   % Algo 1: full Omega_up sweep
%   s.runAlgo1Interval([intval('0.1'), intval('0.2')], 5); 
%   s.runAlgo2All();                      % Algo 2: full Omega_down sweep
%   s.summarizeAlgo1CSV();             % check sup(mu1) < inf(mu2)
%   s.summarizeAlgo2CSV();         % check sup(lam2) < inf(lam3)
%
% Author: Ryoki Endo and Xuefeng Liu

properties
    % ===== Algorithm 1 global parameters (Omega_up) =====
    omega_N (1,1) double = 1000      % bins over [0, pi/3]
    mesh_N  (1,1) double = 32        % FEM mesh parameter
    mesh_LG (1,1) double = 8         % FEM mesh parameter for L-G method
    ord     (1,1) double = 5         % FEM polynomial order
    ep      (1,:) char = '4e-5'      % t in [0, ep]

    % Output CSV for Algorithm 1
    algo1_outFile (1,:) char = 'results/results_algo1.csv'

    % ===== Algorithm 2 parameters (Omega_down) =====
    algo2_InputFile (1,:) char = 'inputs/cell_def.csv'  % Base path for Algo 2 inputs
    algo2_OutFile (1,:) char = 'results/results_algo2.csv'    % Base path for Algo 2 results

    % Behavior
    verbose (1,1) logical = true
    resume  (1,1) logical = true

    % Bounds callback for point/box evidence
    boundsFcn function_handle = @ProofRunner.default_boundsFcn
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
        
        % Ensure results directory exists
        if ~exist('results','dir'), mkdir('results'); end
        
        % Default filenames
        ts = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        if isempty(self.algo1_outFile)
            self.algo1_outFile = fullfile('results',sprintf('quotients_%s.csv',ts));
        end
        if isempty(self.algo2_OutFile)
            self.algo2_OutFile = fullfile('results', 'results_algo2.csv');
        end
    end

    %==================== Setup (INTLAB) ====================%
    function setupAll(self)
        self.setupIntlab();
    end

    function setupIntlab(self)
        try
            my_intlab_config();
        catch ME
            warning('INTLAB setup call failed: %s', ME.message);
        end
        format compact short infsup
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
        %    Corresponds to t = I_hull(0, ep)
        t_interval = I_hull(0, I_intval(self.ep));

        % 3. Partition the angular interval
        edges = linspace(interval(1), interval(2), K+1);
        idxList = 1:K;

        % 4. Resume capability
        doneMask = false(size(idxList));
        if self.resume && isfile(self.algo1_outFile)
            have = self.readResultsIndices(self.algo1_outFile);
            doneMask(ismember(idxList, have)) = true;
            if self.verbose
                fprintf('Resume: %d/%d bins already present -> skipped.\n', sum(doneMask), K);
            end
        end
        todo = idxList(~doneMask);
        total = numel(todo);

        if total == 0
            if self.verbose, fprintf('Nothing to do for Algorithm 1.\n'); end
            return;
        end

        self.ensureCsvHeader(self.algo1_outFile);

        % Warm-up for ETA
        if self.verbose && total >= 1
            d0 = I_hull(edges(todo(1)), edges(todo(1)+1));
            t0 = tic;
            try
                % Call the external function verification_step_1
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
                
                writematrix(row, self.algo1_outFile, 'WriteMode','append');
                
                % detailed print per step
                gap = I_inf(interval_mu(2) - interval_mu(1));
                fprintf('Bin %d: delta=[%.5g,%.5g] gap >= %.5g\n', ...
                   i, I_inf(edges(i)), I_sup(edges(i+1)), gap);
                
            catch ME
                warning('Bin %d failed: %s', i, ME.message);
                writematrix([i, NaN, NaN, NaN, NaN], self.algo1_outFile, 'WriteMode','append');
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
            fprintf('Algorithm 1 Done. Results -> %s\n', self.algo1_outFile);
        end
    end

    function summarizeAlgo1CSV(self, file)
        if nargin<2 || isempty(file), file = self.algo1_outFile; end
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
    function runAlgo2All(self, varargin)
        % RUNALGO2 Executes Algorithm 2 with progress bar and ETA.
        
        self.setupIntlab();

        verification_step_2(self.algo2_InputFile, self.algo2_OutFile);
        
        if self.verbose
            fprintf('Algorithm 2 batch completed.\n');
        end
    end

    function summarizeAlgo2CSV(self, file)
        % summarizeAlgo2CSV reads the Algorithm 2 result CSV 
        % and checks the verified inclusion property.
        %
        % CSV Format: i, x_inf, x_sup, theta_inf, theta_sup, lam2_sup, lam3_inf
        
        % Use default file if none provided
        if nargin < 2 || isempty(file), file = self.algo2_OutFile; end
        
        if ~isfile(file)
            fprintf('Algorithm 2 results file not found: %s\n', file);
            return;
        end

        fprintf('--- Verifying Algorithm 2 Results: %s ---\n', file);
        
        try
            % Read matrix
            data = readmatrix(file);
            
            % Skip header row if it exists and results in NaNs
            if size(data,1) > 0 && any(isnan(data(1,:)))
                data = data(2:end,:);
            end
            
            if isempty(data)
                fprintf('File is empty or could not be read.\n');
                return; 
            end
            
            % --- Update Column Indices ---
            % New Format: 
            % Col 1: i
            % Col 6: sup(lambda_2)
            % Col 7: inf(lambda_3)
            sup_lam2 = data(:, 6);
            inf_lam3 = data(:, 7);
            
            % Check separation: sup(lam2) < inf(lam3)
            fails = sup_lam2 >= inf_lam3;
            n_fail = sum(fails);
            total_checked = size(data, 1);
            
            if n_fail > 0
                fprintf('FAIL: %d violations found.\n', n_fail);
                % Print the indices (i) of the failed rows
                failed_indices = data(fails, 1);
                fprintf('Failed Indices: %s\n', mat2str(failed_indices', 5));
            else
                fprintf('>>> PROOF SUCCESSFUL (Algo 2): sup(lam2) < inf(lam3) verified for all %d rows.\n', total_checked);
            end
                
        catch ME
            warning('Could not read or parse %s: %s', file, ME.message);
        end
    end

    %==================== Partial Evidence / Utils ====================%
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
    
    function [up2, lo3] = boundsOnBox(self, x_lo, x_hi, theta_lo, theta_hi, varargin)
        % BOUNDSONBOX Computes rigorous eigenvalue bounds over a box in (x, theta) space.
        %
        % Parameters:
        %   x_lo, x_hi         : Interval for x-coordinate
        %   theta_lo, theta_hi : Interval for polar angle theta (radians)
        %
        % Logic based on Algorithm 4:
        %   The domain inclusion T^{x_lo, theta_lo} \subset T^p \subset T^{x_hi, theta_hi} holds.
        %   - sup(lambda_2) is obtained from the "inner" (smallest) geometry: (x_lo, theta_lo).
        %   - inf(lambda_3) is obtained from the "outer" (largest) geometry: (x_hi, theta_hi).
        
        p = inputParser;
        addParameter(p,'N_LG', self.mesh_LG);
        addParameter(p,'N_rho', self.mesh_N);
        addParameter(p,'ord', self.ord);
        parse(p, varargin{:}); q = p.Results;
        
        % 1. Construct the "Smallest" Triangle (Inner)
        %    Corresponds to p_{i,j} in the paper -> Yields Upper Bounds
        y_small = x_lo * tan(theta_lo);
        tri_small = [I_intval(0), I_intval(0); I_intval(1), I_intval(0); x_lo, y_small];

        % 2. Construct the "Largest" Triangle (Outer)
        %    Corresponds to p_{i+1,j+1} in the paper -> Yields Lower Bounds
        y_large = x_hi * tan(theta_hi);
        tri_large = [I_intval(0), I_intval(0); I_intval(1), I_intval(0); x_hi, y_large];
        
        % Compute eigenvalues
        % Upper bound of lambda_2 comes from the smallest domain
        lams_sup = self.boundsFcn(tri_small, q.N_LG, q.N_rho, q.ord, 1);
        % Lower bound of lambda_3 comes from the largest domain
        lams_inf = self.boundsFcn(tri_large, q.N_LG, q.N_rho, q.ord, 1);
        
        up2 = I_sup(lams_sup(1)); 
        lo3 = I_inf(lams_inf(2));
        
        if self.verbose
            fprintf('Box Bounds [(x,theta)]: sup(lam2) <= %.17g, inf(lam3) >= %.17g\n', up2, lo3);
        end
    end

    function explain(self)
        fprintf('ProofRunner manages the verification of eigenvalue gaps.\n');
        fprintf('Algorithm 1: Uses verification_step_1.m to verify near-equilateral regions.\n');
        fprintf('Algorithm 2: Uses verification_step_2.m to verify global regions via monotonicity.\n');
    end
end

%==================== Private Methods ====================%
methods (Access = private)
    function printAlgo1Preface(self, interval, K)
        fprintf('--- Algorithm 1: Rigorous difference-quotients ---\n');
        fprintf('Parameters: mesh_N=%d, ord=%d, ep=%s, omega_N=%d\n', ...
            self.mesh_N, self.ord, self.ep, self.omega_N);
        fprintf('Bins: %d\nOutput: %s\n', K, self.algo1_outFile);
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
end

methods (Static)
    function lams = default_boundsFcn(tri, NL, NR, ord, lg)
        lams = calc_eigen_bounds_any_order(tri, NL, NR, ord, lg);
    end
end

end