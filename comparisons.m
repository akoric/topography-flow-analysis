clear; close all; clc;

%% 1. load dem results
S1     = load('W_std_results_v3.mat','T');  
T_U  = S1.T;
S2     = load('dem_results_v3.mat','T');       
Tdem   = S2.T;

U_1_mean = load('DEM1_Moments.mat','U_mean');    
U_2_mean = load('DEM2_Moments.mat','U_mean');    
U_3_mean = load('DEM3_Moments.mat','U_mean');    
U_4_mean = load('DEM4_Moments.mat','U_mean');    

% 2. parameters
paramNames = {'Ra_sw','Rq_sw','Ku_sw','Ly','Ly','Ly'}; 

% 3. find min frame for each topography --> obtain 'minFrames'
nParams    = numel(paramNames); 
demGroups  = unique(Tdem.DEM); 
nGroups    = numel(demGroups);
minFrames  = zeros(nGroups,1);
RowIdx  = zeros(nGroups,1);
frame_plus = 1;
for i = 1:nGroups
    g     = demGroups(i);
    idxs  = find(T_U.DEM == g);
    [~,loc] = min(T_U.Frame(idxs)); 
    loc = loc+frame_plus;
    minFrames(i) = T_U.Frame(idxs(loc));
    RowIdx(i) = idxs(loc);
end

val_U = T_U{RowIdx, paramNames};
val_dem = Tdem{[1 2 3 4], paramNames};
U_1 = U_1_mean.U_mean(1:599,:,minFrames(1));
U_1_max = max(U_1(:));

U_2 = U_2_mean.U_mean(1:599,:,minFrames(2));
U_2_max = max(U_2(:));

U_3 = U_3_mean.U_mean(1:599,:,minFrames(3));
U_3_max = max(U_3(:));

U_4 = U_4_mean.U_mean(1:599,:,minFrames(4));
U_4_max = max(U_4(:));

U_max = [U_1_max U_2_max U_3_max U_4_max]';

% normalization 
% 1. dem     [m]
% 2. U_mean  [m/s]     --> U_mean_max
% 3. U_var   [m2/s2]   --> U_mean_max_2
% 4. U_std   [m/s]     --> U_mean_max
% 5. U_skew  [nondim]
% 6. U_kurt  [nondim]
%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Ra' [same as input]
% 'Rq' [same as input]
% 'Sk' [nondim]
% 'Ku' [nondim]
% 'Lx' [m]
% 'Ly' [m]

Ra_U = val_U(:, 1)./U_max; % normalized by the maximum mean velocity
Rq_U = val_U(:, 2)./U_max; % normalized by the maximum mean velocity
% Ra_U = val_U(:, 1)./U_max.^2;  % for U_var
% Rq_U = val_U(:, 2)./U_max.^2; 
Sk_U = val_U(:, 3);
Ku_U = val_U(:, 4);
Lx_U = val_U(:, 5)/3; % normalized by the longer edge(600*0.005m) of topography domain
Ly_U = val_U(:, 6)/3; % normalized by the longer edge(600*0.005m) of topography domain

Ra_dem = val_dem(:, 1)/3; % normalized by the longer edge(600*0.005m) of topography domain
Rq_dem = val_dem(:, 2)/3; % normalized by the longer edge(600*0.005m) of topography domain
Sk_dem = val_dem(:, 3);
Ku_dem = val_dem(:, 4);
Lx_dem = val_dem(:, 5)/3; % normalized by the longer edge(600*0.005m) of topography domain
Ly_dem = val_dem(:, 6)/3; % normalized by the longer edge(600*0.005m) of topography domain

val_U_non = [Ra_U Rq_U Sk_U Ku_U Lx_U Ly_U];
val_dem_non = [Ra_dem Rq_dem Sk_dem Ku_dem Lx_dem Ly_dem];

%%
param_idx = 1:4;  % represents 'Ra','Rq','Sk','Ku','Lx','Ly'
power_range = -3:3;
max_num_params = 3;  % maximum number of parameters in each combination
R2_threshold = 0.98; % chosen threshold

% record all combinations exceeding the threshold
saved_results = [];

for k = 1:max_num_params
    param_combos = nchoosek(param_idx, k);

    for combo_idx = 1:size(param_combos, 1)
    idxs = param_combos(combo_idx, :);

    % generate all power exponent combinations (e.g. 5^k combinations)
        grid_cell = cell(1, k);
        [grid_cell{:}] = ndgrid(power_range);
        for j = 1:k
            grid_cell{j} = grid_cell{j}(:);
        end
        power_combos = [grid_cell{:}];

        for pi = 1:size(power_combos, 1)
            exponents = power_combos(pi, :);
            y1 = ones(size(val_U_non, 1), 1);
            y2 = ones(size(val_dem_non, 1), 1);

            try
                for j = 1:k
                    y1 = y1 .* val_U_non(:, idxs(j)).^exponents(j);
                    y2 = y2 .* val_dem_non(:, idxs(j)).^exponents(j);
                end

                % filter invalid values
                valid = isfinite(y1) & isfinite(y2);
                y1 = y1(valid); y2 = y2(valid);

                if numel(y1) < 3
                    continue;
                end

                % linear fit and R²
                p = polyfit(y1, y2, 1);
                y2_fit = polyval(p, y1);
                R2 = 1 - sum((y2 - y2_fit).^2) / sum((y2 - mean(y2)).^2);
                % save if above threshold
                if R2 > R2_threshold
                    result.param_ids = idxs;
                    result.powers = exponents;
                    result.R2 = R2;
                    result.coeff = p;
                    saved_results = [saved_results; result];
                end
            catch
                continue;
            end
        end
    end
end

% display total number of combinations
disp(['Found ', num2str(length(saved_results)), ' good combinations above R^2 = ', num2str(R2_threshold)]);

% optionally display the top combinations
[~, sort_idx] = sort([saved_results.R2], 'descend');
for i = 1:min(10, numel(sort_idx))
    r = saved_results(sort_idx(i));
    disp(['Params: ', mat2str(r.param_ids), ', Powers: ', mat2str(r.powers), ', R^2 = ', num2str(r.R2)]);
end


%%
param_idx = 1:6;  % represents 'Ra','Rq','Sk','Ku','Lx','Ly'
power_range = -3:3;
max_num_params = 2;  % tunable: maximum number of parameters used simultaneously

best_R2 = -Inf;
best_combo = [];

% iterate over parameter combination sizes
for k = 1:max_num_params
    param_combos = nchoosek(param_idx, k);

    for combo_idx = 1:size(param_combos, 1)
    idxs = param_combos(combo_idx, :);

    % create Cartesian product of power indices
        ranges = cell(1, k);
        [ranges{:}] = ndgrid(power_range);
        for j = 1:k
            ranges{j} = ranges{j}(:);
        end
        power_combos = [ranges{:}];


        for pi = 1:size(power_combos, 1)
            try
                exponents = power_combos(pi, :);
                y1 = ones(size(val_U_non, 1), 1);
                y2 = ones(size(val_dem_non, 1), 1);

                for j = 1:k
                    y1 = y1 .* val_U_non(:, idxs(j)).^exponents(j);
                    y2 = y2 .* val_dem_non(:, idxs(j)).^exponents(j);
                end

                % filter invalid values
                valid = isfinite(y1) & isfinite(y2);
                y1 = y1(valid); y2 = y2(valid);

                % linear fit and R²
                p = polyfit(y1, y2, 1);
                y2_fit = polyval(p, y1);
                R2 = 1 - sum((y2 - y2_fit).^2) / sum((y2 - mean(y2)).^2);

                if R2 > best_R2
                    best_R2 = R2;
                    best_combo = zeros(1,6);
                    best_combo(idxs) = exponents;
                    best_y1 = y1;
                    best_y2 = y2;
                    best_params = idxs;
                end
            catch
                continue;
            end
        end
    end
end

% output result
disp(['Best exponents: [a b c d e f] = ', mat2str(best_combo)]);
disp(['Best R^2 = ', num2str(best_R2)]);
disp(['Used params: ', mat2str(best_params)]);

% plotting
figure;
plot(best_y1, best_y2, 'ko', 'LineWidth', 2);
xlabel('U combination');
ylabel('DEM combination');
title(['Best linearity: R2 = ', num2str(best_R2)]);


%% multiple frames 
% 'Ra','Rq','Sk','Ku','Lx','Ly'
power_range     = -2:2;
demGroups       = unique(Tdem.DEM);
nGroups         = numel(demGroups);
framePlusRange  = 0:1;
nFP             = numel(framePlusRange);

valid_combos   = [];
worstR2_for_cb = [];
val_dem = Tdem{[1 2 3 4], paramNames};
Ra_dem = val_dem(:, 1)/3; % normalized by the longer edge(600*0.005m) of topography domain
Rq_dem = val_dem(:, 2)/3; % normalized by the longer edge(600*0.005m) of topography domain
Sk_dem = val_dem(:, 3);
Ku_dem = val_dem(:, 4);
Lx_dem = val_dem(:, 5)/3; % normalized by the longer edge(600*0.005m) of topography domain
Ly_dem = val_dem(:, 6)/3; % normalized by the longer edge(600*0.005m) of topography domain

val_dem_non = [Ra_dem Rq_dem Sk_dem Ku_dem Lx_dem Ly_dem];
for a = power_range
for b = power_range
for c = 0
for d = 0
for e = 0
for f = power_range

    R2s = nan(nFP,1);

    % compute R2 for each frame_plus
    for k = 1:nFP
    frame_plus = framePlusRange(k);

    % step 1: recompute minFrames and RowIdx
        minFrames = zeros(nGroups,1);
        RowIdx    = zeros(nGroups,1);
        for i = 1:nGroups
            idxs = find(T_U.DEM == demGroups(i));
            [~, loc0]     = min(T_U.Frame(idxs));
            loc           = min(loc0 + frame_plus, numel(idxs));
            minFrames(i)  = T_U.Frame(idxs(loc));
            RowIdx(i)     = idxs(loc);
        end
    % step 2: extract val_U_non for this frame_plus
        val_U = T_U{RowIdx, paramNames};
        Umax  = zeros(nGroups,1);
        for i = 1:nGroups
            Utmp       = eval(sprintf('U_%d_mean.U_mean(1:599,:,minFrames(%d))', i, i));
            Umax(i)    = max(Utmp(:));
        end

        Ra_U = val_U(:,1)./Umax;
        % Rq_U = val_U(:,2)./Umax.^2; % for U_var 
        % Ra_U = val_U(:,1)./Umax.^2;
        Rq_U = val_U(:,2)./Umax;
        Sk_U = val_U(:,3);
        Ku_U = val_U(:,4);
        Lx_U = val_U(:,5)/3;
        Ly_U = val_U(:,6)/3;
    val_U_non = [Ra_U Rq_U Sk_U Ku_U Lx_U Ly_U];

    % step 3: compute y1, y2 and R2 for this combination
        y1 = val_U_non(:,1).^a .* val_U_non(:,2).^b .* val_U_non(:,3).^c .* ...
             val_U_non(:,4).^d .* val_U_non(:,5).^e .* val_U_non(:,6).^f;
        y2 = val_dem_non(:,1).^a .* val_dem_non(:,2).^b .* val_dem_non(:,3).^c .* ...
             val_dem_non(:,4).^d .* val_dem_non(:,5).^e .* val_dem_non(:,6).^f;

        valid      = isfinite(y1) & isfinite(y2);
        x          = y1(valid);
        Y          = y2(valid);
        if numel(x) < 2
            R2s(k) = 0;
        else
            [p,S,mu] = polyfit(x, Y, 1);
            Yfit     = polyval(p, x, S, mu);
            R2s(k)   = 1 - sum((Y - Yfit).^2) / sum((Y - mean(Y)).^2);
        end
    end

    % if all R2 are above threshold, treat as valid combination
    if all(R2s >= 0.98)
        valid_combos   = [valid_combos; a b c d e f];
        worstR2_for_cb = [worstR2_for_cb; min(R2s)];
    end

end
end
end
end
end
end

% output result
if isempty(valid_combos)
    disp('No parameter combination meets the R^2 threshold for all frame_plus values.');
else
    T = array2table([valid_combos worstR2_for_cb], ...
        'VariableNames',{'a','b','c','d','e','f','minR2'});
    disp(T);

    [~, idxBest] = max(worstR2_for_cb);
    fprintf('Best fixed combination [a b c d e f] = %s, worst R^2 = %.5f\n',...
        mat2str(valid_combos(idxBest,:)), worstR2_for_cb(idxBest));
end


%% plot
a = 0;
b = -3;
c = 1;
d = -2;
e = 0;
f = 0;
%'Ra','Rq','Sk','Ku','Lx','Ly'
y1 = val_U_non(:, 1).^(a).*val_U_non(:, 2).^(b).*val_U_non(:, 3).^(c).*val_U_non(:, 4).^(d).*val_U_non(:, 5).^(e).*val_U_non(:, 6).^(f);
y2 = val_dem_non(:, 1).^(a).*val_dem_non(:, 2).^(b).*val_dem_non(:, 3).^(c).*val_dem_non(:, 4).^(d).*val_dem_non(:, 5).^(e).*val_dem_non(:, 6).^(f);
valid = isfinite(y1) & isfinite(y2);
y1 = y1(valid); 
y2 = y2(valid);

p = polyfit(y2, y1, 1);
y1_fit = polyval(p, y2);
R2 = 1 - sum((y1 - y1_fit).^2) / sum((y1 - mean(y1)).^2);

plot(y2, y1,'o','Color','k','LineWidth',2)
hold on;
y2_fit_range = linspace(min(y2), max(y2), 100);
plot(y2_fit_range, polyval(p, y2_fit_range), 'r-', 'LineWidth', 2);
xlabel('Topo');
ylabel('Vel');
%%
figure
plot(y1)
hold on 
plot(y2)