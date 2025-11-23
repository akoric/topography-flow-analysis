%% parameter_calculation.m
% Calculate spectra and statistical metrics for DEM files
% within specified frame ranges and save the results.

clearvars; close all;

% -------------------------------------------------------------------------
% --- CHANGE THESE NAMES FIRST ---
% -------------------------------------------------------------------------
varList = { ...
    'U_kurt','U_mean','U_skew','U_std','U_var', ...
    'V_kurt',...
    'V_mean',...
    'V_skew'...
    'V_std','V_var', ...
    'W_kurt', ...
    'W_mean','W_skew','W_std','W_var'
    };

saveFolder = '/path/to/save/results_newnew';

% -------------------------------------------------------------------------
% Loop over variables
% -------------------------------------------------------------------------
% 4435 3223
hurstIdx_perDEM = {1:2,1:2,1:2,1:2};

for v = 1:numel(varList)
    varName  = varList{v};                     % e.g. 'U_kurt'
    saveName = [varName '_results.mat'];    % e.g. 'U_kurt_results.mat'

    % DEM file list
    demFiles = {
        'DEM1_Moments.mat', ...
        'DEM2_Moments.mat', ...
        'DEM3_Moments.mat', ...
        'DEM4_Moments.mat' ...
        };

    % Corresponding frame ranges [start, end]
    idxRanges = [ ...
        21, 49;  % DEM1
        27, 49;  % DEM2
        32, 49;  % DEM3
        17, 46   % DEM4
        ];
    % For single-frame processing, you can use:
    % idxRanges = [1,1; 1,1; 1,1; 1,1];

    % Grid resolution (m)
    res = 0.005;

    % -------------------------------------------------------------------------
    % Initialization
    % -------------------------------------------------------------------------
    allResults = [];   % Store all computed results
    demID      = [];   % DEM file index
    frameID    = [];   % Frame index

    % -------------------------------------------------------------------------
    % Main loop (to be implemented)
    % -------------------------------------------------------------------------
    for k = 1:numel(demFiles) % each topograph
        % change S and U to the poperty to be loaded
        S = load(demFiles{k}, varName);
        U = S.(varName);  % Ny×Nx×Nt
        i0 = idxRanges(k,1);
        i1 = min(idxRanges(k,2));

        for i = i0:i1  % each frame
            % Extract frame
            img0 = U(1:600,1:200,i);
            
            nanRows = all(isnan(img0),2);  
            nanCols = all(isnan(img0),1);  
                
            img = img0(~nanRows, ~nanCols);

            % — 1) 2D PSD
            [psd2d, kx, ky]  = compute_psd2d(img, res);
            [k_rad, psd_k]  = radial_average(psd2d, res);

            % — 2) Fractal & Hurst (radial) —
            sh  = estimate_slope(k_rad, psd_k,'fractal');
            fract_k = 2+1-(-sh-2)/2;
            sl  = estimate_slope(k_rad, psd_k, 'hurst', hurstIdx_perDEM{k});
            hurst_k = 1-(sl+1)/2;

            % — 3) co% int
            [Lx, lag1, cov_x] = compute_covariance(img,0,res);
            [Ly, lag2, cov_y] = compute_covariance(img,1,res);
            [r_vals, radial_cov] = radial_average_covariance_matlab(img, res);
            index = find(radial_cov < 0.05, 1);
            if isempty(index)
                warning('Covariance never drops below 0.05; integrating full range.');
                index = numel(radial_cov);
            end
            Lk = trapz(r_vals(1:index), radial_cov(1:index));

            % — 4) 1D PSD (X, Y) + Fractal/Hurst —
            [kx1, psd_x] = compute_psd1d(img,0,res);
            [ky1, psd_y] = compute_psd1d(img,1,res);
            [fract_x, hurst_x] = fractal_and_hurst(kx1, psd_x,1);
            [fract_y, hurst_y] = fractal_and_hurst(ky1, psd_y,1);

            % — 5)
            [Ra, Rq, Rq_s, Sk, Ku] = tribo_params(img,res);

            [Ra_sw, Rq_sw, Sk_sw, Ku_sw] = tribo_params_streamwise(img);

            % — collect results —
            row = [fract_x, hurst_x, fract_y, hurst_y, ...
                fract_k, hurst_k, Ra, Rq, Rq_s, Sk, Ku, ...
                Lk, Lx, Ly, Ra_sw, Rq_sw, Sk_sw, Ku_sw];
            allResults = [allResults; row];
            demID      = [demID;   k];
            frameID    = [frameID; i];
        end
    end

    % construct a table corresponding to "row"
    colNames = {'FracD_x','Hurst_x','FracD_y','Hurst_y', ...
        'FracD_k','Hurst_k','Ra','Rq','Rq_slope','Sk','Ku','Lk','Lx','Ly',...
        'Ra_sw','Rq_sw','Sk_sw','Ku_sw'};
    demNames = cellfun(@(f) erase(f,'.mat'), demFiles,'UniformOutput',false)';

    Tdata = array2table(allResults, 'VariableNames', colNames);

    T = [ ...
        table(demID, frameID, 'VariableNames', {'DEM','Frame'}), ...
        Tdata ...
        ];

    % -------------------------------------------------------------------------
    % Save results (.mat and .csv)
    % -------------------------------------------------------------------------
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder); 
    end

    [~, baseName, ~] = fileparts(saveName);

    outMAT = fullfile(saveFolder, [baseName '.mat']);
    outCSV = fullfile(saveFolder, [baseName '.csv']);

    save(outMAT, 'T');
    writetable(T, outCSV);
    fprintf('Saved: %s and %s\n', outMAT, outCSV);
end
%% ====== functions ======
function [r_vals, radial_cov] = radial_average_covariance_matlab(image, res)
    image_zero_mean = image - mean(image(:));
    cov2d = xcorr2(image_zero_mean, image_zero_mean);
    cov2d = cov2d / (var(image(:),1) * numel(image));

    [ny, nx] = size(cov2d);
    center_y = ceil(ny/2);
    center_x = ceil(nx/2);
    [y, x] = ndgrid(1:ny, 1:nx);
    r = sqrt((x - center_x).^2 + (y - center_y).^2);
    r = round(r);

    r_max = max(r(:));
    radial_cov = zeros(r_max + 1, 1);
    counts = zeros(r_max + 1, 1);

    for i = 1:ny
        for j = 1:nx
            idx = r(i, j) + 1;
            radial_cov(idx) = radial_cov(idx) + cov2d(i, j);
            counts(idx) = counts(idx) + 1;
        end
    end

    radial_cov = radial_cov ./ counts;
    r_vals = (0:r_max)' * res;
end


function [psd2d, kx, ky] = compute_psd2d(img, res)
    v = var(img(:));
    F = fftshift(fft2(img));
    psd2d = abs(F).^2;
    kx = fftshift(fftfreq(size(img,2),res));
    ky = fftshift(fftfreq(size(img,1),res));
    inner = trapz(ky, psd2d, 1);
    inner = inner(:);
    I = trapz(kx, inner);
    if I>0, psd2d = psd2d./I*v; end
end

function [k, radial] = radial_average(psd2d, res)
    [Ny,Nx] = size(psd2d);
    [X,Y] = meshgrid(1:Nx,1:Ny);
    R = round(sqrt((X-Nx/2).^2+(Y-Ny/2).^2));
    rmax = max(R(:));
    radial = zeros(rmax+1,1); cnt = radial;
    for r=0:rmax
        m = (R==r);
        radial(r+1)=sum(psd2d(m));
        cnt(r+1)=nnz(m);
    end
    radial = radial./cnt;
    Lx = Nx*res; Ly = Ny*res;
    k = linspace(0,rmax/max(Lx,Ly),rmax+1)';
end

function slope = estimate_slope(k, psd, mode, idx_opt)
    % log-space (skip k = 0)
    kl = log(k(2:end));
    pl = log(psd(2:end));
    switch lower(mode)
        case 'hurst'
            default_idx = 1:3;
        case 'fractal'
            default_idx = 5:25;
        otherwise
            error('mode error');
    end
    if nargin >= 4 && ~isempty(idx_opt)
        if isscalar(idx_opt)
            idx = 1:idx_opt;
        else
            idx = idx_opt(:).';
        end
    else
        idx = default_idx;
    end
    n = numel(kl);
    idx = idx(idx>=1 & idx<=n);
    if isempty(idx) || numel(idx) < 2
        error('estimate_slope: invalid idx range for mode "%s".', mode);
    end

    p = polyfit(kl(idx), pl(idx), 1);
    slope = p(1);
end


function [Ra,Rq,Rqs,Sk,Ku] = tribo_params(img,res)
    m = mean(img(:));
    Ra = mean(abs(img(:)-m));
    Rq = sqrt(mean((img(:)-m).^2));
    [dx,dy] = gradient(img,res);
    Rqs = sqrt(mean(dx(:).^2+dy(:).^2));
    s = std(img(:));
    Sk = mean(((img(:)-m)/s).^3);
    Ku = mean(((img(:)-m)/s).^4);
end

function [Ra_sw, Rq_sw, Sk_sw, Ku_sw] = tribo_params_streamwise(img)
    Ny = size(img,2);
    Ra_sw = zeros(Ny,1); Rq_sw = zeros(Ny,1); Sk_sw = zeros(Ny,1); Ku_sw = zeros(Ny,1);
    for i = 1:Ny
        line = img(:,i);
        m = mean(line);
        s = std(line);
        Ra_sw_1(i) = mean(abs(line - m));
        Rq_sw_1(i) = sqrt(mean((line - m).^2));
        Sk_sw_1(i) = mean(((line - m)/s).^3);
        Ku_sw_1(i) = mean(((line - m)/s).^4);
    end
    Ra_sw = mean(Ra_sw_1, 'omitnan');
    Rq_sw = mean(Rq_sw_1, 'omitnan');
    Sk_sw = mean(Sk_sw_1, 'omitnan');
    Ku_sw = mean(Ku_sw_1, 'omitnan');

end


function [k1, psd1] = compute_psd1d(img, axis, res)
    if axis==0, data = img; else data = img.'; end
    N  = size(data,2);
    F  = fft(data,[],2);
    psd1 = mean(abs(F).^2./N,1);
    k1   = fftfreq(N,res).';
    half=1:floor(N/2);
    psd1 = psd1(half);
    k1   = k1(half);
end

function [FD,H] = fractal_and_hurst(k1,psd1,dim)
    sh = estimate_slope(k1,psd1,'fractal');
    FD = dim+1 -(-sh-1)/2;
    sl = estimate_slope(k1,psd1,'hurst');
    H  = 1-(sl+1)/2;
end


function [L, lags, covf] = compute_covariance(img, axis, res)
    % Compute covariance function, lags, and integral scale L
    if axis == 0
        covf_all = cell(size(img,1), 1);
        lags_all = cell(size(img,1), 1);    
            for i = 1:size(img,1)
                % mean(img(i,:))
                    if any(isnan(img(i,:)))
                        covf_all{i} = NaN;
                        lags_all{i} = NaN;
                        continue;
                    end
                [covf_all{i}, lags_all{i}] = autocorr(img(i,:) - mean(img(i,:)), 'NumLags', size(img,2)-1);
            end  

            % for j = 1:size(img,2)
            %     sum_val = 0;
            %     for i = 1:size(img,1)
            %         sum_val = sum_val + covf_all{i}(j);
            %     end
            %     covf(j) = sum_val / size(img,1);
            % end
            for j = 1:size(img,2)
                vals = zeros(size(img,1),1);
                valid_count = 0;
                for i = 1:size(img,1)
                    if ~isempty(covf_all{i}) && length(covf_all{i}) >= j && ~isnan(covf_all{i}(j))
                        vals(i) = covf_all{i}(j);
                        valid_count = valid_count + 1;
                    else
                        vals(i) = NaN;
                    end
                end
                if valid_count > 0
                    covf(j) = mean(vals(~isnan(vals)));
                else
                    covf(j) = NaN;
                end
            end


    else
    covf_all = cell(size(img,2), 1);
    lags_all = cell(size(img,2), 1);    
            for i = 1:size(img,2)
                    if any(isnan(img(:,i)))
                        covf_all{i} = NaN;
                        lags_all{i} = NaN;
                        continue;
                    end
                [covf_all{i}, lags_all{i}] = autocorr(img(:,i) - mean(img(:,i)), 'NumLags', size(img,1)-1);
            end  

            % for j = 1:size(img,1)
            %     sum_val = 0;
            %     for i = 1:size(img,2)
            %         sum_val = sum_val + covf_all{i}(j);
            %     end
            %     covf(j) = sum_val / size(img,2);
            % end     
            for j = 1:size(img,1)
                vals = zeros(size(img,2),1);
                valid_count = 0;
                for i = 1:size(img,2)
                    if ~isempty(covf_all{i}) && length(covf_all{i}) >= j && ~isnan(covf_all{i}(j))
                        vals(i) = covf_all{i}(j);
                        valid_count = valid_count + 1;
                    else
                        vals(i) = NaN;
                    end
                end
                if valid_count > 0
                    covf(j) = mean(vals(~isnan(vals)));
                else
                    covf(j) = NaN;
                end
            end

    end

    lags = lags_all{1} * res;

    index = find(covf < 0.05, 1);
    if isempty(index)
        warning('Covariance never drops below 0.05; integrating full range.');
        index = numel(covf);
    end
    if any(isnan(covf(1:index))) || any(isnan(lags(1:index)))
        L = NaN;
    else
        L = trapz(lags(1:index), covf(1:index));
    end
    % L = trapz(lags(1:index), covf(1:index));

end



function f = fftfreq(n,d)
    if mod(n,2)==0
        k = [0:n/2-1, -n/2:-1];
    else
        k = [0:(n-1)/2, -(n-1)/2:-1];
    end
    f = k./(n*d);
end


%% visulization of velocity field (example code below kept commented out)
% figure('Position',[100 100 600 650]);
% hAx = axes('Position',[0.05 0.15 0.9 0.8]);
% hImg = imagesc(U_mean(:,:,1),'Parent',hAx);
% axis equal tight; colorbar; colormap(jet);
% title(hAx,'Frame 1 / ' + string(size(U_mean,3)));
%
% hTxt = uicontrol('Style','text', ...
%     'Units','pixels', ...
%     'Position',[260 30 80 20], ...
%     'String','1', ...
%     'FontSize',10);
%
% hSld = uicontrol('Style','slider', ...
%     'Units','pixels', ...
%     'Position',[150 10 300 20], ...
%     'Min',1,'Max',size(U_mean,3),'Value',1, ...
%     'SliderStep',[1/(size(U_mean,3)-1) , 10/(size(U_mean,3)-1)], ...
%     'Callback',@(src,~) onSlide(src, hImg, hTxt, hAx, U_mean));
%
% function onSlide(src, hImg, hTxt, hAx, U_mean)
%     idx = round(src.Value);
%     set(hImg, 'CData', U_mean(:,:,idx));
%     set(hTxt, 'String', sprintf('%d', idx));
%     title(hAx, sprintf('Frame %d / %d', idx, size(U_mean,3)));
% end
%
% figure('Position',[100 100 600 600]);
% imagesc(bed(:,:));
% axis equal tight;
% colormap(jet);
% colorbar;
% title('Bed Elevation Map');
% xlabel('X index');
% ylabel('Y index');
%
% pan on;
% zoom on;
