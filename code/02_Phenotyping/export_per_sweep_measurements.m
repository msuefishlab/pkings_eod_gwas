%EXPORT_PER_SWEEP_MEASUREMENTS
%Re-measures the P0 phase on EVERY replicate sweep (not just the per-fish
%average) so that per-individual measurement uncertainty can be quantified,
%and applies the two uncertainty-aware P0 detectors in p0_detectors.m.
%
%The existing pipeline (EOD_Analysis -> PlotEODs) averages all sweeps of a
%fish before measuring, discarding within-individual variability and
%attenuating small P0 phases (sweeps are centered on P1, not P0).  This
%script keeps the replicates: EOD_Analysis.m already loads and normalizes
%them into normalized_eods(i).wave (points x sweeps).
%
%PREREQUISITE: run EOD_Analysis first so `normalized_eods` and
%`averaged_eods` exist in the workspace.  Then run this script.
%
%Outputs (output_data/02_Phenotyping/):
%  measurement_data_per_sweep.csv  one row per sweep: legacy scalars + both detectors
%  p0_windows.csv                  one row per fish: aligned, z-scored, P0 shape window
%  p0_windows_per_sweep.csv        one row per SWEEP: aligned, z-scored P0 window
%  p0_baseline.csv                 one row per fish: flat pre-P1 baseline strip (averaged wave)
%  p0_baseline_per_sweep.csv       one row per SWEEP: flat pre-P1 baseline strip
%  p0_aligned_per_sweep.csv        one row per SWEEP: P1-aligned, ptp-normalized trace
%
%p0_aligned_per_sweep.csv holds the raw (peak-to-peak-normalized, NOT z-scored)
%waveform on a fixed P1-relative ms axis, so amplitudes are comparable ACROSS
%sweeps.  This is what evoked-response detectors (F_sp / plus-minus reference)
%need: each replicate sweep is a trial of the same discharge, and the cross-sweep
%variance at each time point is the noise floor.  (The z-scored p0_windows_*
%files destroy that cross-sweep amplitude relationship and cannot be used for it.)
%
%The two baseline strips capture a stretch of flat trace WELL BEFORE the P0 zone
%(no P0 present by construction).  Downstream (p0_statistical_detection.Rmd) the
%SAME P0 detector is slid across these strips to build a per-fish empirical NULL,
%turning the P0-window score into a calibrated p-value -- each fish is its own
%noise control, so nothing has to be thresholded by eye.  The strip is
%mean-subtracted only (NOT z-scored); R z-scores each sliding sub-window so the
%null gets the identical treatment to the P0 window.
%
%The per-sweep windows are what let the R report detect a real-but-tiny P0 by
%its CONSISTENCY across replicates (a dip that recurs at the same place across
%sweeps is real even when each sweep is near the noise floor) rather than from a
%single averaged trace, where a small P1-aligned (not P0-aligned) dip smears out.
%
%The per-fish measurement_data.csv written by PlotEODs.m is left untouched.

root = get_root();
outdir = fullfile(root, 'output_data', '02_Phenotyping');
if ~exist(outdir, 'dir'); mkdir(outdir); end

% Shape window: ~1 ms so a real P0 dip is reliably inside it (a 0.5 ms window,
% anchored on the fragile iZC1, could clip a slightly early dip).
WIN_MS   = 0.0010;      % window duration ending at iZC1 (seconds)
N_target = 256;         % resampled length of the per-fish averaged window
N_sweepw = 128;         % resampled length of the per-sweep windows (smaller file)

% Flat baseline strip for the per-fish empirical null: a stretch WELL BEFORE the
% P0 zone (5.0 -> 2.8 ms before P1, safely left of the P0).  Resampled to N_base so
% its resolution is 128 samples/ms -- identical to the per-sweep P0 window -- so a
% 1 ms sliding sub-window in R is exactly N_sweepw points at native resolution.
BASE_MS0 = 0.0050;      % strip start: ms before P1
BASE_MS1 = 0.0028;      % strip end:   ms before P1 (2.8 ms: the P0 reaches ~1.5 ms
                        % before P1, so 1.5 ms clipped the P0; keep the strip flat)
N_base   = 448;         % 3.5 ms * 128 samples/ms

% P1-aligned raw trace for evoked-response detectors (F_sp / plus-minus reference).
% Fixed window RELATIVE TO P1 (negative = before), ptp-normalized, NOT z-scored, so
% cross-sweep amplitude is preserved.  Spans the P0 zone plus flat baseline.
ALIGN_MS0 = -0.0030;    % window start rel. P1 (s): 3.0 ms before P1
ALIGN_MS1 =  0.0005;    % window end   rel. P1 (s): 0.5 ms after P1
N_align   = 448;        % 3.5 ms * 128 samples/ms
nfish = numel(averaged_eods);

sweep_rows = {};        % accumulate per-sweep measurements
sweepw_rows = {};       % accumulate per-sweep P0 windows
base_rows   = {};       % accumulate per-sweep baseline strips
align_rows  = {};       % accumulate per-sweep P1-aligned raw traces
win_names  = strings(nfish, 1);
win_mat    = nan(nfish, N_target);
base_names = strings(nfish, 1);
base_mat   = nan(nfish, N_base);

% Helper: extract the iZC1-anchored window, z-score, resample to n points.
% Returns [] if the window is too short. (Local function at end of file.)

for i = 1:nfish
    name = string(averaged_eods(i).sample_name);
    sr   = averaged_eods(i).sampRate;

    %% ---- per-sweep measurement + detectors ----
    W = normalized_eods(i).wave;
    [~, nsw] = size(W);
    for j = 1:nsw
        sweep = W(:, j).';                       % row vector for standard_eod_measurement
        if isempty(sweep) || all(sweep == 0) || all(isnan(sweep))
            continue
        end
        try
            md = standard_eod_measurement(sweep, sr, name);
        catch
            continue                             % skip un-measurable sweeps
        end
        det = p0_detectors(sweep, sr, md.iP1, md.iZC1);

        sweep_rows(end+1, :) = { name, j, ...
            md.aP0, md.vP0, md.tP0, md.iZC1, ...
            det.p0_prominence, det.p0_curvature, double(det.p0_turning_present), ...
            det.p0_fit_amp, det.p0_fit_se, det.p0_deltaBIC, det.fit_r2 }; %#ok<SAGROW>

        % Per-sweep P0 shape window (for consistency scoring in R).
        sw_win = p0_window(sweep, md.iZC1, fix(sr * WIN_MS), N_sweepw);
        if ~isempty(sw_win)
            sweepw_rows(end+1, :) = [{name, j}, num2cell(sw_win(:).')]; %#ok<SAGROW>
        end

        % Per-sweep flat baseline strip (feeds the per-fish empirical null in R).
        sw_base = p0_baseline(sweep, md.iP1, fix(sr * BASE_MS0), fix(sr * BASE_MS1), N_base);
        if ~isempty(sw_base)
            base_rows(end+1, :) = [{name, j}, num2cell(sw_base(:).')]; %#ok<SAGROW>
        end

        % Per-sweep P1-aligned raw trace (feeds F_sp / plus-minus in R).
        sw_align = p0_aligned(sweep, md.iP1, ALIGN_MS0, ALIGN_MS1, sr, N_align);
        if ~isempty(sw_align)
            align_rows(end+1, :) = [{name, j}, num2cell(sw_align(:).')]; %#ok<SAGROW>
        end
    end

    %% ---- aligned P0 shape window (from the averaged wave) ----
    try
        mavg = standard_eod_measurement(averaged_eods(i).wave, sr, name);
        av_win = p0_window(averaged_eods(i).wave, mavg.iZC1, fix(sr * WIN_MS), N_target);
        if ~isempty(av_win)
            win_mat(i, :) = av_win;
        end
        av_base = p0_baseline(averaged_eods(i).wave, mavg.iP1, ...
                              fix(sr * BASE_MS0), fix(sr * BASE_MS1), N_base);
        if ~isempty(av_base)
            base_mat(i, :) = av_base;
        end
    catch
        % leave this fish's window/baseline as NaN
    end
    win_names(i)  = name;
    base_names(i) = name;
end

%% ---- write per-sweep table ----
sweep_vars = {'sample_name','sweep_index','aP0','vP0','tP0','iZC1', ...
              'p0_prominence','p0_curvature','p0_turning_present', ...
              'p0_fit_amp','p0_fit_se','p0_deltaBIC','fit_r2'};
Tsweep = cell2table(sweep_rows, 'VariableNames', sweep_vars);
writetable(Tsweep, fullfile(outdir, 'measurement_data_per_sweep.csv'));

%% ---- write per-fish averaged P0 window table ----
Twin = array2table(win_mat, 'VariableNames', compose('p%03d', 1:N_target));
Twin = addvars(Twin, win_names, 'Before', 1, 'NewVariableNames', 'sample_name');
writetable(Twin, fullfile(outdir, 'p0_windows.csv'));

%% ---- write per-sweep P0 window table ----
sw_vars = [{'sample_name', 'sweep_index'}, compose('p%03d', 1:N_sweepw)];
Tsw = cell2table(sweepw_rows, 'VariableNames', sw_vars);
writetable(Tsw, fullfile(outdir, 'p0_windows_per_sweep.csv'));

%% ---- write per-fish averaged baseline strip table ----
Tbase = array2table(base_mat, 'VariableNames', compose('b%03d', 1:N_base));
Tbase = addvars(Tbase, base_names, 'Before', 1, 'NewVariableNames', 'sample_name');
writetable(Tbase, fullfile(outdir, 'p0_baseline.csv'));

%% ---- write per-sweep baseline strip table ----
base_vars = [{'sample_name', 'sweep_index'}, compose('b%03d', 1:N_base)];
Tswbase = cell2table(base_rows, 'VariableNames', base_vars);
writetable(Tswbase, fullfile(outdir, 'p0_baseline_per_sweep.csv'));

%% ---- write per-sweep P1-aligned raw trace table ----
align_vars = [{'sample_name', 'sweep_index'}, compose('a%03d', 1:N_align)];
Talign = cell2table(align_rows, 'VariableNames', align_vars);
writetable(Talign, fullfile(outdir, 'p0_aligned_per_sweep.csv'));

fprintf(['Wrote %d sweep rows (%d sweep windows, %d baseline strips, ' ...
         '%d aligned traces) for %d fish to %s\n'], ...
        height(Tsweep), height(Tsw), height(Tswbase), height(Talign), nfish, outdir);


function w = p0_window(wave, iZC1, width, npts)
%Extract the iZC1-anchored P0 window, z-score it (amplitude-invariant), and
%resample to npts. Returns a 1 x npts row, or [] if the window is too short.
    wave = wave(:);
    i1 = max(1, iZC1 - width);
    i2 = min(numel(wave), iZC1);
    seg = wave(i1:i2);
    if numel(seg) < 3
        w = [];
        return
    end
    seg = seg - mean(seg);
    s = std(seg);
    if s > 0; seg = seg / s; end
    w = interp1(linspace(0, 1, numel(seg)), seg, linspace(0, 1, npts), 'linear');
end


function b = p0_baseline(wave, iP1, off0, off1, npts)
%Extract a flat pre-P1 baseline strip spanning [iP1-off0, iP1-off1] (samples),
%WELL LEFT of the P0 zone, mean-subtract only (R z-scores each sliding sub-window
%to match the P0-window treatment), and resample to npts. Returns [] if too short.
    wave = wave(:);
    i1 = max(1, iP1 - off0);
    i2 = max(1, iP1 - off1);
    if i2 - i1 < 10
        b = [];
        return
    end
    seg = wave(i1:i2);
    seg = seg - mean(seg);
    b = interp1(linspace(0, 1, numel(seg)), seg, linspace(0, 1, npts), 'linear');
end


function a = p0_aligned(wave, iP1, ms0, ms1, sr, npts)
%Extract a fixed P1-relative segment [iP1+ms0*sr, iP1+ms1*sr] (ms0<0 = before P1),
%WITHOUT mean-subtraction or z-scoring, and resample to npts. The sweep is already
%ptp-normalized and locally baselined upstream, so amplitudes are comparable across
%sweeps -- required for the F_sp / plus-minus cross-sweep variance. [] if OOB.
    wave = wave(:);
    i1 = iP1 + round(ms0 * sr);
    i2 = iP1 + round(ms1 * sr);
    if i1 < 1 || i2 > numel(wave) || i2 - i1 < 10
        a = [];
        return
    end
    seg = wave(i1:i2);
    a = interp1(linspace(0, 1, numel(seg)), seg, linspace(0, 1, npts), 'linear');
end


function gitRoot = get_root()
    [status, cmdout] = system('git rev-parse --show-toplevel');
    if status == 0
        gitRoot = strtrim(cmdout);
    else
        error('Not a Git repository or Git is not installed.');
    end
end
