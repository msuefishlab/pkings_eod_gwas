function det = p0_detectors(wave, s_rate, iP1, iZC1)
%P0_DETECTORS  Quantify presence/strength of the P0 phase without relying on
%the fixed-window area (aP0) or single min-voltage (vP0) used by
%standard_eod_measurement.m.
%
%Two independent, uncertainty-aware detectors are returned so they can be
%cross-checked downstream (see code/02_Phenotyping/p0_phenotype_assessment.Rmd):
%
%  Detector 1 -- derivative / turning-point.  A genuine third phase is a real
%  head-negative deflection (a local minimum) in the trace just before P1.
%  We score its topographic PROMINENCE (how far it stands out from the
%  surrounding trace) and its CURVATURE (2nd-derivative magnitude, i.e. how
%  sharply defined it is).  BP -> no dip (~0); TP -> clear dip; wobble -> weak.
%
%  Detector 2 -- parametric decomposition.  We fit the P0..P1 segment as a sum
%  of Gaussians with and without a P0 component and compare by BIC.  The fitted
%  P0 amplitude comes with a standard error, so a CI straddling 0 flags an
%  ambiguous individual.
%
%Inputs:
%  wave   1 x n_pts, peak-to-peak-normalized, P1-centered waveform (as fed to
%         standard_eod_measurement.m -- a single sweep or an average).
%  s_rate sample rate (Hz).
%  iP1    index of P1 (global max), from standard_eod_measurement.
%  iZC1   index of the zero-crossing just before P1, from standard_eod_measurement.
%
%Output struct det with fields:
%  p0_prominence      prominence of the strongest pre-P1 dip (0 if none)
%  p0_curvature       |2nd derivative| at that dip (0 if none)
%  p0_turning_present logical: a significant, negative-going dip exists
%  p0_fit_amp         |fitted P0 Gaussian amplitude| from the 6-param model
%  p0_fit_se          standard error of the P0 amplitude (NaN if fit failed)
%  p0_deltaBIC        BIC(no P0) - BIC(with P0); positive => P0 warranted
%  fit_r2             R^2 of the with-P0 model over the fitted segment

wave = wave(:).';                 % force row
n_pts = numel(wave);

det = struct('p0_prominence',0,'p0_curvature',0,'p0_turning_present',false, ...
             'p0_fit_amp',NaN,'p0_fit_se',NaN,'p0_deltaBIC',NaN,'fit_r2',NaN);

if isempty(iP1) || isempty(iZC1) || iP1 < 3 || iZC1 < 3
    return
end

% Search window: from ~0.8 ms before ZC1 up to ZC1 (the head-negative P0 dip
% precedes the zero crossing on the P1 upswing).
width = max(5, fix(s_rate * 0.0008));
rstart = max(2, iZC1 - width);
rend   = min(n_pts - 1, iZC1);
if rend - rstart < 3
    return
end
region  = rstart:rend;
yregion = wave(region);

% Baseline noise from the pre-EOD samples; require the dip to clear it and to
% actually go negative (a real head-negative phase).
noise = std(wave(1:min(40, n_pts)));
thr   = 3 * noise;

%% ---- Detector 1: turning-point prominence + curvature ----
% Local minima = peaks of the negated trace; prominence measures how much the
% dip stands out from the surrounding baseline (robust to DC offset).
promBest = 0; locBest = [];
try
    [~, locs, ~, proms] = findpeaks(-yregion);
    if ~isempty(proms)
        [promBest, k] = max(proms);
        locBest = region(1) + locs(k) - 1;    % absolute index of the dip
    end
catch
    % findpeaks unavailable: fall back to the simple minimum of the region
    [~, k] = min(yregion);
    locBest = region(1) + k - 1;
    promBest = max(0, mean(yregion) - yregion(k));
end

if ~isempty(locBest)
    isNeg = wave(locBest) < 0;
    det.p0_prominence      = promBest;
    det.p0_turning_present = (promBest > thr) && isNeg;

    % Curvature: 2nd derivative of a lightly smoothed trace at the dip,
    % scaled by s_rate^2 so it is comparable across recordings.
    ys = movmean(wave, 3);
    d2 = gradient(gradient(ys)) * s_rate^2;
    det.p0_curvature = abs(d2(locBest));
end

%% ---- Detector 2: single bounded Gaussian on the DETRENDED P0 region ----
% We fit ONLY the pre-ZC1 region (the P0 dip), not P1.  Fitting P1 and P0
% together in this short window is unidentifiable -- two large Gaussians trade
% off and the amplitude runs away.  Instead we remove the P1 upswing with a
% linear detrend (a line through the segment endpoints), leaving the P0 dip as
% a residual bump, then fit ONE Gaussian to that residual.  Its amplitude is
% well-determined by the depth of the dip and is bounded to the normalized
% waveform's physical range (|amp| <= 1).
fstart = max(1, iZC1 - width);
fend   = min(n_pts, iZC1);
if fend - fstart < 6
    return
end
idx = fstart:fend;
t   = (idx - iP1) / s_rate * 1000;    % ms, P1 at t=0
y   = wave(idx);

pline = polyfit(t([1 end]), y([1 end]), 1);   % line through the endpoints
r     = y - polyval(pline, t);                 % residual = the P0 dip as a bump
nf    = numel(r);

gauss = @(a, mu, sig, tt) a .* exp(-((tt - mu).^2) ./ (2 * max(abs(sig), 1e-3).^2));
% Soft bounds via penalty (keeps this in base MATLAB, no Optimization Toolbox):
% |amp| <= 1 (physical) and sigma in a plausible 0.02-0.4 ms range.
penalty = @(p) 1e6 * ((abs(p(1)) > 1) * (abs(p(1)) - 1)^2 + ...
                      (abs(p(3)) > 0.4) * (abs(p(3)) - 0.4)^2);
objB = @(p) sum((r - gauss(p(1), p(2), p(3), t)).^2) + penalty(p);

[vdip, kdip] = min(r);
p0B  = [min(vdip, -eps), t(kdip), 0.08];
opts = optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000);

try
    pB   = fminsearch(objB, p0B, opts);
    sseB = sum((r - gauss(pB(1), pB(2), pB(3), t)).^2);   % with-P0 model
    sseA = sum(r.^2);                                     % no-P0 (flat residual)

    BIC = @(s, k) nf * log(max(s, eps) / nf) + k * log(nf);
    det.p0_deltaBIC = BIC(sseA, 0) - BIC(sseB, 3);        % >0 => P0 warranted
    det.p0_fit_amp  = min(abs(pB(1)), 1);                 % bounded amplitude
    det.fit_r2 = 1 - sseB / max(sum((r - mean(r)).^2), eps);

    % Standard error of the P0 amplitude via the numerical Jacobian.
    modelB = @(p, tt) gauss(p(1), p(2), p(3), tt);
    J = numeric_jacobian(modelB, pB, t);
    mse = sseB / max(nf - 3, 1);
    JtJ = J.' * J;
    if rcond(JtJ) > 1e-12
        C = mse * inv(JtJ); %#ok<MINV>
        det.p0_fit_se = sqrt(abs(C(1, 1)));
    end
catch
    % leave fit fields at their NaN defaults
end
end

function J = numeric_jacobian(model, p, t)
%Central-difference Jacobian of model(p,t) w.r.t. each parameter.
n = numel(t);
m = numel(p);
J = zeros(n, m);
for i = 1:m
    h = 1e-6 * max(abs(p(i)), 1e-3);
    pp = p; pm = p;
    pp(i) = pp(i) + h;
    pm(i) = pm(i) - h;
    J(:, i) = (model(pp, t) - model(pm, t)).' / (2 * h);
end
end
