function results = vascular_postprocess(vel, pres, inletPatch, outletPatches, wallPatch)
% VASCULAR_POSTPROCESS
% Post-processing for the vascular tree Navier–Stokes solve.
%
% Computes:
%   - flux at inlet, each outlet, and optionally wall:
%         Q = ∫_Γ v·n ds
%   - boundary lengths
%   - mean pressure on inlet and on each outlet patch:
%         p̄ = (∫_Γ p ds) / (∫_Γ ds)
%   - global outlet-mean pressure and ΔP = p_in - p_out
%   - mass-conservation diagnostics

% -----------------------------
% Geometry / connectivity
% -----------------------------
b   = vel.b;      % boundary array
T   = vel.t;      % element connectivity for velocity
X   = vel.x;      % coordinates

patchCol = 5;     % your mesh stores patch ID in column 5

nOutlet = numel(outletPatches);
useWall = (nargin >= 5) && ~isempty(wallPatch);

% Storage
Q_out   = zeros(nOutlet,1);
L_out   = zeros(nOutlet,1);

Q_in    = 0.0;
L_in    = 0.0;

Q_wall  = 0.0;
L_wall  = 0.0;

% pressure integrals
Pin_num = 0.0; Pin_den = 0.0; % inlet
Pout_num_k = zeros(nOutlet,1); % per outlet
Pout_den_k = zeros(nOutlet,1);

Pout_num_all = 0.0; % combined outlets
Pout_den_all = 0.0;

domain_centroid = mean(X,1);  % for outward normal

% Loop over boundary segments
for k = 1:size(b,1)

    elem  = b(k,1);
    nA    = b(k,2);
    nB    = b(k,3);
    patch = b(k,patchCol);

    if nA <= 0 || nB <= 0
        continue;
    end

    isInlet  = (patch == inletPatch);
    isOutlet = any(patch == outletPatches);
    isWall   = useWall && (patch == wallPatch);

    if ~(isInlet || isOutlet || isWall)
        continue;
    end

    % 1) edge geometry ---
    xA = X(nA,:); xB = X(nB,:);
    tvec = xB - xA;
    L = norm(tvec);
    if L < 1e-12
        continue;
    end

    % candidate outward normal
    n = [tvec(2), -tvec(1)] / L;
    mid = 0.5*(xA + xB);
    if dot(n, mid - domain_centroid) < 0
        n = -n;
    end

    % local face index inside this triangle
    te = T(elem,:);
    ia = find(te == nA, 1);
    ib = find(te == nB, 1);
    face = local_face_triangle(ia, ib);

    % 2) quadrature weights on this physical edge 
    fw      = vel.e.f{face}.gw;      % weights on reference edge
    masterL = master_edge_len(face); % |reference edge|
    w_ds    = fw * (L / masterL);    % physical ds weights

    % 3) velocity at quad points 
    fyV = vel.e.f{face}.y;
    tve = vel.t(elem,:);
    v1l = vel.u(2*(tve-1) + 1);
    v2l = vel.u(2*(tve-1) + 2);
    vx_q = fyV * v1l;
    vy_q = fyV * v2l;
    vn_q = vx_q*n(1) + vy_q*n(2);

    % 4) pressure at quad points
    fyP = pres.e.f{face}.y;
    tpe = pres.t(elem,:);
    pl  = pres.u(tpe);
    p_q = fyP * pl;

    % 5) integrate on this edge
    flux_edge   = sum(w_ds .* vn_q);
    p_int_edge  = sum(w_ds .* p_q);
    L_edge      = L;
    w_sum_edge  = sum(w_ds);

    if isInlet
        Q_in    = Q_in    + flux_edge;
        L_in    = L_in    + L_edge;
        Pin_num = Pin_num + p_int_edge;
        Pin_den = Pin_den + w_sum_edge;
    end

    if isOutlet
        idx = find(outletPatches == patch, 1);
        Q_out(idx)      = Q_out(idx)      + flux_edge;
        L_out(idx)      = L_out(idx)      + L_edge;
        Pout_num_k(idx) = Pout_num_k(idx) + p_int_edge;
        Pout_den_k(idx) = Pout_den_k(idx) + w_sum_edge;

        Pout_num_all = Pout_num_all + p_int_edge;
        Pout_den_all = Pout_den_all + w_sum_edge;
    end

    if isWall
        Q_wall = Q_wall + flux_edge;
        L_wall = L_wall + L_edge;
    end
end

% 6) Derived quantities
p_in_mean        = Pin_num / (Pin_den + eps);
p_out_mean_all   = Pout_num_all / (Pout_den_all + eps);
p_out_mean_each  = Pout_num_k ./ (Pout_den_k + eps);

Qout_total       = sum(Q_out);

absQ             = abs(Q_out);
absSum           = sum(absQ);
if absSum > 0
    outletFrac = absQ / absSum;
else
    outletFrac = zeros(size(absQ));
end

mass_balance_all = Q_in + Qout_total + (useWall * Q_wall);  % signed
mass_error_abs   = abs(Q_in) - abs(Qout_total);
mass_error_relpc = 100 * abs(mass_error_abs) / (abs(Q_in) + eps);

% Pack results
results.inletPatch      = inletPatch;
results.outletPatches   = outletPatches(:);
if useWall
    results.wallPatch   = wallPatch;
else
    results.wallPatch   = [];
end

results.Q_in            = Q_in;
results.Q_out           = Q_out;
results.Q_wall          = Q_wall;

results.L_in            = L_in;
results.L_out           = L_out;
results.L_wall          = L_wall;

results.p_in_mean       = p_in_mean;
results.p_out_mean_all  = p_out_mean_all;
results.p_out_mean_each = p_out_mean_each;

results.deltaP          = p_in_mean - p_out_mean_all;
results.outletFractions = outletFrac;

results.mass_balance_all = mass_balance_all;
results.mass_error_abs   = mass_error_abs;
results.mass_error_relpc = mass_error_relpc;

% aliases so your older printing code still works
results.p_out_mean   = p_out_mean_all;
results.mass_balance = mass_balance_all;

end

% helpers

function face = local_face_triangle(i,j)
pair = sort([i j]);
if     isequal(pair,[1 2]), face = 1;
elseif isequal(pair,[1 3]), face = 2;
elseif isequal(pair,[2 3]), face = 3;
else
    error('Invalid triangle edge (%d,%d)', i, j);
end
end

function L = master_edge_len(face)
% reference P1 triangle: edges 1 & 2 length 1, edge 3 length sqrt(2)
if face == 3
    L = sqrt(2);
else
    L = 1.0;
end
end
