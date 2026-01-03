%  AFEM solve -- (FEM Tutorials)
%    This function assembles an FEM problem and solves it.

quad=10;
etype = 'triangle';
load('vesseltree0.mat');
%load('vesseltree1.mat');
%load('vesseltree2.mat');
%load('vesseltree3.mat');

nonlinear_tol = 1e-6; % setting the nonlinear tolerance
iter = 0; % initializing the interation counter
residual_vector_norm = 1; % initializing the residual vector norm

% Making Domain function space ...
Omega2.e=fem_get_basis(Omega2.p, quad, etype);
Omega.e=fem_get_basis(Omega.p, quad, etype);

% Making variable function space ...
vel=Omega2;
vel.name = 'vel'; % field name ...

% Setting our initial guess for velocity to be zero ...
vel.u = zeros(vel.dm * size(vel.x,1), 1);
n_vel = size(vel.u,1);

% Making variable function space ...
pres=Omega;
pres.name = 'pres'; % field name ...
pres.dm = 1; % scalar field ...

% Setting our initial guess for pressre to be zero ...
pres.u = zeros(size(pres.x,1),1);
n_pres = size(pres.u,1);

% Starting the Newton-Raphson iteration ...
while (iter == 0) || (residual_vector_norm > nonlinear_tol)
    % Making a FEM object storing my unknown variables (must do this every iteration) ...
    vars.vel = vel;
    vars.pres = pres;

    % Making the residual vector and applying boundary conditions (velocity) ...
    R1 = fem_assemble_block_residual(@navier_stokes_momentum_eqn, Omega, vel, vars);
    for i = 1:size(vel.b,1)

        % Drichlet patch boundaries 2, 3, 4, 5
        if(vel.b(i,5) == 2 || vel.b(i,5) == 3 || vel.b(i,5) == 4 || vel.b(i,5) == 5 )
            continue
        end

        % Looking at the inlet boundary 1
        if (vel.b(i,5) == 1 )
            for j = 2:4
                nn = vel.dm * (vel.b(i,j) - 1);
                R1(nn+1:nn+2) = -[-50 - vel.u(nn+1); 0 - vel.u(nn+2)];
            end
        end

        % Looking at the inlet boundary 6
        if (vel.b(i,5) == 6)
            for j = 2:4
                nn = vel.dm * (vel.b(i,j) - 1);
                R1(nn+1:nn+2) = -[0 - vel.u(nn+1); 0 - vel.u(nn+2)];
            end

        end
    end
    % Making the residual vector and applying boundary conditions (on species B) ...
    R2 = fem_assemble_block_residual(@navier_stokes_mass_eqn, Omega, pres, vars);

    % Making the global residual vector + computing the norm ...
    R = [R1; R2];
    residual_vector_norm = norm(R,2);
    disp(['Current residual (iter=' num2str(iter) '): ' num2str(residual_vector_norm)])
    if(residual_vector_norm < nonlinear_tol)
        continue
    end

    % Creating the discrete matrix operators corresponding to operators a, b, and c
    disp(['        constructing the Jacobian blocks ...'])
    A  = fem_assemble_block_matrix_perturbation(@navier_stokes_momentum_eqn, Omega, vel, vel, vars);
    B  = fem_assemble_block_matrix_perturbation(@navier_stokes_momentum_eqn, Omega, vel, pres, vars);
    C  = fem_assemble_block_matrix_perturbation(@navier_stokes_mass_eqn, Omega, pres, vel, vars);
    D = sparse(size(pres.u,1),size(pres.u,1));

    % Editing block matrices for dirichlet conditions
    for i = 1:size(vel.b,1)
        % Boundary 1
        if(vel.b(i, 5) == 1)
            for j = 2:4
                nn = vel.dm * (vel.b(i, j) - 1);
                % Set velocity Dirichlet in A and B
                A(nn+1:nn+2,:) = 0;
                A(nn+1:nn+2, nn+1:nn+2) = eye(2);
                B(nn+1:nn+2,:) = 0;
            end
        elseif(vel.b(i, 5) == 6)
            for j = 2:4
                nn = vel.dm * (vel.b(i, j) - 1);
                % Set velocity Dirichlet in A and B
                A(nn+1:nn+2,:) = 0;
                A(nn+1:nn+2, nn+1:nn+2) = eye(2);
                B(nn+1:nn+2,:) = 0;
            end
        else
            continue
        end
    end

    % Composing the Jacobian from our Jacobian blocks ...
    disp(['        assembly of the Global Jacobian ...'])
    J = [ A B; C D];

    % Apply Newton Raphson ...
    disp(['        solving for the NR update ...'])
    U = J \ R;
    vel.u(1:n_vel) = vel.u(1:n_vel) - U(1:n_vel);
    pres.u(1:n_pres) = pres.u(1:n_pres) - U(n_vel+1:end);

    % Update the iteration counter ...
    iter = iter + 1;

    disp(['  ']) % skip a line so output is prettier
end

% Plotting the velocity field ...
quiver(vel.x(:,1), vel.x(:,2), vel.u(1:vel.dm:n_vel), vel.u(2:vel.dm:n_vel))

% plotting the pressure field ...
figure(2); trisurf(pres.t(:,1:3), pres.x(:,1), pres.x(:,2), pres.u,'Facecolor','interp','LineStyle','none')
hold; trisurf(pres.t(:,1:3), pres.x(:,1), pres.x(:,2), pres.u,'Facecolor','interp','LineStyle','none')

% Answering the project question:
fprintf('\n');

inlet_id   = 1;
outlet_ids = [2, 3, 4, 5];
wall_id    = 6;   
results = vascular_postprocess(vel, pres, inlet_id, outlet_ids, wall_id);

% 1) Calculating Flow
fprintf('1) FLOW (2D flux units: mm^2/s; multiply by thickness to get mm^3/s)\n');
fprintf('Inlet (Boundary %d):  Length = %.4f mm,  Flux Q_in   = % .6f\n', ...
        results.inletPatch, results.L_in, results.Q_in);

if ~isempty(results.wallPatch)
    fprintf('Walls (Boundary %d):  Length = %.4f mm,  Flux Q_wall = % .6e \n', ...
            results.wallPatch, results.L_wall, results.Q_wall);
else
    fprintf('\n');
end
Qout_total = 0;
for k = 1:numel(results.outletPatches)
    bid = results.outletPatches(k);
    fprintf('Outlet (Boundary %d): Length = %.4f mm,  Flux Q = % .6f\n', ...
            bid, results.L_out(k), results.Q_out(k));
    Qout_total = Qout_total + results.Q_out(k);
end
fprintf('Total outlet flux (sum of 2–5): Q_out_total = %.6f\n\n', Qout_total);

% 2) Outflow distribution
abs_out = abs(results.Q_out);
abs_sum = sum(abs_out);

fprintf('2) OUTFLOW DISTRIBUTION (percent of total outlet flux)\n');

for k = 1:numel(results.outletPatches)
    bid = results.outletPatches(k);
    pct = 100 * abs_out(k) / (abs_sum + eps);

    fprintf('Outlet (Boundary %d): Length = %.4f mm,  Flux |Q| = % .6f   (%5.2f%%)\n', ...
        bid, results.L_out(k), abs_out(k), pct);
end

fprintf('Total outlet flux (sum of 2–5): |Q_out| = %.6f   (100.00%%)\n\n', abs_sum);

% 3) Pressure stats
fprintf('3) PRESSURE STATISTICS\n');
fprintf('Inlet (Boundary %d):   mean p_in    = %.6f\n', ...
        results.inletPatch, results.p_in_mean);

for k = 1:numel(results.outletPatches)
    bid = results.outletPatches(k);
    fprintf('Outlet (Boundary %d):  mean p_out_%d = %.6f\n', ...
            bid, bid, results.p_out_mean_each(k));
end

fprintf('Global outlet mean pressure (length-weighted over 2–5): p_out = %.6f\n', ...
        results.p_out_mean_all);

fprintf('Total pressure drop ΔP = p_in - p_out = %.6f\n\n', results.deltaP);

% 4) Error analysis: Mass conservation check
fprintf('4) MASS CONSERVATION CHECK\n');
fprintf('|Q_in|   = %.6f\n', abs(results.Q_in));
fprintf('|Q_out|  = %.6f\n', abs(Qout_total));
fprintf('Difference |Q_in|-|Q_out| = %.6f  (%.2f%%)\n', ...
        results.mass_error_abs, results.mass_error_relpc);
fprintf('Signed net flux over all open boundaries = %.6e \n', ...
        results.mass_balance_all);
