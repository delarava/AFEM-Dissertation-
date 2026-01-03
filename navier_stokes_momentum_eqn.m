function [Re] = navier_stokes_momentum_eqn(e,testsp, teste, vars);
% navier_stokes_momentum_eqn -- (FEM Dissertation)
% Compute element residual for Navier-Stokes momentum equation (2D)
%
% Inputs:
%   e        - current element index
%   testsp   - test function space struct (contains .t, .dm, ...)
%   teste    - element / basis object for test functions evaluated on the
%              local element (expects fields: y, ydx, ydy)
%   vars     - a structure containing variables ...
%                  contains vars.[name]e - which contains basis functions
%                                          evaluated for the local element
%
% Output:
%   Re - element residual column vector (testsp.dm * ne x 1)

% Problem parameters (physical)
% Values for these are given in the coursework sheet
rho = 1e-3;
mu = 4e-3;

% Getting the local nodal values of velocity and pressure ...
vl(:,1) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 1);
vl(:,2) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 2);
pl(:,1) = vars.pres.u(vars.pres.t(e,:));

% evaulating the weighted sum of velocity and pressure variables @ quadrature points ...
v = [vars.vele.y(:,:) * vl(:,1) vars.vele.y(:,:) * vl(:,2)];
p = vars.prese.y(:,:) * pl(:,1);

% Velocity gradients
dv1_dx1 = vars.vele.dy(:,:,1) * vl(:,1);  
dv1_dx2 = vars.vele.dy(:,:,2) * vl(:,1);  
dv2_dx1 = vars.vele.dy(:,:,1) * vl(:,2);  
dv2_dx2 = vars.vele.dy(:,:,2) * vl(:,2);  

% Getting local row / local column sizes ...
ne= size(testsp.t,2);
Re = zeros(ne, 1);
for i = 1:ne
    % getting the local element residual indices ordered v_1^{n=1}, v_2^{n=1}, v_1^{n=2}, v_2^{n=2} ...
    vei = (testsp.dm * (i - 1) + 1):(testsp.dm * i);

    % compute advection: ρ(v·∇)v
    advection1 = rho * ( v(:,1).*dv1_dx1 + v(:,2).*dv1_dx2);
    advection2 = rho * ( v(:,1).*dv2_dx1 + v(:,2).*dv2_dx2);

    % Finally, computing weak form ρ(v⋅∇)v⋅w + μ∇v:∇w − p(∇⋅w)
    % Note: ∇v:∇w = ∂v₁/∂x₁ ∂w₁/∂x₁ + ∂v₁/∂x₂ ∂w₁/∂x₂ + ∂v₂/∂x₁ ∂w₂/∂x₁ + ∂v₂/∂x₂ ∂w₂/∂x₂
    Re(vei) =  [dot(teste.gw, advection1 .* teste.y(:,i) ...
        + mu * (dv1_dx1 .* teste.dy(:,i,1) + dv1_dx2 .* teste.dy(:,i,2)) ...
        - p .* teste.dy(:, i, 1));

        dot(teste.gw, advection2 .* teste.y(:,i) ...
        + mu * (dv2_dx1 .* teste.dy(:,i,1) + dv2_dx2 .* teste.dy(:, i, 2)) ...
        - p .* teste.dy(:, i,2))];
end
end
