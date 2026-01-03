function [Re] = navier_stokes_mass_eqn(e,testsp, teste, vars);
%  navier_stokes_mass_eqn -- (FEM Dissertation)
%
% This function computes the element matrix for the 2D discrete Laplacian
%
% Inputs --
%
%        e       - the current element (integer)
%
%        testsp  - the space of test functions ...
%
%        teste   - the element / basis object for test functions
%                  evaluated for the local element ...
%
%        vars    - a structure containing variables ...
%                  contains vars.[name]e - which contains basis functions
%                                          evaluated for the local element
%
%
% Outputs --
%
%        Re - the element residual for Darcy mass equation ...
%
% by David Nordsletten
% Jan 2015
%

% Problem parameters

% Getting the local nodal values of velocity and pressure ...
v(:,1) = vars.vel.u(vars.vel.dm * (vars.vel.t(e,:)-1) + 1);
v(:,2) = vars.vel.u(vars.vel.dm * (vars.vel.t(e,:)-1) + 2);

% evaulating the weighted sum of divergence of velocity @ quadrature points ...
divv = vars.vele.dy(:,:,1) * v(:,1) + vars.vele.dy(:,:,2) * v(:,2);

% Getting local row / local column sizes ...
ne = size(testsp.t,2);
Re = zeros(ne, 1); % (testsp.dm was just 1 so just do (ne,1)) 
% Pressure space has one degree of freedom per node

for i = 1:ne
    % computing ∫ (div v) * q dΩ
    Re(i) = dot(teste.gw, divv .* teste.y(:,i));
end


