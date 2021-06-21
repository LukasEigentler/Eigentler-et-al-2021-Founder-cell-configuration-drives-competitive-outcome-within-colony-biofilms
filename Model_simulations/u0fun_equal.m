function uinit = u0fun_equal(location)
%U0FUN_EQUAL defines a piecewise spatially homogeneous initial condition
% 
% U0_FUN_EQUAL(location) defines the model's initial condition on the
% spatial mesh based on the locations of the mesh's nodes given by the
% location vector

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 21/06/2021

global kcap
M = length(location.x);
xnorm = sqrt(location.x.^2+location.y.^2);
uinit = zeros(2,M);
uinit(1,xnorm<2) = 0.5*kcap;
uinit(2,xnorm<2) = 0.5*kcap;
end
