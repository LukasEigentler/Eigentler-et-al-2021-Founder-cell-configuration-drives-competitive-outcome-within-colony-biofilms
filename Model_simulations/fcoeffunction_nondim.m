function f = fcoeffunction_nondim(location,state)
%FOCOEFFFUNCTION_NONDIM provides the f matrix for PDEToolbox.

% $Author: Lukas Eigentler
% Last updated: 08/06/2021
global kcap kmax c sigma
N = 2; % Number of equations
nr = length(location.x); % Number of columns
f = zeros(N,nr); % Allocate f
f(1,:) = state.u(1,:).*(1-(state.u(1,:)+state.u(2,:))/kcap)-sigma*state.u(1,:).*state.u(2,:);
f(2,:) = kmax*state.u(2,:).*(1-(state.u(1,:)+state.u(2,:))/kcap)-sigma*c*state.u(1,:).*state.u(2,:);
end