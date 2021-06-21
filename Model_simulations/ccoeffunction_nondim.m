function cmat = ccoeffunction_nondim(location,state)
%COCOEFFFUNCTION_NONDIM provides the c matrix for PDEToolbox.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 21/06/2021

global d kcap N time_start
% time_check = toc(time_start);
if toc(time_start) > 10000 && length(location.x)>1
    return
end
cmat = zeros(4*N^2,length(location.x));

if length(location.x)>4 % ensure diffusion coefficient is not negative
    kk = find(state.u(1,:)+state.u(2,:)>kcap);
    state.u(1,kk) = kcap*state.u(1,kk)./(state.u(1,kk)+state.u(2,kk));
    state.u(2,kk) = kcap*state.u(2,kk)./(state.u(1,kk)+state.u(2,kk));
end
cmat(1,:) =  (1-(state.u(1,:)+state.u(2,:))/kcap); 
cmat(4,:) =  (1-(state.u(1,:)+state.u(2,:))/kcap); 
cmat(4*N+5,:) =  d*(1-(state.u(1,:)+state.u(2,:))/kcap); 
cmat(4*N+8,:) =  d*(1-(state.u(1,:)+state.u(2,:))/kcap); 
end
