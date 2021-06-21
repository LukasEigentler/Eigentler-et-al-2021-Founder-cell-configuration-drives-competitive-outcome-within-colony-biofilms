function uinit = u0fun_both_non_uniform(location)
%U0FUN_BOTH_NON_UNIFORM defines initial condition using microcolonies
% 
% U0FUN_BOTH_NON_UNIFORM(location) defines the model's initial condition on the
% spatial mesh based on the locations of the mesh's nodes given by the
% location vector

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 21/06/2021

global kcap Cxloc_in_B1 Cyloc_in_B1 Cxloc_in_B2 Cyloc_in_B2
M = length(location.x);
uinit = zeros(2,M);
uinit(1,:) = 0;
uinit(2,:) = 0;
for cc = 1:length(Cxloc_in_B1) % loop through all B1 patches
    d = sqrt((location.x-Cxloc_in_B1(cc)).^2 + (location.y-Cyloc_in_B1(cc)).^2); % find closes mesh node
    [~,mind_ind] = min(d);
    uinit(1,mind_ind) = kcap; % assign B1 = kcap at the node
end
for cc = 1:length(Cxloc_in_B2) % loop through all B2 patches
    d = sqrt((location.x-Cxloc_in_B2(cc)).^2 + (location.y-Cyloc_in_B2(cc)).^2); % find closes mesh node
    [~,mind_ind] = min(d);
    uinit(2,mind_ind) = kcap; % assign B2 = kcap at the node
end
end
