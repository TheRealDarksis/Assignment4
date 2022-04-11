function vcvs(nd1,nd2,ni1,ni2,val)
% Adds the stamp of a dependent voltage-controlled 
% voltage-source(VCVS)to the matrices in circuit 
% representation.
%
%   ni1 O-------(+)        |----------o nd1
%                          |
%                         /+\
%                      | /   \    Vnd1-Vnd2 = val*(Vni1-Vni2)
%                Ivcvs | \   /
%                      V  \-/ 
%                          |
%   ni2 O-------(-)        |----------o nd2
%
%  (1) "nd1 & nd2" are the nodes across the dependent
%                  voltage source.
%  (2) "ni1 & ni2" are the nodes corresponding to the 
%                  controller voltage
%
%   nd1: (+) node   \
%   nd2: (-) node   |----->  Vnd1-Vnd2 = val*(Vni1-Vni2)
%   ni1: (+) node   |
%   ni2: (-) node   /
%---------------------------------------------------------------

% The body of the function will go here!




global G C b  

d = size(G,1); % current size of the MNA
xr = d+1;      % new (extera)row/column
b(xr) = 0;      % add new row
G(xr,xr) = 0;   % add new row/column
C(xr,xr) = 0;

if (nd1 ~= 0)
    G(nd1,xr) =  1;
    G(xr,nd1) =  1;
end

if (nd2 ~= 0)
    G(nd2,xr) =  1;
    G(xr,nd2) =  1;
end
if (ni1 ~= 0)
    G(xr,ni1) =  G(xr,ni1) - val; 
end
if (ni2 ~= 0)
    G(xr,ni2) =  G(xr,ni2) + val;    
end

end
