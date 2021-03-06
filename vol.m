function xr = vol(n1,n2,val)
% Adds the stamp of an independent voltage source with a value
% of "val" (Volts) connected between nodes n1 and n2 to the 
% matrices in circuit representation.
%
%                   val 
%                  /  \
%      n1 O-------(+  -)--------O n2    where Vsrc= val (volts)
%                  \  /
%             Isrc ---->
%---------------------------------------------------------------

% The body of the function will go here!
global G C b   %define global variables

d = size(C,1); % current size of the MNA
xr = d+1;      % new (extra)row/column
b(xr) = 0;      % add new row
G(xr,xr) = 0;   % add new row/column
C(xr,xr) = 0;
 
if (n1 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
    b(xr) = val;
end
if (n2 ~= 0)
    G(n2,xr) = -1;
    G(xr,n2) = -1;
end

end
