function cur(n1,n2,val)
% Adds the stamp of an independent current source with a 
% value of "val" (Amperes) connected between nodes n1 and 
% n2 to the the source vector b in circuit representation.
%
%                   val
%                  /  \
%      n1 O-------(--->)-------O n2   where J=val (Amperes)
%                  \  /
%
%      n1: The node at the tail of the current arrow!
%      n2:  "    "   "   " head  "    "   "      "  !
%     val: The value of the current source (Amp)
%----------------------------------------------------------
global F   %define global variable

if (n1 ~= 0)
    F(n1,1) = -val;
end

if (n2 ~= 0)
    F(n2,1) = val;
end

if (n1 ~= 0) && (n2 ~= 0)
    F(n1,n2) = F(n1,n2) - val;
    F(n2,n1) = F(n2,n1) - val;
end

end %func

