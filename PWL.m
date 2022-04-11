function val=PWL(t)
% piecewise linear signal
%----------------------------------------------------------
% Definitions:
[t1, t2] = deal(0, 0.03);
[A1, A2] = deal(0, 1);

    if t<t2
        val = A1;  
    else  
        val = A2;   %<----Fill this
    end

end 
