function [ sum ] = sum_sig_ext( j,m,w,v )
%sum_sig_ext Compute the finite sum for the extinction
%   j is the finite upper bond for the sum

global a b;
sum = 0; % Accumulator variable
for i = 1:numel(j)
    sum = sum + (2*i+1)*(real(a(i,m,w,v)+b(i,m,w,v)));
end

end

