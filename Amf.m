function y = Amf(A, x, direction)
%
% Evaluates the following based on direction
%
%  direction = 0 : return size(A)
%  direction = 1 : returns A*x
%  direction = 2 : returns A'*x
%

if direction == 0
    y = size(A);
elseif direction == 1
        y = A*x;
elseif direction == 2
        y = A'*x;
else
    warning('Amf: No such direction');
end

    