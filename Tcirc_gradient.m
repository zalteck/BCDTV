function [dfhT, dfvT] = Tcirc_gradient(f)
% Calculate the transpose of the circular gradient of f
    [n,m] = size(f);

    dfhT = [ f(:,1) - f(:,m), f(:,2:m) - f(:,1:m-1) ]; % Horizontal 
    dfvT = [ f(1,:) - f(n,:); f(2:n,:) - f(1:n-1,:) ]; % Vertical

end