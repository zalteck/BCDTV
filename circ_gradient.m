function [dfh, dfv] = circ_gradient(f)

    [n,m] = size(f);

    dfh = [f(:,1:m-1) - f(:,2:m), f(:,m) - f(:,1)]; % Horizontal 
    dfv = [f(1:n-1,:) - f(2:n,:); f(n,:) - f(1,:)]; % Vertical

end