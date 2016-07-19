function y = powercone(x, mu, n)

% POWERCONE Returns y such that x.^n / mu.^(n-1) <= y
%
% Currently only accepts scalar, positive integer n. x, mu should be > 0 in
% general.
%
% Helper function to compute exponential cone
%
% Modified By:
% 1. Joel

% if(~isscalar(n) || floor(n) ~= n || n < 1)
%     error('rome_var:powercone:NeedPosIntScalarPower', 'Power must be a positive integer');
% end

if(any(floor(n(:)) ~= n(:)) || any(n(:) < 1))
    error('rome_var:powercone:NeedPosIntScalarPower', 'Power must be a positive integer');
end

% added 27 Mar 2010
n_orig = n;

% degenerate case
if(n == 1)
    y = x;
    return;
end

% size check
if(isscalar(x))
    out_sz = size(mu);
else
    % should go here if mu is a scalar. 
    % If not, or if size is inconsistent, will have an error later
    out_sz = size(x);
end

% make output variable
y = rome_model_var(out_sz);  % do we need this NNOC?
% y = rome_model_var(out_sz, 'Cone', rome_constants.NNOC);  % do we need this NNOC?
rhs = y;

% check n
if(isscalar(n))
    % scalar case:
    while(n > 1)
        % even case
        if(mod(n,2) == 0)
            % make new variable and bound
            rhs = hypoquad(mu, rhs);
            % odd case
        else
            % change lhs / rhs
            rhs = hypoquad(x, rhs);
        end

        % update n
        n = ceil(n / 2);
    end
else
    % vector case:
    while(any(n > 1))
        n_big = (n > 1);
        n_even = find((mod(n, 2) == 0) & n_big);
        n_odd  = find((mod(n, 2) == 1) & n_big);

        if(~isempty(n_even))
            rhs(n_even) = hypoquad(mu(n_even), rhs(n_even));
        end
        if(~isempty(n_odd))
            rhs(n_odd)  = hypoquad(x(n_odd),  rhs(n_odd));
        end

        % update n
        n = ceil(n / 2);
    end
end

% final 
rome_constraint(x <= rhs);

if(mod(n_orig, 2) == 0)
    rome_constraint(x >= -rhs); % Added 27 Mar 2010
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
