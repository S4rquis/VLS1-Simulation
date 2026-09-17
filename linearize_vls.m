function [A, B] = linearize_vls(t0, x0, u0, param, eps)

    if nargin < 5
        eps = 1e-5;
    end

    nx_full = length(x0);
    nu      = length(u0);
    idx_red = [6, 11, 8, 5, 12, 9, 10, 7, 4];  % [w,q,theta, v,r,psi, p,phi,u]

    % Matriz A reduzida
    A = zeros(9, 9);
    for j = 1:nx_full
        dx_p = x0;  dx_m = x0;
        dx_p(j) = x0(j) + eps;
        dx_m(j) = x0(j) - eps;
        fp = vls6dof_rhs(t0, dx_p, u0, param);
        fm = vls6dof_rhs(t0, dx_m, u0, param);
        dfdxj = (fp - fm) / (2*eps);

        if any(idx_red == j)
            col = find(idx_red == j);
            A(:, col) = dfdxj(idx_red);
        end

    end

    % Matriz B reduzida
    B = zeros(9, nu);
    for j = 1:nu
        du_p = u0;  du_m = u0;
        du_p(j) = u0(j) + eps;
        du_m(j) = u0(j) - eps;
        fp = vls6dof_rhs(t0, x0, du_p, param);
        fm = vls6dof_rhs(t0, x0, du_m, param);
        B(:,j) = (fp(idx_red) - fm(idx_red)) / (2*eps);
    end

end