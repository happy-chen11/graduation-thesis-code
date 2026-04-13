function [x, iter, err_history] = jacobi_method_p621(A, b, x0, p, tol, max_iter)
    n = length(A);
    A_p = A;
    for i = 1:n
        for j = 1:n
            if i ~= j && abs(i-j) == 1
                A_p(i,j) = A(i,j) * p;
            end
        end
    end
    b_p = b;
    D = diag(diag(A_p));
    L = tril(A_p, -1);
    U = triu(A_p, 1);
    T_J = -D \ (L + U);
    c_J = D \ b_p;
    x = x0;
    err_history = zeros(max_iter, 1);
    for iter = 1:max_iter
        x_new = T_J * x + c_J;
        err = norm(x_new - x, inf);
        err_history(iter) = err;
        if err < tol
            x = x_new;
            break;
        end
        x = x_new;
    end
    err_history = err_history(1:iter);
end
