clear; clc;
n = 3;
A = gallery('tridiag', n, -1, 2, -1);
A_full = full(A); 
fprintf('原矩阵 A（弱对角占优）：\n');
disp(A_full);
fprintf('\n原矩阵对角占优性验证：\n');
fprintf('行号 | 对角元 | 非对角元和 | 对角占优比 | 是否严格\n');
fprintf('------------------------------------------------\n');
for i = 1:n
    diag_val = abs(A_full(i,i));
    off_diag_sum = sum(abs(A_full(i, [1:i-1, i+1:n])));
    ratio = off_diag_sum / diag_val;
    is_strict = diag_val > off_diag_sum;
    if is_strict
        strict_str = '是';
    else
        strict_str = '否';
    end
    fprintf(' %2d  |   %6.2f  |     %6.2f    |   %8.4f  |    %s\n', ...
        i, diag_val, off_diag_sum, ratio, strict_str);
end
fprintf('结论：内点行（第2行）取等号，为弱对角占优矩阵\n\n');
fprintf('========== 方法一：小扰动法 ==========\n\n');
epsilon = 0.1;
A_pert = A_full + epsilon * eye(n);
fprintf('取 ε = %.1f，扰动后矩阵 A_ε：\n', epsilon);
disp(A_pert);
fprintf('\n扰动后矩阵严格对角占优性验证：\n');
fprintf('行号 | 对角元 | 非对角元和 | 对角占优比 | 是否严格\n');
fprintf('------------------------------------------------\n');
for i = 1:n
    diag_val = abs(A_pert(i,i));
    off_diag_sum = sum(abs(A_pert(i, [1:i-1, i+1:n])));
    ratio = off_diag_sum / diag_val;
    is_strict = diag_val > off_diag_sum;
    if is_strict
        strict_str = '是';
    else
        strict_str = '否';
    end
    fprintf(' %2d  |   %6.2f  |     %6.2f    |   %8.4f  |    %s\n', ...
        i, diag_val, off_diag_sum, ratio, strict_str);
end
fprintf('结论：小扰动法成功将弱对角占优矩阵转化为严格对角占优矩阵\n\n');
fprintf('========== 方法二：正对角相似变换法 ==========\n\n');
fprintf('尝试几何级数形式 D = diag(1, r, 1)\n\n');
r_values = [0.5, 0.8, 1.0, 1.2, 1.5];
for r = r_values
    D = diag([1, r, 1]);
    B = D \ A_full * D;
    fprintf('r = %.1f 时，变换后矩阵 B：\n', r);
    disp(B);
    fprintf('严格对角占优性验证：\n');
    all_strict = true;
    for i = 1:n
        diag_val = abs(B(i,i));
        off_diag_sum = sum(abs(B(i, [1:i-1, i+1:n])));
        if diag_val > off_diag_sum
            fprintf('  第%d行：%.2f > %.2f ?\n', i, diag_val, off_diag_sum);
        else
            fprintf('  第%d行：%.2f ≤ %.2f ?\n', i, diag_val, off_diag_sum);
            all_strict = false;
        end
    end
    if all_strict
        fprintf('  ? 成功转化为严格对角占优矩阵！\n\n');
    else
        fprintf('  ? 未能转化为严格对角占优矩阵\n\n');
    end
end
fprintf('结论：对称三对角矩阵无法通过正对角相似变换转化为严格对角占优矩阵。\n');
fprintf('      原因：第2行条件 2 > r + 1/r 无解（AM-GM不等式）。\n\n');
fprintf('========== 不同扰动量的影响 ==========\n\n');
epsilon_list = [1e-1, 1e-2, 1e-4, 1e-6, 1e-8];
fprintf('ε取值\t\t最大对角占优比\t是否严格\n');
fprintf('----------------------------------------\n');
for eps_val = epsilon_list
    A_test = A_full + eps_val * eye(n);
    max_ratio = 0;
    for i = 1:n
        diag_val = abs(A_test(i,i));
        off_diag_sum = sum(abs(A_test(i, [1:i-1, i+1:n])));
        ratio = off_diag_sum / diag_val;
        max_ratio = max(max_ratio, ratio);
    end
    is_strict = max_ratio < 1;
    if is_strict
        strict_str = '是';
    else
        strict_str = '否';
    end
    fprintf('%.0e\t\t%.6f\t\t%s\n', eps_val, max_ratio, strict_str);
end
fprintf('\n结论：取 ε = 1e-6 即可满足严格对角占优条件\n');