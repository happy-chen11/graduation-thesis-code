clear; clc; close all;
n = 5;                         
h = 1/(n+1);                  
A = gallery('tridiag', n, -1, 2, -1) / h^2; 
b = ones(n, 1) * h^2;      
fprintf('========================================\n');
fprintf(' 有限差分矩阵对角占优性验证\n');
fprintf('========================================\n\n');
fprintf('1. 离散后的线性方程组 Au = b\n');
fprintf('   参数设置：\n');
fprintf('   n = %d（内点个数）\n', n);
fprintf('   h = 1/(n+1) = %.4f\n', h);
fprintf('   热源项 f(x) = 1\n\n');
fprintf('   系数矩阵 A（%d × %d）：\n', n, n);
A_display = full(A);
disp(A_display);
fprintf('\n   右端项 b：\n');
disp(b');
fprintf('\n   离散后的线性方程组：\n');
for i = 1:n
    fprintf('   ');
    for j = 1:n
        if A_display(i,j) ~= 0
            if j > 1 && A_display(i,j) > 0
                fprintf(' + ');
            elseif j > 1 && A_display(i,j) < 0
                fprintf(' - ');
            end
            if abs(A_display(i,j)) ~= 1 || j == i
                fprintf('%.0f u_%d', abs(A_display(i,j)), j);
            else
                fprintf('u_%d', j);
            end
        end
    end
    fprintf(' = %.4f\n', b(i));
end
fprintf('\n');
fprintf('2. 对角占优性验证：\n');
fprintf('行号 | 对角元 | 非对角元和 | 对角占优比 | 是否严格\n');
fprintf('------------------------------------------------\n');
for i = 1:n
    diag_val = abs(A_display(i,i));
    off_diag_sum = sum(abs(A_display(i, [1:i-1, i+1:n])));
    ratio = off_diag_sum / diag_val;
    is_strict = diag_val > off_diag_sum;
    strict_str = '是';
    if ~is_strict
        strict_str = '否';
    end
    fprintf(' %2d  |  %6.2f  |     %6.2f    |   %8.4f  |    %s\n', ...
        i, diag_val, off_diag_sum, ratio, strict_str);
end
fprintf('结论：矩阵为不可约对角占优矩阵（非严格对角占优）\n');
fprintf('       - 边界点（第1行和第%d行）：严格对角占优\n', n);
fprintf('       - 内点（第2至%d行）：弱对角占优（等号成立）\n', n-1);
fprintf('       - 矩阵不可约（三对角结构）\n\n');
fprintf('3. 对角占优程度量化：\n');
fprintf('   引入参数 p ∈ (0,1]，构造参数化矩阵 A(p)：\n');
fprintf('   A(p) = diag(2,2,...,2) + p * (次对角线)\n');
fprintf('   对角占优比 δ = p\n');
fprintf('   - 当 p = 1 时，为标准有限差分矩阵（弱对角占优）\n');
fprintf('   - 当 p < 1 时，对角占优程度增强\n\n');
fprintf('4. 谱半径计算：\n');
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
T_J = -D \ (L + U);
T_GS = -(D + L) \ U;
rho_J = max(abs(eig(full(T_J))));
rho_GS = max(abs(eig(full(T_GS))));
fprintf('   雅可比迭代谱半径 ρ_J = %.6f\n', rho_J);
fprintf('   高斯-赛德尔迭代谱半径 ρ_GS = %.6f\n\n', rho_GS);
rho_J_theory = cos(pi/(n+1));
rho_GS_theory = rho_J_theory^2;
fprintf('理论值验证：\n');
fprintf('   ρ_J 理论值 = cos(π/%d) = %.6f\n', n+1, rho_J_theory);
fprintf('   ρ_GS 理论值 = %.6f\n\n', rho_GS_theory);
fprintf('5. 收敛性分析：\n');
tol = 1e-6;
iter_J = ceil(log(tol) / log(rho_J));
iter_GS = ceil(log(tol) / log(rho_GS));
fprintf('   至 %.0e 容差所需步数：\n', tol);
fprintf('   雅可比迭代：约 %d 步\n', iter_J);
fprintf('   高斯-赛德尔迭代：约 %d 步\n', iter_GS);
fprintf('   收敛速度比（GS/J）= %.2f 倍\n', iter_J/iter_GS);
fprintf('   结论：两种迭代法均收敛，高斯-赛德尔迭代收敛速度更快\n\n');
fprintf('6. 对角占优程度对收敛速度的影响：\n');
p_values = [0.2, 0.4, 0.6, 0.8, 1.0];
fprintf('  p   | 对角占优比 |   ρ_J   |  ρ_GS  | 雅可比步数 | GS步数 | 速度比\n');
fprintf('--------------------------------------------------------------\n');
for p = p_values
    delta = p;
    rho_J_p = p * rho_J_theory;
    rho_GS_p = p^2 * rho_GS_theory;
    iter_J_p = ceil(log(tol) / log(rho_J_p));
    iter_GS_p = ceil(log(tol) / log(rho_GS_p));
    speed_ratio = iter_J_p / iter_GS_p;
    fprintf(' %.2f |   %.4f   | %.4f | %.4f |    %3d    |   %3d   |   %.2f\n', ...
        p, delta, rho_J_p, rho_GS_p, iter_J_p, iter_GS_p, speed_ratio);
end
fprintf('结论：对角占优程度越强（p越小），谱半径越小，收敛越快\n');
fprintf('      高斯-赛德尔迭代收敛速度约为雅可比迭代的2倍\n\n');
fprintf('7. 网格加密影响分析：\n');
n_list = [5, 10, 20, 50];
fprintf('   n   |    h    |   ρ_J   |  ρ_GS  | 雅可比步数 | GS步数\n');
fprintf('--------------------------------------------------------\n');
for n_test = n_list
    h_test = 1/(n_test+1);
    cos_val = cos(pi/(n_test+1));
    rho_J_test = cos_val;
    rho_GS_test = cos_val^2;
    iter_J_test = ceil(log(tol) / log(rho_J_test));
    iter_GS_test = ceil(log(tol) / log(rho_GS_test));
    fprintf('  %3d  | %.5f | %.4f | %.4f |    %3d    |   %3d\n', ...
        n_test, h_test, rho_J_test, rho_GS_test, iter_J_test, iter_GS_test);
end
fprintf('结论：网格加密时，谱半径趋近于1，收敛速度变慢\n');
fprintf('      但高斯-赛德尔迭代的相对优势保持不变\n\n');

fprintf('8. 误差衰减规律验证：\n');
p_selected = [0.2, 0.5, 0.8, 1.0];
tol_iter = 1e-6;
max_iter = 200;
x0 = zeros(n, 1);
err_J_all = cell(length(p_selected), 1);
err_GS_all = cell(length(p_selected), 1);
iter_J_all = zeros(length(p_selected), 1);
iter_GS_all = zeros(length(p_selected), 1);
rho_J_p_all = zeros(length(p_selected), 1);
rho_GS_p_all = zeros(length(p_selected), 1);
for idx = 1:length(p_selected)
    p = p_selected(idx);
    rho_J_p_all(idx) = p * rho_J_theory;
    rho_GS_p_all(idx) = p^2 * rho_GS_theory;
    [~, iter_J_all(idx), err_J] = jacobi_method_p611(A, b, x0, p, tol_iter, max_iter);
    err_J_all{idx} = err_J;
    [~, iter_GS_all(idx), err_GS] = gauss_seidel_method_p611(A, b, x0, p, tol_iter, max_iter);
    err_GS_all{idx} = err_GS;
    fprintf('   p=%.1f (δ=%.1f): 雅可比 %d 步, 高斯-赛德尔 %d 步\n', ...
        p, p, iter_J_all(idx), iter_GS_all(idx));
end
fprintf('\n');
figure('Position', [100, 100, 1200, 400]);
subplot(1, 3, 1);
p_plot = linspace(0.1, 1.0, 100);
rho_J_plot = p_plot * rho_J_theory;
rho_GS_plot = p_plot.^2 * rho_GS_theory;
plot(p_plot, rho_J_plot, 'b-', 'LineWidth', 2);
hold on;
plot(p_plot, rho_GS_plot, 'r--', 'LineWidth', 2);
xlabel('对角占优比 δ (或参数 p)', 'FontSize', 11);
ylabel('谱半径 ρ(T)', 'FontSize', 11);
title('(a) 谱半径与对角占优程度的关系', 'FontSize', 12, 'FontWeight', 'bold');
legend('雅可比迭代', '高斯-赛德尔迭代', 'Location', 'northwest', 'FontSize', 9);
grid on;
xlim([0, 1]);
ylim([0, 1]);
subplot(1, 3, 2);
iter_J_plot = ceil(log(tol) ./ log(rho_J_plot));
iter_GS_plot = ceil(log(tol) ./ log(rho_GS_plot));
semilogy(p_plot, iter_J_plot, 'b-', 'LineWidth', 2);
hold on;
semilogy(p_plot, iter_GS_plot, 'r--', 'LineWidth', 2);
xlabel('对角占优比 δ (或参数 p)', 'FontSize', 11);
ylabel('收敛步数（至 1e-6）', 'FontSize', 11);
title('(b) 收敛步数与对角占优程度的关系', 'FontSize', 12, 'FontWeight', 'bold');
legend('雅可比迭代', '高斯-赛德尔迭代', 'Location', 'northeast', 'FontSize', 9);
grid on;
xlim([0, 1]);
subplot(1, 3, 3);
speed_ratio_plot = -log(rho_GS_plot) ./ (-log(rho_J_plot));
plot(p_plot, speed_ratio_plot, 'g-', 'LineWidth', 2);
xlabel('对角占优比 δ (或参数 p)', 'FontSize', 11);
ylabel('收敛速度比 (GS/J)', 'FontSize', 11);
title('(c) 高斯-赛德尔相对于雅可比的速度比', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, 1]);
ylim([1.5, 2.5]);
line([0, 1], [2, 2], 'Color', 'k', 'LineStyle', '--');
text(0.5, 2.03, '理论值 2.0', 'FontSize', 9, 'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
axes('Position', [0, 0.96, 1, 0.04], 'Visible', 'off');
text(0.5, 0.5, '图1：对角占优程度对雅可比和高斯-赛德尔迭代收敛速度的影响分析', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.6]);
figure('Position', [100, 100, 800, 500]);
n_plot = 2:100;
rho_J_n = cos(pi ./ (n_plot + 1));
rho_GS_n = rho_J_n.^2;
plot(n_plot, rho_J_n, 'b-', 'LineWidth', 2);
hold on;
plot(n_plot, rho_GS_n, 'r--', 'LineWidth', 2);
xlabel('内点个数 n', 'FontSize', 12);
ylabel('谱半径 ρ(T)', 'FontSize', 12);
title('图2：谱半径随网格加密的变化（标准有限差分矩阵）', 'FontSize', 12, 'FontWeight', 'bold');
legend('雅可比迭代', '高斯-赛德尔迭代', 'Location', 'southeast', 'FontSize', 10);
grid on;
xlim([2, 100]);
ylim([0.8, 1]);
figure('Position', [100, 100, 1200, 800]);
for idx = 1:length(p_selected)
    p = p_selected(idx);
    subplot(2, 2, idx);
    semilogy(1:iter_J_all(idx), err_J_all{idx}, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    hold on;
    semilogy(1:iter_GS_all(idx), err_GS_all{idx}, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('迭代步数 k', 'FontSize', 10);
    ylabel('误差 ||e^{(k)}||_{\infty}', 'FontSize', 10);
    title(sprintf('(a%d) p = %.1f (δ = %.1f)', idx, p, p), 'FontSize', 11, 'FontWeight', 'bold');
    legend('雅可比迭代', '高斯-赛德尔迭代', 'Location', 'northeast', 'FontSize', 8);
    grid on;
end
axes('Position', [0, 0.96, 1, 0.04], 'Visible', 'off');
text(0.5, 0.5, '图3：不同对角占优程度下的误差衰减曲线', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.6]);
figure('Position', [100, 100, 1200, 800]);
for idx = 1:length(p_selected)
    p = p_selected(idx);
    subplot(2, 2, idx);
    err_J = err_J_all{idx};
    err_GS = err_GS_all{idx}; 
    if length(err_J) >= 2
        rate_J = err_J(2:end) ./ err_J(1:end-1);
        rate_GS = err_GS(2:end) ./ err_GS(1:end-1);
        iter_J_plot = 2:length(err_J);
        iter_GS_plot = 2:length(err_GS);     
        h1 = plot(iter_J_plot, rate_J, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
        hold on;
        h2 = plot(iter_GS_plot, rate_GS, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');  
        x_max = max(length(err_J), length(err_GS));
        line([2, x_max], [rho_J_p_all(idx), rho_J_p_all(idx)], ...
            'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
        line([2, x_max], [rho_GS_p_all(idx), rho_GS_p_all(idx)], ...
            'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);   
        y_min = min(min(rate_J), min(rate_GS)) - 0.1;
        y_max = max(max(rate_J), max(rate_GS)) + 0.15;
        ylim([max(0, y_min), y_max]); 
        text(x_max - 1, rho_J_p_all(idx) + 0.03, sprintf('ρ_J = %.4f', rho_J_p_all(idx)), ...
            'Color', 'b', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        text(x_max - 1, rho_GS_p_all(idx) - 0.04, sprintf('ρ_{GS} = %.4f', rho_GS_p_all(idx)), ...
            'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        if ~isempty(rate_J)
            text(iter_J_plot(end) - 0.5, rate_J(end) + 0.04, sprintf('实际 ≈ %.4f', rate_J(end)), ...
                'Color', 'b', 'FontSize', 7, 'HorizontalAlignment', 'right');
        end
        if ~isempty(rate_GS)
            text(iter_GS_plot(end) - 0.5, rate_GS(end) - 0.05, sprintf('实际 ≈ %.4f', rate_GS(end)), ...
                'Color', 'r', 'FontSize', 7, 'HorizontalAlignment', 'right');
        end 
        xlabel('迭代步数 k', 'FontSize', 10);
        ylabel('误差衰减率 ||e^{(k)}||/||e^{(k-1)}||', 'FontSize', 10);
        title(sprintf('(b%d) p = %.1f (δ = %.1f)', idx, p, p), 'FontSize', 11, 'FontWeight', 'bold'); 
        legend([h1, h2], {'雅可比衰减率', '高斯-赛德尔衰减率'}, 'Location', 'best', 'FontSize', 8);
        grid on;
        xlim([2, x_max]);
    end
end
axes('Position', [0, 0.96, 1, 0.04], 'Visible', 'off');
text(0.5, 0.5, '图4：误差衰减率与谱半径对比（验证误差衰减率趋近于谱半径）', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.6]);
figure('Position', [100, 100, 1000, 500]);
p_demo = 0.5;
[~, iter_J_demo, err_J_demo] = jacobi_method_p611(A, b, x0, p_demo, 1e-8, 200);
[~, iter_GS_demo, err_GS_demo] = gauss_seidel_method_p611(A, b, x0, p_demo, 1e-8, 200);
rho_J_demo = p_demo * rho_J_theory;
rho_GS_demo = p_demo^2 * rho_GS_theory;
subplot(1, 2, 1);
semilogy(1:iter_J_demo, err_J_demo, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on;
semilogy(1:iter_GS_demo, err_GS_demo, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('迭代步数 k', 'FontSize', 11);
ylabel('误差 ||e^{(k)}||_{\infty}', 'FontSize', 11);
title(sprintf('误差衰减曲线 (p = %.1f, δ = %.1f)', p_demo, p_demo), 'FontSize', 12, 'FontWeight', 'bold');
legend('雅可比迭代', '高斯-赛德尔迭代', 'Location', 'northeast', 'FontSize', 9);
grid on;
subplot(1, 2, 2);
if length(err_J_demo) >= 2 && length(err_GS_demo) >= 2
    rate_J = err_J_demo(2:end) ./ err_J_demo(1:end-1);
    rate_GS = err_GS_demo(2:end) ./ err_GS_demo(1:end-1); 
    h1 = plot(2:iter_J_demo, rate_J, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    hold on;
    h2 = plot(2:iter_GS_demo, rate_GS, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r'); 
    x_max = max(iter_J_demo, iter_GS_demo);
    line([2, x_max], [rho_J_demo, rho_J_demo], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
    line([2, x_max], [rho_GS_demo, rho_GS_demo], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    y_min = min(min(rate_J), min(rate_GS)) - 0.1;
    y_max = max(max(rate_J), max(rate_GS)) + 0.15;
    ylim([max(0, y_min), y_max]);
    text(x_max - 1, rho_J_demo + 0.035, sprintf('ρ_J = %.4f', rho_J_demo), ...
        'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(x_max - 1, rho_GS_demo - 0.045, sprintf('ρ_{GS} = %.4f', rho_GS_demo), ...
        'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(iter_J_demo - 0.8, rate_J(end) + 0.045, sprintf('实际 ≈ %.4f', rate_J(end)), ...
        'Color', 'b', 'FontSize', 8, 'HorizontalAlignment', 'right');
    text(iter_GS_demo - 0.8, rate_GS(end) - 0.055, sprintf('实际 ≈ %.4f', rate_GS(end)), ...
        'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'right');
    xlabel('迭代步数 k', 'FontSize', 11);
    ylabel('误差衰减率 ||e^{(k)}||/||e^{(k-1)}||', 'FontSize', 11);
    title(sprintf('误差衰减率与谱半径对比 (p = %.1f)', p_demo), 'FontSize', 12, 'FontWeight', 'bold');
    legend([h1, h2], {'雅可比衰减率', '高斯-赛德尔衰减率'}, 'Location', 'best', 'FontSize', 9);
    grid on;
    xlim([2, x_max]);
end
axes('Position', [0, 0.96, 1, 0.04], 'Visible', 'off');
text(0.5, 0.5, '图5：误差衰减率与谱半径对比（p = 0.5）- 验证误差衰减率趋近于谱半径', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.6]);
fprintf('\n程序运行完成！\n');