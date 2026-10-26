% demo_nsga2_zdt1.m
% 【算法库】NSGA-II 多目标优化: 标准 ZDT1 基准, Pareto 前沿 vs 理论曲线
% 依赖: 本技能 scripts/ 目录(nsga2.m, sim_fig.m, sim_export.m)
% 运行: octave-cli demo_nsga2_zdt1.m --no-gui -q
% 产出: nsga2_zdt1_front.png + 终端打印前沿质量指标

1;  % 脚本标记
% (自包含: 出图工具内联, nsga2.m 同目录)
function h = sim_fig()
  if exist('OCTAVE_VERSION', 'builtin') ~= 0
    h = figure('visible', 'off');
    graphics_toolkit(h, 'gnuplot');
  else
    h = figure('Visible', 'off');
  end
end

function sim_export(h, outfile)
  print(h, outfile, '-dpng', '-r150');
  fprintf('图已导出: %s (%d bytes)\n', outfile, stat(outfile).size);
end

%% ---- 1. ZDT1 基准(9 变量) ----
% f1 = x1;  g = 1 + 9/8*sum(x(2:9));  f2 = g*(1 - sqrt(f1/g))
% 理论 Pareto 前沿(g=1): f2 = 1 - sqrt(f1), f1 ∈ [0,1]
obj = @(x) [x(1);
            (1 + 9/8*sum(x(2:9))) * (1 - sqrt(x(1) / (1 + 9/8*sum(x(2:9)))))];

fprintf('NSGA-II 求解 ZDT1(9 变量, 50 个体, 150 代)...\n');
rand('state', 5);
[px, pf] = nsga2(obj, 9, zeros(9,1), ones(9,1), ...
                 struct('pop', 50, 'gens', 150));

%% ---- 2. 前沿质量指标 ----
f1 = pf(:,1);  f2_theo = 1 - sqrt(f1);
err_med = median(abs(pf(:,2) - f2_theo));
span = max(f1) - min(f1);
g_vals = 1 + 9/8 * sum(px(:,2:9), 2);
fprintf('\n---- 前沿质量 ----\n');
fprintf('前沿点数     = %d\n', size(px, 1));
fprintf('f1 覆盖范围  = [%.3f, %.3f] (理论 [0,1])\n', min(f1), max(f1));
fprintf('贴合理论曲线 = 中位偏差 %.4f\n', err_med);
fprintf('g 值收敛     = 最大 %.4f (理论最优 1)\n', max(g_vals));

%% ---- 3. 出图: 前沿 vs 理论 ----
h = sim_fig();
f1s = linspace(0, 1, 400);
plot(f1s, 1 - sqrt(f1s), 'r-', 'linewidth', 1.5); hold on;
plot(f1, pf(:,2), 'bo', 'markersize', 4);
grid on;
xlabel('f1'); ylabel('f2');
title('NSGA-II on ZDT1: Pareto Front vs Theory');
legend('theory: 1 - sqrt(f1)', 'NSGA-II', 'location', 'southeast');
sim_export(h, 'nsga2_zdt1_front.png');

fprintf('\n结论: 前沿应铺满 [0,1] 且贴合理论曲线; 若塌缩为单点, 多样性失效。\n');
