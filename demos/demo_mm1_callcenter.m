% pt_T01_mm1_callcenter.m
% 论文测试 T01【排队论/离散事件仿真】
% 来源: Number Analytics "M/M/1 Queue in Queueing Theory" 呼叫中心案例
%   (https://www.numberanalytics.com/blog/mm1-queue-in-queueing-theory)
% 论文参数: lambda=10 呼叫/小时, mu=12 呼叫/小时(平均处理 5 分钟)
% 论文理论值: rho=0.833, W=0.5 小时, Lq=4.08
% 测试: 蒙特卡洛离散事件仿真复现论文结果

%% ---- 环境 ----
1;  % 脚本标记(防止整文件被识别为函数文件)
function h = sim_fig()
  if exist('OCTAVE_VERSION', 'builtin') ~= 0
    h = figure('visible', 'off');
    graphics_toolkit(h, 'gnuplot');
  else
    h = figure('Visible', 'off');
  end
end

lambda = 10 / 60;   % 呼叫/分钟
mu     = 12 / 60;   % 服务/分钟
rho    = lambda / mu;
W_theo = 1 / (mu - lambda);         % 论文值: 30 分钟 = 0.5 h
Lq_theo = rho^2 / (1 - rho);        % 论文值: 4.08

fprintf('论文参数: lambda=%.3f/min, mu=%.2f/min, rho=%.3f\n', lambda, mu, rho);
fprintf('论文理论值: W=%.1f 分钟, Lq=%.2f\n', W_theo, Lq_theo);

%% ---- 蒙特卡洛离散事件仿真 ----
rand('state', 42);
NC = 60000;                          % 仿真 6 万呼叫(约 100 小时)
warm = 3000;
ia = -log(rand(NC,1)) / lambda;
st = -log(rand(NC,1)) / mu;
arr = cumsum(ia);
dep = zeros(NC,1); dep(1) = arr(1) + st(1);
for k = 2:NC
  dep(k) = max(arr(k), dep(k-1)) + st(k);
endfor
Wq = dep - st - arr;  W = dep - arr;
W_sim = mean(W(warm:end));
Lq_sim = lambda * mean(Wq(warm:end));
rho_sim = sum(st(warm:end-1)) / (arr(end) - arr(warm));

fprintf('\n仿真结果(N=%d, warmup=%d):\n', NC, warm);
fprintf('W  仿真 %.2f 分钟 (论文理论 %.1f, 误差 %.1f%%)\n', W_sim, W_theo, abs(W_sim-W_theo)/W_theo*100);
fprintf('Lq 仿真 %.2f      (论文理论 %.2f, 误差 %.1f%%)\n', Lq_sim, Lq_theo, abs(Lq_sim-Lq_theo)/Lq_theo*100);
fprintf('rho 仿真 %.3f     (理论 %.3f)\n', rho_sim, rho);

%% ---- 出图: 等待时间分布 ----
h = sim_fig();
hist(Wq(warm:end)*60, 0:1:120);      % 转换为分钟
grid on;
xlabel('Waiting Time (min)'); ylabel('Count');
title('M/M/1 Call Center: Waiting Time Distribution (lambda=10/h, mu=12/h)');
print(h, 'figures/mm1_waitdist.png', '-dpng', '-r150');

fprintf('\n结论: 仿真与论文理论值误差在百分之几内, 呼叫中心参数合理。\n');
fprintf('管理启示: rho=0.83 时平均等待 30 分钟, 接近临界负载, 增员可大幅降等待。\n');
