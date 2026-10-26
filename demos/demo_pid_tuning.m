% demo_pid_tuning.m
% P / PI / PID 三种控制律的整定对比：三阶被控对象闭环阶跃响应
% 依赖: Octave control 包（或 MATLAB Control System Toolbox）
% 运行: octave-cli demo_pid_tuning.m --no-gui -q
% 产出: pid_tuning.png + 终端打印三控制器指标对比

%% ---- 环境 ----
if exist('OCTAVE_VERSION', 'builtin') ~= 0
  pkg("load", "control");
end

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

function info = step_metrics(t, y, tol_percent)
  yss = y(end);
  info.overshoot = (max(y) - yss) / abs(yss) * 100;
  tol = yss * tol_percent / 100;
  idx = find(abs(y - yss) > tol, 1, 'last');
  if isempty(idx)
    info.settling_time = t(1);
  else
    info.settling_time = t(idx);
  end
  idx_hi = find(y >= 0.9 * yss, 1, 'first');
  idx_lo = find(y >= 0.1 * yss, 1, 'first');
  info.rise_time = t(idx_hi) - t(idx_lo);
  info.steady_state = yss;
end

%% ---- 1. 被控对象: G(s) = 1/(s^3 + 6s^2 + 11s + 6), 极点 -1,-2,-3 ----
G = tf(1, [1 6 11 6]);
fprintf('被控对象 G(s) = 1/(s^3+6s^2+11s+6)\n');

%% ---- 2. 三种控制律 ----
% 纯 P: 有稳态误差
C_P   = tf(2, 1);
% PI: 消除稳态误差, 但相位滞后加剧超调
C_PI  = tf([2 1.5], [1 0]);
% PID: 全面折中
C_PID = tf([4 10 8], [1 0]);

controllers = {C_P, C_PI, C_PID};
names = {'P: Kp=2', 'PI: Kp=2, Ki=1.5', 'PID: Kp=10, Ki=8, Kd=4'};
colors = {'r-', 'g-', 'b-'};

%% ---- 3. 闭环阶跃对比 ----
t = linspace(0, 12, 3000);
h = sim_fig(); hold on;
fprintf('\n%-22s %8s %10s %10s %12s\n', '控制器', '超调%%', '调节时间s', '上升时间s', '稳态误差%%');
for k = 1:3
  T = feedback(controllers{k} * G, 1);
  [y, tt] = step(T, t);
  m = step_metrics(tt, y, 2);
  plot(tt, y, colors{k}, 'linewidth', 1.5);
  fprintf('%-22s %8.2f %10.3f %10.3f %12.3f\n', ...
          names{k}, m.overshoot, m.settling_time, m.rise_time, abs(1 - m.steady_state) * 100);
endfor
plot(t, ones(size(t)), 'k--');
grid on; xlabel('Time (s)'); ylabel('Amplitude');
title('Closed-loop Step: P vs PI vs PID');
legend('P', 'PI', 'PID', 'setpoint');
sim_export(h, 'pid_tuning.png');

fprintf('\n结论: P 有稳态误差; PI 消除稳态误差但超调增大; PID 兼顾速度与超调。\n');
