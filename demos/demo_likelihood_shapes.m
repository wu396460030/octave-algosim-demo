% demo_likelihood_shapes.m
% 【算法库】似然面四形态: 光滑单峰 / 多峰非凸 / 平坦山脊 / 高维简并
% 依赖: 无（单文件自包含, 引擎函数内联）
% 运行: octave-cli demo_likelihood_shapes.m --no-gui -q
% 产出: 4 张图(figures/lik_*.png) + 终端打印形态诊断表
%
% 为什么看这个 / Why this matters:
%   极大似然估计不是"跑出一个数"就结束了。似然面的**形状**决定了这个数能
%   不能信、置信区间怎么算、要不要换求解器:
%     smooth      光滑单峰 — MLE + Hessian 协方差即可, 区间可信
%     multimodal  多峰非凸 — 局部优化只锁定一个峰, 须多起点 + 声明标签
%     ridge       平坦山脊 — 单参数区间无限宽, 只能报参数组合
%     degenerate  高维简并 — Hessian 秩亏, 某些方向根本没有信息
%
% 判据 / Criterion:
%   曲率比 = 观测 Fisher 信息(峰处负 Hessian)最大/最小特征值之比
%     < 30   smooth |  30 ~ 1e4  ridge |  >= 1e4  degenerate
%   峰数 >= 2 且两峰落在互相分离的 95% 置信域内 -> multimodal

1;  % 脚本标记 / script marker

% ---- 出图工具内联(sim_fig / sim_export) ----
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
  s = dir(outfile);
  fprintf('图已导出: %s (%d bytes)\n', outfile, s.bytes);
end

function outfile = figpath(name)
  % 统一输出到同目录 figures/ (不存在则创建)
  outdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
  if ~exist(outdir, 'dir')
    mkdir(outdir);
  end
  outfile = fullfile(outdir, [name '.png']);
end

function t = chi2_thr_local(conf, d)
  if exist('chi2inv', 'file') == 2
    try
      t = chi2inv(conf, d);  return;
    catch
    end
  end
  tab = [0.90 2.706 4.605; 0.95 3.841 5.991; 0.99 6.635 9.210];
  [~, k] = min(abs(tab(:,1) - conf));
  t = tab(k, min(d,2) + 1);
end

function v = feval_ll(f, th)
  v = f(th);
  v = v(1);
  if ~isfinite(v), v = -1e300; end
end

function [lab, n] = label_regions_local(mask)
  % 4-连通域标记(等价 bwlabel, 不依赖 image 包)
  [n2,n1] = size(mask);
  lab = zeros(n2,n1);  n = 0;
  st = zeros(n2*n1, 2);
  for i0 = 1:n2
    for j0 = 1:n1
      if mask(i0,j0) && lab(i0,j0)==0
        n = n + 1;
        top = 1;  st(1,:) = [i0,j0];  lab(i0,j0) = n;
        while top >= 1
          i = st(top,1);  j = st(top,2);  top = top - 1;
          nb = [i-1 j; i+1 j; i j-1; i j+1];
          for q = 1:4
            ii = nb(q,1);  jj = nb(q,2);
            if ii>=1 && ii<=n2 && jj>=1 && jj<=n1 && mask(ii,jj) && lab(ii,jj)==0
              lab(ii,jj) = n;  top = top + 1;  st(top,:) = [ii,jj];
            end
          end
        end
      end
    end
  end
end

function H = num_hess_local(f, th, h)
  % 中心差分 Hessian
  th = th(:);  d = numel(th);
  if numel(h) == 1, h = h*ones(d,1); end
  h = abs(h(:));
  for i = 1:d
    if h(i) < 1e-10, h(i) = max(abs(th(i)),1)*1e-5; end
  end
  H = zeros(d,d);
  f0 = f(th);
  for i = 1:d
    for j = i:d
      ei = zeros(d,1);  ei(i) = h(i);
      ej = zeros(d,1);  ej(j) = h(j);
      if i == j
        H(i,i) = (f(th+ei) - 2*f0 + f(th-ei)) / h(i)^2;
      else
        H(i,j) = (f(th+ei+ej) - f(th+ei-ej) - f(th-ei+ej) + f(th-ei-ej)) / (4*h(i)*h(j));
        H(j,i) = H(i,j);
      end
    end
  end
end

function S = lik_surface(loglik, grid, opts)
  % LIK_SURFACE 网格似然面 + 峰识别 + 顶点精化 + 形态诊断
  if nargin < 3, opts = struct(); end
  if ~isfield(opts, 'conf'), opts.conf = 0.95; end
  d = numel(grid.n);
  thr = chi2_thr_local(opts.conf, d);
  if d > 2
    S = lik_surface_hd_local(loglik, grid, opts);
    return;
  end
  a = linspace(grid.lo(1), grid.hi(1), grid.n(1));
  if d == 2
    b = linspace(grid.lo(2), grid.hi(2), grid.n(2));
    [A,B] = meshgrid(a,b);
    L = zeros(numel(b), numel(a));
    for j = 1:numel(a)
      for i = 1:numel(b)
        L(i,j) = feval_ll(loglik, [A(i,j); B(i,j)]);
      end
    end
    da = a(2)-a(1);  db = b(2)-b(1);
  else
    A = a;  B = 0;
    L = zeros(1, numel(a));
    for j = 1:numel(a), L(j) = feval_ll(loglik, a(j)); end
    da = a(2)-a(1);  db = 0;
  end

  Lmax = max(L(:));
  D = 2*(Lmax - L);
  [ir, ic] = find(L == Lmax, 1);
  if d == 1, th0 = a(ir); else, th0 = [A(ir,ic); B(ir,ic)]; end

  % 局部极大(3x3 邻域严格极大)
  P = -inf(size(L)+2);
  P(2:end-1, 2:end-1) = L;
  ismax = true(size(L));  strict = false(size(L));
  for di = 0:2
    for dj = 0:2
      if di == 1 && dj == 1, continue; end
      nb = P(1+di:di+size(L,1), 1+dj:dj+size(L,2));
      ismax = ismax & (L >= nb);
      strict = strict | (L > nb);
    end
  end
  [pr, pc] = find(ismax & strict);
  if isempty(pr), [~,k] = max(L(:)); [pr,pc] = ind2sub(size(L),k); end

  % 逐峰顶点精化(精化点似然不升则回退网格点)
  peaks = struct('idx',{}, 'th',{}, 'L',{}, 'th_v',{}, 'cov',{}, 'delta',{});
  for q = 1:numel(pr)
    i = pr(q);  j = pc(q);
    p.th = [A(i,j); B(i,j)];  p.L = L(i,j);  p.idx = [i,j];
    p.th_v = p.th;  p.cov = NaN(2,2);  p.delta = 2*(Lmax - L(i,j));
    r0 = max(i-1,1):min(i+1,size(L,1));
    c0 = max(j-1,1):min(j+1,size(L,2));
    if numel(r0)==3 && numel(c0)==3 && d == 2
      Xb = A(r0,c0);  Yb = B(r0,c0);  Lb = L(r0,c0);
      u1 = (Xb(:)-A(i,j))/da;  u2 = (Yb(:)-B(i,j))/db;   % 9x1
      u1 = u1(:);  u2 = u2(:);  yy = Lb(:);
      Mm = [ones(numel(yy),1), u1, u2, u1.^2, u1.*u2, u2.^2];
      pp = Mm \ yy;
      bb = [pp(2)/da; pp(3)/db];
      Aq = [pp(4)/da^2, pp(5)/(2*da*db); pp(5)/(2*da*db), pp(6)/db^2];
      if rcond(Aq) > 1e-14
        tv = p.th + (-0.5*(Aq\bb));
        if feval_ll(loglik, tv) > L(i,j)
          p.th_v = tv;
          Cc = inv(-2*Aq);  p.cov = (Cc+Cc')/2;
        end
      end
    end
    peaks(end+1) = p;
  end

  % Hessian 谱(步长自适应: 逐次减半直到曲率比收敛, 且保证 -H 正定)
  % Rosenbrock 山脊实测: h/4 -> 2498, h/128 -> 2508.000(解析值), 故必须收敛判据
  H = [];
  r_prev = NaN;
  for kk = 0:10
    Hk = num_hess_local(@(t) feval_ll(loglik,t), th0, [da; db]/(4*2^kk));
    evk = sort(eig(-Hk));
    if evk(1) <= 0
      continue;                      % -H 非正定(步长过大), 继续缩小
    end
    H = Hk;
    if numel(evk) >= 2
      rk = evk(end)/evk(1);
    else
      rk = 1;
    end
    if ~isnan(r_prev) && abs(rk - r_prev) / max(abs(rk), 1e-300) < 1e-4
      break;                         % 收敛
    end
    r_prev = rk;
  end
  if isempty(H)
    H = num_hess_local(@(t) feval_ll(loglik,t), th0, [da; db]/4);
  end
  [V,Ev] = eig(-H);
  ev = sort(diag(Ev));
  if numel(ev) >= 2 && ev(1) > 0
    ratio = ev(end)/ev(1);
  else
    ratio = Inf;
  end

  % 可分辨峰数: 95% 置信域连通判据
  in95 = (D <= thr);
  [lab, nlab] = label_regions_local(in95);   %#ok<ASGLU>
  regs = [];
  for q = 1:numel(peaks)
    if peaks(q).L >= Lmax - thr/2
      regs(end+1) = lab(peaks(q).idx(1), peaks(q).idx(2));
    end
  end
  regs = unique(regs(regs > 0));
  n_sig = max(numel(regs), 1);

  if n_sig >= 2
    shape = 'multimodal';
  elseif ratio >= 1e4
    shape = 'degenerate';
  elseif ratio >= 30
    shape = 'ridge';
  else
    shape = 'smooth';
  end

  S.A = A;  S.B = B;  S.L = L;  S.D = D;  S.Lmax = Lmax;  S.d = d;
  S.th_hat = th0(:);
  for q = 1:numel(peaks)
    if peaks(q).idx(1)==ir && peaks(q).idx(2)==ic && all(isfinite(peaks(q).th_v))
      S.th_hat = peaks(q).th_v(:);  break;
    end
  end
  S.cov = pinv(-H);  S.cov = (S.cov+S.cov')/2;
  S.peaks = peaks;  S.npeaks = n_sig;
  S.diag = struct('eig', ev, 'ratio', ratio, 'V', V, 'thr', thr, ...
                  'frac95', sum(in95(:))/numel(L));
  S.shape = shape;
end

function S = lik_surface_hd_local(loglik, grid, opts)
  % d>2: 多起点 MLE + Hessian 谱 + 零空间(不做全网格)
  d = numel(grid.n);
  th0 = (grid.lo(:) + grid.hi(:))/2;
  nstart = 1;
  if isfield(opts, 'nstart'), nstart = opts.nstart; end
  o = optimset('MaxFunEvals', 20000, 'MaxIter', 20000, ...
               'TolX', 1e-10, 'TolFun', 1e-12, 'Display', 'off');
  best = -inf;  thb = th0;
  rand('state', 7);
  for s = 1:nstart
    if s == 1
      x0 = th0;
    else
      x0 = grid.lo(:) + (grid.hi(:)-grid.lo(:)).*rand(d,1);
    end
    [x, fv] = fminsearch(@(z) -feval_ll(loglik, z), x0, o);
    if -fv > best, best = -fv;  thb = x; end
  end
  h0 = (grid.hi(:) - grid.lo(:)) / max(grid.n(:));
  H = num_hess_local(@(t) feval_ll(loglik,t), thb, h0);
  [V,Ev] = eig(-H);
  [ev, ord] = sort(diag(Ev));
  V = V(:, ord);
  tol = max(abs(ev(end)), 1e-300) * 1e-8;
  ni = find(abs(ev) < tol);
  rank_h = d - numel(ni);
  if rank_h < d
    shape = 'degenerate';
    nullsp = V(:, ni);
    ratio = Inf;
  else
    nullsp = zeros(d, 0);
    if ev(1) > 1e-12*max(abs(ev(end)),1e-300)
      ratio = ev(end)/ev(1);
      if ratio >= 1e4, shape = 'degenerate';
      elseif ratio >= 30, shape = 'ridge';
      else, shape = 'smooth'; end
    else
      ratio = Inf;  shape = 'degenerate';
    end
  end
  S = struct();
  S.th_hat = thb(:);  S.Lmax = best;  S.d = d;
  S.cov = pinv(-H);  S.cov = (S.cov+S.cov')/2;
  S.peaks = struct('th_v', thb(:), 'L', best, 'delta', 0);
  S.npeaks = 1;
  S.rank = rank_h;  S.nnull = d - rank_h;  S.nullspace = nullsp;
  S.diag = struct('eig', ev, 'ratio', ratio, 'V', V);
  S.shape = shape;
end

function v = prof_obj(loglik, xf, idx, tfix)
  % 把自由参数 xf(行或列) 拼回完整参数(第 idx 位固定为 tfix), 返回 -loglik
  % 注意: fminsearch 传进来的是行向量, 必须先转列再纵向拼接
  xf = xf(:);
  v = -feval_ll(loglik, [xf(1:idx-1); tfix; xf(idx:end)]);
end

function [Lb, thb] = max_free_local(loglik, th_seed, idx, tfix, lo, hi, o, nstart, npert)
  % 固定第 idx 位为 tfix, 对其余参数做多起点极大化
  d = numel(th_seed);
  free = 1:d;  free(idx) = [];
  if isempty(free)
    thb = th_seed;  thb(idx) = tfix;  Lb = feval_ll(loglik, thb);
    return;
  end
  f = @(xf) prof_obj(loglik, xf, idx, tfix);
  nf = numel(free);
  seeds = zeros(nf*max(nstart,1)+1, nf);
  k = 1;  seeds(k,:) = th_seed(free)';
  rng = (hi(free) - lo(free))';
  for s = 1:max(nstart-1,0)
    for p = 1:npert
      k = k + 1;
      seeds(k,:) = max(lo(free)', min(hi(free)', ...
                   th_seed(free)' + rng.*(2*rand(1,nf)-1)*0.3));
    end
  end
  Lb = -inf;  thb = th_seed;
  for k2 = 1:size(seeds,1)
    [x, fv] = fminsearch(f, seeds(k2,:), o);
    if -fv > Lb
      Lb = -fv;  thb = th_seed;  thb(free) = x';
    end
  end
end

function t = cross_local(tv, Dp, thr, k, side)
  % 线性内插 D=thr 的穿越点
  n = numel(tv);  t = tv(k);
  k2 = k + side;
  if k2 < 1 || k2 > n, return; end
  d1 = Dp(k);  d2 = Dp(k2);
  if d2 <= d1, return; end
  if d1 <= thr && d2 > thr
    t = tv(k) + ((thr-d1)/(d2-d1)) * (tv(k2)-tv(k));
  end
end

function P = lik_profile(loglik, th_hat, idx, lo, hi, npt, opts)
  % 剖面似然 + 可辨识性三分类判决
  if nargin < 7, opts = struct(); end
  if ~isfield(opts, 'nstart'), opts.nstart = 5; end
  if ~isfield(opts, 'conf'), opts.conf = 0.95; end
  thr = chi2_thr_local(opts.conf, 1);    % Wilks: D ~ chi2_1, 95% = 3.841
  lo = lo(:);  hi = hi(:);
  o = optimset('MaxFunEvals', 20000, 'MaxIter', 20000, ...
               'TolX', 1e-10, 'TolFun', 1e-12, 'Display', 'off');
  rand('state', 11);
  [Lmax, th_g] = max_free_local(loglik, th_hat, idx, th_hat(idx), lo, hi, o, 6, 3);
  tv = linspace(lo(idx), hi(idx), npt)';
  Lp = zeros(npt,1);
  for k = 1:npt
    Lp(k) = max_free_local(loglik, th_g, idx, tv(k), lo, hi, o, opts.nstart, 2);
  end
  Dp = 2*(Lmax - Lp);
  in95 = (Dp <= thr);
  if any(in95)
    tlo = cross_local(tv, Dp, thr, find(in95,1,'first'), -1);
    thi = cross_local(tv, Dp, thr, find(in95,1,'last'), 1);
  else
    kmin = find(Dp == min(Dp), 1);
    tlo = tv(kmin);  thi = tv(kmin);
  end
  tole = (hi(idx)-lo(idx))*1e-3;
  open_lo = (tlo <= lo(idx)+tole);
  open_hi = (thi >= hi(idx)-tole);
  flat = all(Dp <= thr);
  if flat
    cls = 'structural';  int95 = [NaN; NaN];
  elseif open_lo || open_hi
    cls = 'practical';   int95 = [tlo; thi];
  else
    cls = 'identifiable'; int95 = [tlo; thi];
  end
  P.t = tv;  P.D = Dp;  P.int95 = int95;  P.class = cls;
  P.thr = thr;  P.flat = flat;  P.open = [open_lo, open_hi];
end

%% ===================== 主流程 =====================

%% ---- 1. 光滑单峰: 线性回归似然面 ----
% 锚点: Fisher (1922) doi:10.1098/rsta.1922.0009
% 验证: 面上 MLE = OLS 解析解; Hessian 反演协方差 = CRLB
fprintf('=== 1. smooth (Fisher 1922) ===\n');
rand('state', 21);  randn('state', 21);
n = 40;
x = linspace(-2, 2, n)';
y = 1.5 + 0.8 * x + 0.5 * randn(n, 1);
X = [ones(n,1), x];
beta_ols = X \ y;
ll1 = @(th) -n/2 * log(sum((y - X*[th(1); th(2)]).^2) / n);

g1.lo = [0.5; 0.2];  g1.hi = [2.5; 1.4];  g1.n = [121; 121];
S1 = lik_surface(ll1, g1);
fprintf('  MLE (surface)  = [%.4f %.4f]\n', S1.th_hat(1), S1.th_hat(2));
fprintf('  OLS (analytic) = [%.4f %.4f]   max diff %.2e\n', ...
        beta_ols(1), beta_ols(2), max(abs(S1.th_hat - beta_ols)));
fprintf('  shape = %s, curvature ratio = %.2f\n', S1.shape, S1.diag.ratio);

%% ---- 2. 多峰非凸: 等权双高斯混合 ----
% 锚点: Richardson & Green (1997) doi:10.1111/1467-9868.00095
% 验证: 两个等高全局峰 = 标签互换对; 严格满足 L(m1,m2) = L(m2,m1)
fprintf('\n=== 2. multimodal (Richardson & Green 1997) ===\n');
rand('state', 23);  randn('state', 23);
y2 = zeros(60, 1);
for i = 1:60
  if rand < 0.5
    y2(i) = -2 + 0.7*randn;
  else
    y2(i) =  2 + 0.7*randn;
  end
end
sg = 0.7;
phi = @(m) exp(-0.5*((y2 - m)/sg).^2) / (sg*sqrt(2*pi));
ll2 = @(th) sum(log(0.5*phi(th(1)) + 0.5*phi(th(2)) + 1e-300));
g2.lo = [-4; -4];  g2.hi = [4; 4];  g2.n = [161; 161];
S2 = lik_surface(ll2, g2);
Lpk = [S2.peaks.L];  [~, ord] = sort(Lpk, 'descend');
fprintf('  global mode A  = (%.3f, %.3f)\n', ...
        S2.peaks(ord(1)).th(1), S2.peaks(ord(1)).th(2));
fprintf('  global mode B  = (%.3f, %.3f)  height diff = %.2e\n', ...
        S2.peaks(ord(2)).th(1), S2.peaks(ord(2)).th(2), ...
        abs(S2.peaks(ord(1)).L - S2.peaks(ord(2)).L));
fprintf('  shape = %s, resolvable peaks = %d\n', S2.shape, S2.npeaks);

%% ---- 3. 平坦山脊: 负 Rosenbrock ----
% 锚点: Rosenbrock (1960) doi:10.1093/comjnl/3.3.175
% 验证: 曲率比 vs 解析 Hessian [[802,-400],[-400,200]] 特征值比 2508
fprintf('\n=== 3. ridge (Rosenbrock 1960) ===\n');
ll3 = @(th) -(100*(th(2) - th(1)^2)^2 + (1 - th(1))^2);
g3.lo = [-2; -1];  g3.hi = [2; 3];  g3.n = [161; 161];
S3 = lik_surface(ll3, g3);
eva = sort(eig([802 -400; -400 200]));
fprintf('  peak = (%.4f, %.4f)  (true optimum (1,1))\n', S3.th_hat(1), S3.th_hat(2));
fprintf('  curv ratio: engine %.1f vs analytic %.1f\n', S3.diag.ratio, eva(2)/eva(1));
fprintf('  shape = %s\n', S3.shape);

%% ---- 4. 高维简并: 3 参数仅 2 个线性组合可辨识 ----
% 锚点: Gutenkunst et al. (2007) doi:10.1371/journal.pcbi.0030189
%       剖面判决方法 Raue et al. (2009) doi:10.1093/bioinformatics/btp358
% 验证: Hessian 秩 2 < 3; 沿零空间方向似然严格平坦
fprintf('\n=== 4. degenerate (Gutenkunst 2007 / Raue 2009) ===\n');
A = [1 0 1; 0 1 1];
A(1,:) = A(1,:) / norm(A(1,:));
A(2,:) = A(2,:) / norm(A(2,:));
th_true = [0.8; -0.3; 0.2];
b = A * th_true;
s = [0.30; 0.10];
ll4 = @(th) -0.5 * sum(((A*th(:) - b) ./ s).^2);
g4.lo = [-3;-3;-3];  g4.hi = [3;3;3];  g4.n = [11;11;11];
S4 = lik_surface(ll4, g4, struct('nstart', 6));
fprintf('  shape = %s, rank = %d / %d, null dim = %d\n', ...
        S4.shape, S4.rank, S4.d, S4.nnull);
if S4.nnull > 0
  v0 = S4.nullspace(:,1);
  fprintf('  |L(theta) - L(theta - 5*null_dir)| = %.2e  (flat)\n', ...
          abs(ll4(S4.th_hat) - ll4(S4.th_hat - 5*v0)));
end
P4 = lik_profile(ll4, S4.th_hat, 3, g4.lo, g4.hi, 31);
fprintf('  profile(theta3) = %s  (structural = 不可辨识)\n', P4.class);

%% ---- 5. 汇总表 ----
fprintf('\n==== Shape summary / 形态汇总 ====\n');
fprintf('%-12s %-11s %14s %6s\n', 'case', 'shape', 'curv.ratio', 'peaks');
fprintf('%-12s %-11s %14.2f %6d\n', 'smooth', S1.shape, S1.diag.ratio, S1.npeaks);
fprintf('%-12s %-11s %14.2f %6d\n', 'multimodal', S2.shape, S2.diag.ratio, S2.npeaks);
fprintf('%-12s %-11s %14.1f %6d\n', 'ridge', S3.shape, S3.diag.ratio, S3.npeaks);
fprintf('%-12s %-11s %14s %6s\n', 'degenerate', S4.shape, 'Inf (rank<d)', 1);

%% ---- 6. 出图 ----
% 图1 光滑单峰
h = sim_fig();
contour(S1.A, S1.B, S1.D, [0 1 2 4 6 10 20 40 60 100 160 250], 'b-', 'linewidth', 1);
hold on;
k95 = 2.448;
if exist('chi2inv', 'file') == 2
  try, k95 = sqrt(chi2inv(0.95, 2)); catch, end
end
[Ve,Ee] = eig(S1.cov);
thv = linspace(0, 2*pi, 100);
ell = S1.th_hat*ones(1,100) + k95*Ve*diag(sqrt(diag(Ee)))*[cos(thv); sin(thv)];
plot(ell(1,:), ell(2,:), 'r--', 'linewidth', 1.2);
plot(beta_ols(1), beta_ols(2), 'rp', 'markersize', 12, 'linewidth', 1.5);
grid on;  xlabel('intercept a');  ylabel('slope b');
title('1. Smooth Unimodal Likelihood (Fisher 1922)');
legend('Delta contour', '95% region', 'MLE', 'location', 'northwest');
sim_export(h, figpath('lik_1_smooth'));

% 图2 多峰
h = sim_fig();
contour(S2.A, S2.B, S2.D, [0 1 2 4 8 16 32 64 128 256 512 1024], 'b-', 'linewidth', 1);
hold on;
plot([-2 2], [2 -2], 'rp', 'markersize', 12, 'linewidth', 1.5);
grid on;  xlabel('mu_1');  ylabel('mu_2');
title('2. Bimodal: Label Switching (Richardson & Green 1997)');
legend('Delta contour', 'two equal global modes', 'location', 'northwest');
sim_export(h, figpath('lik_2_multimodal'));

% 图3 山脊
h = sim_fig();
contour(S3.A, S3.B, S3.D, [0 0.5 1 2 4 8 16 32 64 128 256 512 1024], 'b-', 'linewidth', 1);
hold on;
xx = linspace(-1.2, 1.6, 200);
plot(xx, xx.^2, 'r-', 'linewidth', 1.5);
plot(1, 1, 'rp', 'markersize', 12, 'linewidth', 1.5);
grid on;  xlabel('theta_1');  ylabel('theta_2');
title('3. Curved Flat Ridge (Rosenbrock 1960)');
legend('Delta contour', 'valley floor', 'optimum', 'location', 'northwest');
sim_export(h, figpath('lik_3_ridge'));

% 图4 形态判据阈值线
h = sim_fig();
semilogy([1 2], [S1.diag.ratio, S3.diag.ratio], 'bo-', 'linewidth', 1.5, 'markersize', 8);
hold on;
semilogy([0.5 2.5], [30 30], 'r--', 'linewidth', 1.2);
semilogy([0.5 2.5], [1e4 1e4], 'k-.', 'linewidth', 1.2);
grid on;
set(gca, 'xtick', [1 2], 'xticklabel', {'smooth', 'ridge'});
xlim([0.5 2.5]);
ylabel('curvature ratio (max/min eigenvalue)');
title('4. Shape Diagnosis Thresholds (30 / 1e4)');
legend('observed', 'ridge onset (30)', 'degenerate onset (1e4)', 'location', 'northwest');
sim_export(h, figpath('lik_4_thresholds'));

fprintf('\n结论: 同一套极大似然流程, 换个数据就换一种面形;\n');
fprintf('      先诊断形态(曲率比 + 可分辨峰数), 再决定报告方式与求解器。\n');
