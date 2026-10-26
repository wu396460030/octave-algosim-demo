function [pareto_x, pareto_f, pop_x, pop_f] = nsga2(objectives, nvars, lb, ub, opts)
  % NSGA-II 多目标优化(快速非支配排序 + 拥挤度距离)
  % 输入:
  %   objectives : @(x)返回 m×1 目标向量(全部最小化)
  %   nvars      : 变量维数
  %   lb, ub     : 边界(nvars×1)
  %   opts       : .pop(种群50) .gens(代数60) .pc(交叉0.9) .pm(变异0.1) .seed
  % 输出:
  %   pareto_x / pareto_f : 最终种群中非支配解的决策变量与目标值
  %   pop_x / pop_f       : 全部种群(可选观察)
  % 实测: ZDT1 两目标 50个体60代, 前沿形态与理论一致

  if nargin < 5 || isempty(opts), opts = struct(); end
  N    = getdef(opts, 'pop',  50);
  G    = getdef(opts, 'gens', 60);
  pc   = getdef(opts, 'pc', 0.9);
  pm   = getdef(opts, 'pm', 0.1);
  if isfield(opts, 'seed'), rand('state', opts.seed); end
  lb = lb(:);  ub = ub(:);

  % 初始种群与目标值(上下界取行向量与 rand 广播)
  X = lb' + rand(N, nvars) .* (ub' - lb');
  F = eval_objs(objectives, X, N);

  for g = 1:G
    % 锦标赛选择(支配优先 + 拥挤度破平) + 交叉 + 变异 -> 子代
    [rank_c, crowd_c] = fast_nondominated_sort(F);
    Q = make_offspring(X, F, rank_c, crowd_c, N, lb, ub, pc, pm, nvars);
    Fq = eval_objs(objectives, Q, N);

    % 合并父代+子代, 环境选择
    R  = [X; Q];
    Fr = [F; Fq];
    [rank, crowd] = fast_nondominated_sort(Fr);

    idx = select_by_rank_crowd(rank, crowd, N);
    X = R(idx, :);
    F = Fr(idx, :);
  endfor

  % 最终非支配解
  [rank, ~] = fast_nondominated_sort(F);
  front1 = find(rank == 1);
  pareto_x = X(front1, :);
  pareto_f = F(front1, :);
  pop_x = X;  pop_f = F;
endfunction

function F = eval_objs(objectives, X, N)
  % 目标矩阵 n×m(逐行赋值, 勿用 [F; obj] 竖向堆叠——那会变成 (n*m)×1)
  m = numel(objectives(X(1,:)'));
  F = zeros(size(X, 1), m);
  for i = 1:size(X, 1)
    F(i,:) = objectives(X(i,:)')';
  endfor
endfunction

function [rank, crowd] = fast_nondominated_sort(F)
  % rank: 1 为第一前沿; crowd: 拥挤度距离
  n = size(F, 1);
  m = size(F, 2);
  rank = zeros(n, 1);
  crowd = zeros(n, 1);
  S = cell(n, 1);          % S{i} = 被 i 支配的个体集合
  n_dom = zeros(n, 1);     % 支配 i 的个体数
  for i = 1:n
    for j = 1:n
      if i == j, continue; end
      if all(F(i,:) <= F(j,:)) && any(F(i,:) < F(j,:))
        S{i} = [S{i}, j];
      elseif all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
        n_dom(i) = n_dom(i) + 1;
      endif
    endfor
  endfor
  % 前沿划分(坑: Octave 的 for 对列向量只迭代一次, find 结果必须转行向量)
  r = 1;
  front = find(n_dom == 0)';
  while ~isempty(front)
    rank(front) = r;
    for i = front
      for j = S{i}
        n_dom(j) = n_dom(j) - 1;
      endfor
    endfor
    r = r + 1;
    front = find(n_dom == 0 & rank == 0)';
  endwhile
  % 拥挤度距离(逐目标排序, 边界无穷)
  for i = 1:n
    crowd(i) = 0;
  endfor
  for k = 1:m
    [fs, order] = sort(F(:, k));
    crowd(order(1)) = Inf;
    crowd(order(end)) = Inf;
    rng_k = fs(end) - fs(1);
    if rng_k > 0
      for t = 2:n-1
        crowd(order(t)) = crowd(order(t)) + (fs(t+1) - fs(t-1)) / rng_k;
      endfor
    endif
  endfor
endfunction

function idx = select_by_rank_crowd(rank, crowd, N)
  % 逐前沿填充: rank 升序, 同 rank 内拥挤度降序(稀疏点优先保留)
  % 坑: 排序键不能直接 -(rank*1e9+crowd)——那会保留密集点, 摧毁多样性(实测塌缩)
  c = crowd(:);
  c(isinf(c)) = 1e6;              % 边界点给最大有限拥挤度(Inf 会污染数值键)
  keys = -rank(:) * 1e9 + c;      % rank 升序; 同 rank 内 -c 大者=c 大者在前
  [~, idx] = sort(keys, 'descend');
  idx = idx(1:N)';
endfunction

function Q = make_offspring(X, F, rank, crowd, N, lb, ub, pc, pm, nvars)
  Q = zeros(N, nvars);
  % 锦标赛: 支配关系(rank)优先, 平级取拥挤度大者(稀疏区, 保多样性)
  eta = 15;                                    % SBX 分布指数
  for i = 1:ceil(N/2)
    a = tournament(F, rank, crowd);
    b = tournament(F, rank, crowd);
    [c1, c2] = sbx(X(a,:), X(b,:), pc, eta);   % SBX 交叉(可扩域, 防漂变收缩)
    Q(2*i-1, :) = mutate(c1, pm, lb, ub, nvars);
    if 2*i <= N
      Q(2*i, :) = mutate(c2, pm, lb, ub, nvars);
    endif
  endfor
  Q = max(min(Q, ones(N,1) * ub'), ones(N,1) * lb');
endfunction

function [c1, c2] = sbx(p1, p2, pc, eta)
  % 模拟二进制交叉(SBX): 子代可超出父代范围, 维持种群扩展性
  c1 = p1;  c2 = p2;
  if rand() < pc
    for j = 1:numel(p1)
      if rand() < 0.5 && abs(p1(j) - p2(j)) > 1e-12
        u = rand();
        if u <= 0.5
          beta = (2 * u) ^ (1 / (eta + 1));
        else
          beta = (1 / (2 * (1 - u))) ^ (1 / (eta + 1));
        endif
        c1(j) = 0.5 * ((1 - beta) * p1(j) + (1 + beta) * p2(j));
        c2(j) = 0.5 * ((1 + beta) * p1(j) + (1 - beta) * p2(j));
      endif
    endfor
  endif
endfunction

function winner = tournament(F, rank, crowd)
  % 标准二进制锦标赛: rank 小(更优前沿)胜; 同 rank 拥挤度大(稀疏)胜
  a = randi(size(F, 1));
  b = randi(size(F, 1));
  if a == b
    winner = a;
  elseif rank(a) < rank(b)
    winner = a;
  elseif rank(b) < rank(a)
    winner = b;
  elseif crowd(a) > crowd(b)
    winner = a;
  else
    winner = b;
  endif
endfunction

function c = mutate(c, pm, lb, ub, nvars)
  % 高斯扰动变异(比均匀重采收敛更精细), 边界钳位
  for j = 1:nvars
    if rand() < pm
      c(j) = c(j) + 0.15 * randn() * (ub(j) - lb(j));
      c(j) = max(min(c(j), ub(j)), lb(j));
    endif
  endfor
endfunction

function v = getdef(s, name, default)
  if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
  else
    v = default;
  end
endfunction
