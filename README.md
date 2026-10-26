# MATLAB/Octave 仿真与算法助手（开源演示版）
# MATLAB/Octave Simulation & Algorithm Assistant (Open Demo)

> 论文级参数、可直接运行、结果可与理论对照的工程仿真代码集。
> Paper-grade parameters, directly runnable, results you can check against theory.

![NSGA-II Pareto Front](demos/figures/nsga2_front.png)

## 这是什么 / What is this

本仓库是 **「MATLAB/Octave 仿真与算法助手」技能的开源演示版**：从完整版的 24 个方法族、
20 篇论文反向测试（20/20 通过）与 13 项回归测试中，精选 3 个自包含示例供你直接运行体验。

This repo is the **open demo of the "MATLAB/Octave Simulation & Algorithm Assistant" skill**:
3 self-contained examples distilled from the full version (24 method families, 20/20 paper-based
reverse tests, 13/13 regression tests).

- ✅ **Octave 11.x 与 MATLAB R2016b+ 双兼容** / Works on both Octave and MATLAB
- ✅ 每个脚本单文件自包含，复制即跑 / Each script is a single self-contained file
- ✅ 参数带论文/标准出处，结果有理论对照 / Paper-sourced parameters, theory-anchored results

## 快速开始 / Quick Start

```bash
# Octave（开源免费）/ Octave (open source)
octave-cli demos/demo_pid_tuning.m --no-gui

# MATLAB 中直接运行脚本即可 / In MATLAB just run the script
```

## 三个精选示例 / Three Featured Demos

| Demo | 方法族 Family | 亮点 Highlight | 效果图 Figure |
|---|---|---|---|
| `demo_pid_tuning.m` | 控制 PID Control | P/PI/PID 三控制器闭环对比，指标表输出 | ![PID](demos/figures/pid_tuning.png) |
| `demo_mm1_callcenter.m` | 排队论+蒙特卡洛 Queueing + MC | 呼叫中心案例复现论文理论值，**误差 0.8%** | ![M/M/1](demos/figures/mm1_waitdist.png) |
| `demo_nsga2_zdt1.m` + `nsga2.m` | 多目标优化 Multi-objective | NSGA-II 前沿铺满 [0,1]，贴合理论曲线偏差 0.0000 | ![NSGA-II](demos/figures/nsga2_front.png) |

运行方式 / Run: `octave-cli demos/demo_nsga2_zdt1.m --no-gui`（其余同理 / likewise）

## 完整版有什么 / What's in the Full Version

开源演示只是冰山一角。完整技能覆盖 **24 个方法族**，关键词直达正确方法：
The open demos are just the tip. The full skill covers **24 method families**, routed by keyword:

- **智能优化 Intelligent optimization**: GA / PSO / 模拟退火 SA / 蚁群 ACO / NSGA-II / 网格搜索 / 穷举
- **图论与组合 Graph & combinatorial**: 最短路 Dijkstra / 网络流 Max-flow / 二分图匹配 Hungarian / DP / 回溯 / 分支定界
- **仿真 Simulation**: 蒙特卡洛 / 离散事件 / 多智能体 Boids / 代理模型 RBF-Kriging / 数字孪生（卡尔曼同化）
- **学习 Learning**: BP 神经网络 / LMS / RLS / 在线自适应
- **数学规划 Math programming**: 线性 glpk / 整数 / 二次 qp / 多目标
- **拟合与数值 Fitting & numerics**: leasqr / lsqnonlin / 插值样条 / ODE / 图像处理

### 质量证据 / Quality evidence

- 20 篇论文/标准案例反向测试 **20/20 通过**（呼叫中心误差 0.8%、CWRU 轴承 BSF 突出 8.7×、
  LP 影子价格验证精确命中、GA/SA 均达 TSP 穷举最优……）
- 自写算法库回归测试 **13/13 PASS**（每个算法配已知解判定）
- 20/20 paper-based reverse tests passed; 13/13 regression tests on the in-house algorithm library.

## 获取完整版 / Get the Full Version

完整版上架于 WorkBuddy 技能市场（含 24 方法族路由、15 个实测算法、全部参考文档）。
The full skill is published on the WorkBuddy skill marketplace.

> 🔗 **WorkBuddy 客户端 → 技能市场 → 搜索「MATLAB/Octave 仿真与算法助手」**
> （链接上架后更新 / link to be updated after publishing）

## 环境要求 / Requirements

- GNU Octave 11.x（开源）或 MATLAB R2016b+
- 控制类示例需 `control` 包 / 控制示例 need the `control` package: `pkg install -forge control`
- 其余示例仅依赖核心函数 / Other demos use core functions only

## License

MIT © 2026 why（见 [LICENSE](LICENSE)）。完整版技能的版权归作者所有，未经许可不得二次分发。
MIT for this demo repo. The full skill is proprietary; redistribution requires permission.
