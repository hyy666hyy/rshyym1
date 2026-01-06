function [sum_rate, W_opt, C_opt] = RSMA_VLC_optimization(H_est, H_error, params, P_t, noise_power)
% RSMA-VLC系统和速率预编码优化算法 (算法3.1)
% 基于松弛变量和连续凸近似(SCA)
%
% 口径统一说明（避免混用“噪声=1”与“SNR(dB)”）：
% - 本函数不再直接使用/传入 SNR(dB)。主程序在每个 SNR(dB) 点上先计算:
%   * 噪声方差 noise_power（由系统参数得到）
%   * 对应的总发射功率约束 P_t
% - 预编码向量 W 的总功率约束: ||w0||^2 + sum_k ||wk||^2 <= P_t
% - 所有 SINR 的噪声项统一使用 noise_power（而不是写死 1）
% - 不完美 CSIT 的误差强度用“等效线性 SNR”驱动，并与 calculate_channel.m 的 ref_SNR 对齐

Nr = params.Nr;
Nt = params.Nt;
M = params.M;
max_iter = params.max_iter;
tol = params.tol;

% 噪声功率（方差）
if nargin < 5 || isempty(noise_power)
    noise_power = 1;
end

% 统一口径：用等效线性SNR = (P_t * E[|h|^2]) / noise_power
% 用于CSIT误差随SNR变化的缩放（与 calculate_channel.m 的 ref_SNR 对齐）
avg_h2 = mean(abs(H_est(:)).^2);
effective_snr_linear = (P_t * avg_h2) / max(noise_power, eps);
P_e_current = effective_snr_linear^(-0.6);
P_e_ref = params.ref_SNR^(-0.6);
error_scale = sqrt(P_e_current / max(P_e_ref, eps));

% 生成真实信道样本
H_real = zeros(Nr, Nt, M);
for m = 1:M
    H_real(:, :, m) = H_est + error_scale * H_error(:, :, m);
end

% 初始化预编码矩阵
% 使用随机初始化
W_current = cell(1, Nr+1);  % w0, w1, ..., wNr
W_current{1} = randn(Nt, 1);  % 公共流 (VLC: 实值)
for k = 1:Nr
    W_current{k+1} = randn(Nt, 1);  % 私有流 (VLC: 实值)
end

% 归一化功率
total_power = norm(W_current{1})^2 + sum(arrayfun(@(k) norm(W_current{k+1})^2, 1:Nr));
scale_factor = sqrt(P_t / max(total_power, eps));
W_current{1} = W_current{1} * scale_factor;
for k = 1:Nr
    W_current{k+1} = W_current{k+1} * scale_factor;
end

% 初始化公共速率分配
C_current = ones(Nr, 1) * 0.01;

% 初始化松弛变量（关键：g_p/g_c 必须与当前 W、噪声同量纲，否则 SCA 下界会极度保守导致速率≈0）
[a_p, f_p, g_p, a_c, f_c, g_c, C_current] = initialize_slack_variables(H_est, W_current, C_current, params, noise_power);

% 迭代优化
sum_rate_prev = 0;
iteration = 0;
converged = false;

fprintf('  开始SCA迭代优化...\n');
while ~converged && iteration < max_iter
    iteration = iteration + 1;

    % 解决凸优化子问题 (问题3.20)
    [W_new, C_new, a_p_new, f_p_new, g_p_new, a_c_new, f_c_new, g_c_new, ~] = ...
        solve_sca_subproblem(H_est, W_current, C_current, a_p, f_p, g_p, a_c, f_c, g_c, ...
                           params, P_t, noise_power);
    
    % 更新变量
    W_current = W_new;  % 预编码矩阵
    C_current = C_new;  % 公共速率分配
    a_p = a_p_new; f_p = f_p_new; g_p = g_p_new;
    a_c = a_c_new; f_c = f_c_new; g_c = g_c_new;
    
    % 计算实际系统和速率
    current_sum_rate = calculate_sum_rate(H_real, W_current, C_current, params, noise_power);
    
    % 检查收敛
    if abs(current_sum_rate - sum_rate_prev) < tol
        converged = true;
        fprintf('    迭代 %d: 收敛，系统和速率 = %.4f\n', iteration, current_sum_rate);
    else
        sum_rate_prev = current_sum_rate;
        if mod(iteration, 10) == 0
            fprintf('    迭代 %d: 系统和速率 = %.4f\n', iteration, current_sum_rate);
        end
    end
end

% 最终结果
sum_rate = current_sum_rate;
W_opt = W_current;
C_opt = C_current;

fprintf('  优化完成，迭代次数: %d，系统和速率: %.4f\n', iteration, sum_rate);

end

%% 辅助函数
function [a_p, f_p, g_p, a_c, f_c, g_c, C_out] = initialize_slack_variables(H_est, W, C_in, params, noise_power)
% 初始化松弛变量
% 目的：给 SCA 的线性化点提供“合理尺度”的 (w_old, g_old)。
% 若 g_old 初值取 1 而实际噪声/干扰是 1e-12 量级，会让 |h^T w|^2 / g 的下界≈0，导致 a≈0、C≈0，算法很快“收敛到0”。

Nr = params.Nr;
B = params.B;

w0 = W{1};
w_private = W(2:end);

g_p = zeros(Nr, 1);
g_c = zeros(Nr, 1);
a_p = zeros(Nr, 1);
a_c = zeros(Nr, 1);
f_p = ones(Nr, 1);
f_c = ones(Nr, 1);

R_c_all = zeros(Nr, 1);

for k = 1:Nr
    hk = H_est(k, :)';
    % 干扰（私有流）
    interf_p = 0;
    total_priv = 0;
    for j = 1:Nr
        pj = abs(hk' * w_private{j})^2;
        total_priv = total_priv + pj;
        if j ~= k
            interf_p = interf_p + pj;
        end
    end

    % g 的合理初值：干扰+噪声
    g_p(k) = interf_p + noise_power;
    g_c(k) = total_priv + noise_power;

    % a 的合理初值：当前点的速率（作为下界变量的起点）
    sig_p = abs(hk' * w_private{k})^2;
    sinr_p = sig_p / max(g_p(k), eps);
    a_p(k) = B * log2(1 + sinr_p);

    sig_c = abs(hk' * w0)^2;
    sinr_c = sig_c / max(g_c(k), eps);
    R_c_all(k) = B * log2(1 + sinr_c);
    a_c(k) = R_c_all(k);

    % f 的一致初值：f = 2^(a/B)
    f_p(k) = exp((log(2)/B) * a_p(k));
    f_c(k) = exp((log(2)/B) * a_c(k));
end

% 让初始 C 可行：sum(C) <= min_k a_c(k)
R_c_common = min(R_c_all);
C_total = sum(C_in);
if C_total > R_c_common && C_total > 0
    C_out = C_in * (R_c_common / C_total);
else
    C_out = C_in;
end

end

function [tau_p, tau_c] = compute_taylor_expansion(H_real, W, g_p_old, g_c_old, params)
% 计算一阶泰勒展开项 (公式3.15和3.18)

Nr = params.Nr;
% Nt = params.Nt;          % 维度信息（此函数内部未直接使用，仅用于对照论文/调试）
% M  = size(H_real, 3);    % 样本数（此函数内部使用 mean(H_real,3) 已隐式用到）

% 预编码向量：W = {w0, w1, ..., wNr}
% w0: 公共流预编码；wk: 用户k私有流预编码
w0 = W{1};
w_private = W(2:end);

% 初始化泰勒展开项
tau_p = zeros(Nr, 1);
tau_c = zeros(Nr, 1);

for k = 1:Nr
    % 计算平均信道
    H_avg = mean(H_real, 3);
    hk = H_avg(k, :)';
    
    % 计算私有流的泰勒展开 (公式3.15)
    numerator = 2 * (w_private{k}' * (hk * hk') * w_private{k});
    denom = g_p_old(k);
    tau_p(k) = numerator / denom - (abs(hk' * w_private{k}) / denom)^2;
    
    % 计算公共流的泰勒展开 (公式3.18)
    numerator_c = 2 * (w0' * (hk * hk') * w0);
    denom_c = g_c_old(k);
    tau_c(k) = numerator_c / denom_c - (abs(hk' * w0) / denom_c)^2;
end

end

function sum_rate = calculate_sum_rate(H_real, W, C, params, noise_power)
% 计算RSMA-VLC系统的实际和速率
% 输入:
%   H_real - 真实信道矩阵 (Nr x Nt x M)
%   W - 预编码矩阵 cell {w0, w1, ..., wNr}
%   C - 公共速率分配向量 (Nr x 1)
%   params - 参数结构体
%   noise_power - 噪声功率
% 输出:
%   sum_rate - 系统和速率 (bits/s/Hz)

Nr = params.Nr;          % 用户数
M = size(H_real, 3);     % 蒙特卡洛样本数
B = params.B;            % 系统带宽

% 提取预编码向量
w0 = W{1};               % 公共流预编码
w_private = W(2:end);    % 私有流预编码 cell数组

% 初始化总速率
total_sum_rate = 0;

% ===================== 核心计算循环 =====================
for m = 1:M  % 遍历所有信道样本
    H_m = H_real(:, :, m);  % 第m个信道样本
    
    % 计算所有用户的公共速率（解码公共流s0）
    R_c_all = zeros(Nr, 1);  % 各用户能解码的公共速率
    
    for k = 1:Nr  % 遍历所有用户
        hk = H_m(k, :)';  % 用户k的信道向量
        
        % ========== 步骤1: 计算公共流SINR ==========
        % 公式(3.5): γ_{k,c} = |h_k^T w_0|² / (Σ_{i=1}^{Nr} |h_k^T w_i|² + σ²)
        
        % 公共信号功率
        signal_c = abs(hk' * w0)^2;  % |h_k^T w_0|²
        
        % 干扰功率（所有私有流都是干扰）
        interference_c = 0;
        for j = 1:Nr
            interference_c = interference_c + abs(hk' * w_private{j})^2;
        end
        
        % 公共流SINR
        SINR_c = signal_c / (interference_c + noise_power);
        
        % 公共速率 (公式3.6)
        R_c_all(k) = B * log2(1 + SINR_c);
    end
    
    % ========== 步骤2: 确定可实现的公共速率 ==========
    % 公式(3.7): R_c = min{R_{1,c}, R_{2,c}, ..., R_{N_r,c}}
    R_c_common = min(R_c_all);  % 所有用户都能解码的公共速率
    
    % ========== 步骤3: 检查公共速率分配 ==========
    % 总分配的公共速率不能超过可实现的公共速率
    C_total = sum(C);
    if C_total > R_c_common
        % 按比例缩放
        scale_factor = R_c_common / C_total;
        C_scaled = C * scale_factor;
    else
        C_scaled = C;
    end
    
    % ========== 步骤4: 计算各用户的私有速率 ==========
    user_rates = zeros(Nr, 1);  % 各用户总速率
    
    for k = 1:Nr
        hk = H_m(k, :)';
        
        % 公式(3.9): γ_{k,p} = |h_k^T w_k|² / (Σ_{j≠k} |h_k^T w_j|² + σ²)
        
        % 私有信号功率
        signal_p = abs(hk' * w_private{k})^2;  % |h_k^T w_k|²
        
        % 干扰功率（其他用户的私有流）
        interference_p = 0;
        for j = 1:Nr
            if j ~= k
                interference_p = interference_p + abs(hk' * w_private{j})^2;
            end
        end
        
        % 私有流SINR
        SINR_p = signal_p / (interference_p + noise_power);
        
        % 私有速率 (公式3.10)
        R_p = B * log2(1 + SINR_p);
        
        % 用户总速率 (公式3.11)
        user_rates(k) = C_scaled(k) + R_p;
    end
    
    % ========== 步骤5: 累加该样本的系统和速率 ==========
    total_sum_rate = total_sum_rate + sum(user_rates);
    
end

% ========== 步骤6: 计算平均系统和速率 ==========
sum_rate = total_sum_rate / M;

end

function [W_opt, C_opt, a_p_opt, f_p_opt, g_p_opt, a_c_opt, f_c_opt, g_c_opt, obj_value] = ...
    solve_sca_subproblem(H_est, W_current, C_current, a_p_old, f_p_old, g_p_old, a_c_old, f_c_old, g_c_old, ...
                       params, P_t, noise_power)
% 解决SCA凸优化子问题 (对应问题3.20)

Nr = params.Nr;
Nt = params.Nt;
B = params.B;
max_power = P_t;

% 上一次迭代点（用于一阶近似/切线下界）
w0_old = W_current{1};
w_private_old = zeros(Nt, Nr);
for k = 1:Nr
    w_private_old(:, k) = W_current{k+1};
end

cvx_begin quiet
    % 变量声明
    variable w0(Nt)
    variable w_private(Nt, Nr)
    variable C(Nr) nonnegative
    variable a_p(Nr) nonnegative
    variable f_p(Nr) nonnegative
    variable g_p(Nr) nonnegative
    variable a_c(Nr) nonnegative
    variable f_c(Nr) nonnegative
    variable g_c(Nr) nonnegative
    
    % 目标函数：最大化系统和速率 (公式3.20a)
    maximize( sum(C) + sum(a_p) )
    
    subject to
        % 功率约束
        % CVX 不允许对凸表达式再做 .^2（会触发 DCP 规则错误）。
        % 这里等价写成“所有元素平方和”形式：||w0||_2^2 + ||w_private||_F^2 <= P_t
        sum_square(w0) + sum_square(w_private(:)) <= max_power;

        % 约束(3.13d): LED线性工作区域
        % 论文形式：\sum_{j=0}^{Nr} |w_j^T e_i| <= min{I_DC-I_L, I_U-I_DC}
        % 这里 e_i 选择第 i 个LED，对应 |w_j(i)|。
        for i = 1:Nt
            sum(abs([w0(i); w_private(i, :)'])) <= min(params.I_DC - params.I_L, params.I_U - params.I_DC);
        end
        
        % 私有速率约束 (公式3.14b, 3.14d, 3.16)
        for k = 1:Nr
            % 约束(3.14b): f_p >= 2^(a_p/B)
            % CVX 中避免使用 2^x（会走 power/exp 内部路径且更容易触发 DCP/依赖问题）
            % 等价写法：2^(a_p/B) = exp( log(2) * a_p / B )
            exp((log(2)/B) * a_p(k)) <= f_p(k);
            
            % 约束(3.14d): g_p >= 干扰+噪声
            hk = H_est(k, :)';
            A = hk * hk';
            interference_power = 0;
            for j = 1:Nr
                if j ~= k
                    interference_power = interference_power + quad_form(w_private(:, j), A);
                end
            end
            g_p(k) >= interference_power + noise_power;

            % 约束(3.16): 一阶下界近似后的 SINR 约束
            % 目标：f_p - 1 <= |h_k^T w_k|^2 / g_{k,p}
            % 右侧为凸函数的“下界切线”（仿射），在(w_old,g_old)处：
            %   |h^T w|^2/g >= 2*w_old^T A w / g_old - (w_old^T A w_old)/g_old^2 * g
            g_old = max(g_p_old(k), 1e-9);
            num_old = w_private_old(:, k)' * A * w_private_old(:, k);
            (2 * (w_private_old(:, k)' * A * w_private(:, k))) / g_old - (num_old / (g_old^2)) * g_p(k) >= f_p(k) - 1;
        end
        
        % 公共速率约束 (公式3.17a, 3.17b, 3.17c, 3.17e, 3.19)
        for k = 1:Nr
            % 约束(3.17a): 公共速率分配约束
            sum(C) <= a_c(k);
            
            % 约束(3.17b): 速率下界
            % 这里简化处理，实际应该是a_c <= R_k,c
            a_c(k) >= 0;
            
            % 约束(3.17c): f_c >= 2^(a_c/B)
            exp((log(2)/B) * a_c(k)) <= f_c(k);
            
            % 约束(3.17e): g_c >= 总干扰+噪声
            hk = H_est(k, :)';
            A = hk * hk';
            total_interference = 0;
            for j = 1:Nr
                total_interference = total_interference + quad_form(w_private(:, j), A);
            end
            g_c(k) >= total_interference + noise_power;

            % 约束(3.19): 一阶下界近似后的公共 SINR 约束
            % 目标：f_c - 1 <= |h_k^T w_0|^2 / g_{k,c}
            g_old = max(g_c_old(k), 1e-9);
            num_old = w0_old' * A * w0_old;
            (2 * (w0_old' * A * w0)) / g_old - (num_old / (g_old^2)) * g_c(k) >= f_c(k) - 1;
        end
        
cvx_end

% 打包输出
W_opt = cell(1, Nr+1);
W_opt{1} = w0;
for k = 1:Nr
    W_opt{k+1} = w_private(:, k);
end

C_opt = C;
a_p_opt = a_p; f_p_opt = f_p; g_p_opt = g_p;
a_c_opt = a_c; f_c_opt = f_c; g_c_opt = g_c;
obj_value = cvx_optval;

end