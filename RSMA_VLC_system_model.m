function [sum_rate, user_rates] = RSMA_VLC_system_model(H, W, C, params)
% 计算RSMA-VLC系统的速率
% 输入: H - 信道矩阵 (Nr x Nt x M)
%       W - 预编码矩阵 {w0, w1, ..., wNr}
%       C - 公共速率分配 (Nr x 1)
%       params - 参数
% 输出: sum_rate - 系统和速率
%       user_rates - 各用户速率

Nr = params.Nr;
M = size(H, 3);

% 提取预编码向量
w0 = W{1};  % 公共流预编码
w_private = W(2:end);  % 私有流预编码

% 初始化速率存储
R_kp = zeros(Nr, M);  % 私有速率
R_kc = zeros(Nr, M);  % 公共速率

for m = 1:M
    H_m = H(:, :, m);
    
    for k = 1:Nr
        hk = H_m(k, :)';
        
        % 计算私有速率SINR (公式3.9)
        signal_power = abs(hk' * w_private{k})^2;
        
        % 干扰功率 (其他用户的私有流)
        interference_power = 0;
        for j = 1:Nr
            if j ~= k
                interference_power = interference_power + abs(hk' * w_private{j})^2;
            end
        end
        
        % 噪声功率（方差）：与主程序/优化器统一口径
        if isfield(params, 'noise_power') && ~isempty(params.noise_power)
            noise_power = params.noise_power;
        else
            noise_power = 1;
        end
        
        SINR_kp = signal_power / (interference_power + noise_power);
        R_kp(k, m) = params.B * log2(1 + SINR_kp);
        
        % 计算公共速率SINR (公式3.5)
        signal_power_c = abs(hk' * w0)^2;
        interference_power_c = interference_power + signal_power;  % 所有私有流都是干扰
        
        SINR_kc = signal_power_c / (interference_power_c + noise_power);
        R_kc(k, m) = params.B * log2(1 + SINR_kc);
    end
end

% 取平均速率
R_kp_avg = mean(R_kp, 2);
R_kc_avg = mean(R_kc, 2);

% 公共速率约束 (公式3.7)
R_c_common = min(R_kc_avg);  % 所有用户都能解码的公共速率

% 检查公共速率分配是否合理
if sum(C) > R_c_common
    warning('公共速率分配超过可解码限制，进行缩放');
    scale_factor = R_c_common / sum(C);
    C = C * scale_factor;
end

% 计算总速率 (公式3.11)
user_rates = C + R_kp_avg;

% 系统和速率 (公式3.12)
sum_rate = sum(user_rates);

end