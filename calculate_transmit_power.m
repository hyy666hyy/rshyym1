function [P_t, noise_power] = calculate_transmit_power(target_SNR_dB, H_est, params)
% 计算VLC系统发射功率

% 目标SNR（线性）
SNR_linear = 10^(target_SNR_dB/10);

% 噪声功率（根据系统参数计算）
% 对于VLC，噪声包括：
% 1. 散粒噪声 (shot noise)
% 2. 热噪声 (thermal noise)
% 3. 背景光噪声
B = params.B;  % 带宽
q = 1.6e-19;   % 电子电荷
R = params.roe; % 光电响应度
I_bg = 1e-3;   % 背景光电流（示例值）

% 散粒噪声方差
sigma_shot2 = 2 * q * R * I_bg * B;

% 热噪声方差（简化）
sigma_thermal2 = 1e-12;  % 示例值

% 总噪声功率
noise_power = sigma_shot2 + sigma_thermal2;

% 平均信道增益
avg_h = mean(abs(H_est(:)).^2);

% 避免 avg_h≈0 导致 P_t 变为 Inf/NaN（例如用户/LED 超出FOV导致信道接近0）
avg_h_floor = 1e-15;
if ~isfinite(avg_h) || avg_h < avg_h_floor
	warning('平均信道增益 avg_h=%.3e 过小/无效，已夹紧到 %.1e 以避免 P_t 变为 Inf。请检查场景/FOV/位置设置。', avg_h, avg_h_floor);
	avg_h = avg_h_floor;
end

% 所需发射功率
P_t = SNR_linear * noise_power / avg_h;

end