function [H, H_est, H_error] = calculate_channel(params)
% VLC信道模型计算
% 输入: params - 参数结构体
% 输出: H - 真实信道矩阵 (Nr x Nt x M)
%       H_est - 估计信道矩阵 (Nr x Nt)
%       H_error - 信道误差矩阵 (Nr x Nt x M)

Nr = params.Nr;
Nt = params.Nt;
M = params.M;

% 初始化信道矩阵
H = zeros(Nr, Nt, M);
H_est = zeros(Nr, Nt);

% 计算朗伯发射阶数
m = -log(2) / log(cos(params.phi_half_rad));

% 计算每个用户和每个LED之间的信道
for k = 1:Nr  % 遍历用户
    user_pos = params.user_positions(k, :);
    
    for i = 1:Nt  % 遍历LED
        led_pos = params.led_positions(i, :);
        
        % 计算距离
        d = norm(led_pos - user_pos);
        
        % 计算发射角 (LED到用户的连线与LED法线的夹角)
        % 假设LED垂直向下照射
        led_normal = [0, 0, -1];  % LED法线方向 (垂直向下)
        vec_user = user_pos - led_pos;  % LED到用户的向量
        
        % 计算夹角
        cos_phi = abs(dot(led_normal, vec_user) / (norm(led_normal) * norm(vec_user)));
        phi = acos(cos_phi);
        
        % 计算接收角 (假设PD垂直向上)
        pd_normal = [0, 0, 1];  % PD法线方向 (垂直向上)
        vec_led = led_pos - user_pos;  % 用户到LED的向量
        cos_psi = abs(dot(pd_normal, vec_led) / (norm(pd_normal) * norm(vec_led)));
        psi = acos(cos_psi);
        
        % 检查是否在视场角内
        if psi <= params.FOV_rad
            % 计算光集中器增益
            g_psi = (params.nr^2) / (sin(params.FOV_rad)^2);
        else
            g_psi = 0;
        end
        
        % 计算信道增益 (公式2.1)
        h = ((m+1)*params.A / (2*pi*d^2)) * ...
            (cos(phi)^m) * params.Tf * params.roe * g_psi * cos(psi);
        
        H_est(k, i) = h;
    end
end

% 生成信道误差 (用于不完美CSIT)
% 假设误差方差与SNR的-0.6次方成正比
P_e = params.ref_SNR^(-0.6);
H_error = sqrt(P_e/2) * (randn(Nr, Nt, M) + 1j*randn(Nr, Nt, M));

% 生成真实信道 (估计信道 + 误差)
for m_idx = 1:M
    H(:, :, m_idx) = H_est + H_error(:, :, m_idx);
end

% 由于VLC信道是实值的，取实部
H_est = real(H_est);
H = real(H);
H_error = real(H_error);

% 打印信道矩阵到终端
fprintf('\n========== 信道矩阵信息 ==========\n');
fprintf('估计信道矩阵 H_est (%d x %d):\n', Nr, Nt);
disp(H_est);

fprintf('真实信道矩阵 H (%d x %d x %d):\n', Nr, Nt, M);
for m_idx = 1:1
    fprintf('  第 %d 个样本:\n', m_idx);
    disp(H(:, :, m_idx));
end

fprintf('信道误差矩阵 H_error (%d x %d x %d):\n', Nr, Nt, M);
for m_idx = 1:1
    fprintf('  第 %d 个样本:\n', m_idx);
    disp(H_error(:, :, m_idx));
end
fprintf('==================================\n\n');

end