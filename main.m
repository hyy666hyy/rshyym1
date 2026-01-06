%% RSMA-VLC系统和速率预编码优化仿真主程序
% 复现论文第三章图3.3和图3.4的RSMA方案
% 作者：王谦，南昌大学

clc; clear; close all;
addpath(genpath(pwd));

%% 仿真参数设置
fprintf('========== RSMA-VLC系统和速率仿真 ==========\n');

% 基本参数
params.Nt = 4;               % LED数量 (发射天线)
params.Nr = 2;               % 用户数量
params.B = 1;                % 系统带宽 (Hz)
params.snr_dB = 0:5:30;      % 信噪比范围 (dB)
params.num_snr = length(params.snr_dB);

% 房间尺寸 (7x7x5 m^3)
params.room_size = [7, 7, 5];  % [长, 宽, 高]

% LED安装高度 (天花板)
params.led_height = 4.5;     % m

% 用户高度 (桌面)
params.user_height = 0.8;    % m (假设桌子高度)

% VLC系统参数 (表3.2)
params.A = 1e-4;             % PD检测面积 (m^2)
params.Tf = 1;               % 光滤波器增益
params.roe = 0.54;           % 光电转换效率常数 (A/W)
params.phi_half = 60;        % LED半功率角 (度)
params.phi_half_rad = deg2rad(params.phi_half);
params.FOV = 60;             % 接收视场角 (度)
params.FOV_rad = deg2rad(params.FOV);
params.nr = 1.5;             % 折射系数

% LED线性工作区参数
params.I_DC = sqrt(6);       % 直流偏置 (A)
params.I_U = 5;              % 线性区域上界 (A)
params.I_L = 1;              % 线性区域下界 (A)

% 算法参数
params.max_iter = 500;       % 最大迭代次数
params.tol = 1e-5;           % 收敛精度
params.M = 100;             % 蒙特卡洛样本数 (用于不完美CSIT)

% 添加参考SNR值用于信道误差生成
% 使用第一个SNR值作为参考
params.ref_SNR = 10^(params.snr_dB(1)/10);  % 参考SNR值


% 定义仿真场景 (参考图3.2)
fprintf('选择仿真场景:\n');
fprintf('1. 欠载情况 - 四个LED服务两个用户\n');
fprintf('2. 过载情况 - 两个LED服务三个用户\n');
scene_choice = input('请输入场景选择 (1或2): ');

if scene_choice == 1
    params.Nt = 4;
    params.Nr = 2;
    fprintf('选择场景1: 四个LED服务两个用户\n');
    
    % 定义LED和用户位置 (示例位置)
    params.led_positions = [3, 5.0; 5.0, 5.0; 3, 3; 5.0, 3.0]; % (x,y)
    params.user_positions = [3.5, 4; 4.5, 4.0]; % 用户1和用户2
    
elseif scene_choice == 2
    params.Nt = 2;
    params.Nr = 3;
    fprintf('选择场景2: 两个LED服务三个用户\n');
    
    % 定义LED和用户位置
    params.led_positions = [3.0, 4.0; 5.0, 4.0];
    params.user_positions = [3.5, 4.0; 4.0, 4.0; 4.5, 4.0]; 
else
    error('无效的场景选择');
end

% 添加高度信息
params.led_positions = [params.led_positions, params.led_height * ones(size(params.led_positions,1),1)]; 
params.user_positions = [params.user_positions, params.user_height * ones(size(params.user_positions,1),1)];

%% 计算信道矩阵
fprintf('计算VLC信道矩阵...\n');
[H, H_est, H_error] = calculate_channel(params);

%% 初始化结果存储
sum_rates_rsma = zeros(1, params.num_snr);

%% 主仿真循环
fprintf('\n开始RSMA-VLC系统和速率仿真...\n');
for snr_idx = 1:params.num_snr
    fprintf('处理SNR = %d dB (%d/%d)...\n', params.snr_dB(snr_idx), snr_idx, params.num_snr);
    
    % 当前目标SNR
    target_SNR_dB = params.snr_dB(snr_idx);

    % 计算所需发射功率
    % 使用VLC噪声模型 + 平均信道增益进行估算
    [P_t, noise_power] = calculate_transmit_power(target_SNR_dB, H_est, params);
    params.noise_power = noise_power;

    % 当前SNR对应的功率约束
    params.Pt = P_t;  % 总发射功率
    
    % 运行RSMA-VLC优化算法
    [sum_rate, W_opt, C_opt] = RSMA_VLC_optimization(H_est, H_error, params, P_t, noise_power);
    
    % 存储结果
    sum_rates_rsma(snr_idx) = sum_rate;
    
    fprintf('  系统和速率 = %.4f bits/s/Hz\n', sum_rate);
end

%% 保存结果
save('rsma_vlc_results.mat', 'params', 'sum_rates_rsma');

%% 绘制结果
plot_results(params, sum_rates_rsma);

fprintf('\n仿真完成！\n');