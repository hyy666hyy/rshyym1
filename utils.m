function m = lambertian_order(phi_half)
% 计算朗伯发射阶数
% phi_half: 半功率角 (弧度)
m = -log(2) / log(cos(phi_half));
end

function g = optical_concentrator_gain(psi, FOV, n_r)
% 计算光集中器增益
% psi: 入射角 (弧度)
% FOV: 接收视场角 (弧度)
% n_r: 折射系数

if psi <= FOV
    g = (n_r^2) / (sin(FOV)^2);
else
    g = 0;
end
end

function check_led_constraints(W, params)
% 检查LED约束是否满足
w0 = W{1};
total_signal = w0;
for k = 2:length(W)
    total_signal = total_signal + W{k};
end

for i = 1:params.Nt
    current_value = total_signal(i);
    max_value = min(params.I_DC - params.I_L, params.I_U - params.I_DC);
    
    if current_value > max_value
        warning('LED约束不满足: LED %d, 当前值: %.4f, 最大允许值: %.4f', ...
            i, current_value, max_value);
    end
    
    if current_value < 0
        warning('LED非负约束不满足: LED %d, 当前值: %.4f', i, current_value);
    end
end
end