function plot_results(params, sum_rates_rsma)
% 绘制RSMA-VLC系统和速率仿真结果

figure('Position', [100, 100, 800, 600]);

% 绘制系统和速率 vs SNR
plot(params.snr_dB, sum_rates_rsma, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b');
hold on;
grid on;

% 添加参考线 (理论极限)
% 这里可以添加其他方案的曲线进行比较

xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('系统和速率 (bits/s/Hz)', 'FontSize', 12, 'FontWeight', 'bold');

if params.Nt == 4 && params.Nr == 2
    title('RSMA-VLC系统和速率 (欠载情况: 4LEDs × 2用户)', 'FontSize', 14, 'FontWeight', 'bold');
elseif params.Nt == 2 && params.Nr == 3
    title('RSMA-VLC系统和速率 (过载情况: 2LEDs × 3用户)', 'FontSize', 14, 'FontWeight', 'bold');
else
    title('RSMA-VLC系统和速率', 'FontSize', 14, 'FontWeight', 'bold');
end

legend('RSMA-VLC', 'Location', 'northwest', 'FontSize', 10);

% 设置坐标轴
xlim([min(params.snr_dB), max(params.snr_dB)]);
set(gca, 'FontSize', 11, 'GridLineStyle', '--', 'GridAlpha', 0.3);

% 保存图像
saveas(gcf, 'rsma_vlc_sum_rate.png');
fprintf('结果图像已保存为 rsma_vlc_sum_rate.png\n');

end