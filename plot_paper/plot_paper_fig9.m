% Plot Y direction pixel error distribution curve, single plot
% required parameters: n_ranges, ranges
dist = pixErrDistY;

figure
x_ = -40:0.1:60;
legend_val = {'range 1: top left';
              'range 2: top center';
              'range 3: top right';
              'range 4: bottom left';
              'range 5: bottom center';
              'range 6: bottom right';};
%color = jet(7);
%color = [color(2:4,:); color(end-3:end,:)];
for i = 1:n_ranges
    % plot PDF
    y_ = pdf(dist{i}, x_);
    plot(x_, y_, 'LineWidth', 2);
    hold on
    % plot mu
    ver_line = dist{i}.mu;
    % xline(ver_line, '--k');   % only in MATLAB 2018b or higher
end
title('Fitted PDF for horizontal ranges');
ylabel('probability density');
legend(legend_val,'Location', 'NorthEast');
xlabel('\gamma pixel error');
grid on
