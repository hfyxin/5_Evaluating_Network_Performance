% Plot 1 histogram and 1 PDF on the same graph

% First, need the distribution array and PDF already in workspace, e.g.
data = pixErrX{1};      % vector
dist = pixErrDistX{1};  % Normal Distribution

% Plot histogram
figure
yyaxis left
histogram(data,'EdgeColor','k'); % 'FaceColor','b'
xlabel('\gamma pixel error (px)');
ylabel('# of detections');
grid on
hold on
ax = gca;
ax.YColor = 'black';   % y axis color

% PDF function
x_ = -40:0.1:60;
y_ = pdf(dist, x_);
yyaxis right
plot(x_, y_, 'LineWidth', 2, 'Color', 'b');
ylabel('probability density');
ax.YColor = 'black';    % y axis color

xlim([-20 20]);     % x axis limit
legend({'Data points';'Fitted PDF'},'Location', 'NorthEast');