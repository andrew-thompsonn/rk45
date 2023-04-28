clc; clear all; close all;

%%
RE = 6378e3;

data = importdata("Optimum_1_100p000_20p000");
figure(); hold on;

plot(data(:, 2),  data(:, 3),  'r', 'Linewidth', 2);
plot(data(:, 10), data(:, 11), '--', 'Color', [0.4 0.4 0.4], 'Linewidth', 2);

theta = linspace(0, 2*pi);
x = RE*cos(theta);
y = RE*sin(theta);
plot(x, y, 'b', 'LineWidth', 3);

plot(data(end, 10), data(end, 11), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3.5);
plot(0, 0, 'bo', 'MarkerFaceColor', 'b', 'Markersize', 4.5);

legend(["Spacecraft", "Moon", "Earth"]);

xlabel("X [m]");
ylabel("Y [m]");

grid minor;
axis equal;

%% ./exe_three_body 1 100 5

% rk45  - 10,442 steps (156.6 seconds) | TIME_STEP = 10s | TOL = 2E3 | 
% euler - 27,888 steps (71.20 seconds) | TIME_STEP = 10s |

% rk45  - 10,443 steps (154.0 seconds) | TIME_STEP = 5s  | TOL = 2E3 |
% euler - 27,890 steps (215.8 seconds) | TIME_STEP = 5s  |

%% ./exe_three_body 1 100 5 (All files compiled with level 2 gcc optimization)

% rk45  - 10,443 steps (57.42 seconds) | TIME_STEP = 5s  | TOL = 2E3 |
% euler - 27,890 steps (76.29 seconds) | TIME_STEP = 5s  |