y0 = [2500; 0; 0; 0];
t_range = [0 50];

[t, y] = ode45(@case1_ode, t_range, y0);

figure;
plot(t, y(:, 1), 'r', t, y(:, 2), 'g', t, y(:, 3), 'b', t, y(:, 4), 'm', 'LineWidth',2);
legend('Z','A', 'A1', 'B');
xlabel('Time (t)');
ylabel('Mass by cell state');