y0 = [0; 2000; 2000; 0];
t_range = [0 50];

[t, y] = ode45(@case2_ode, t_range, y0);

figure;
plot(t, y(:, 1), 'r', t, y(:, 2), 'g', 'LineWidth',2);
legend('A1','B');
xlabel('Time (t)');
ylabel('Mass by cell state');