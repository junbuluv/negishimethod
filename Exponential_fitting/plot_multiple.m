

figure()
subplot(2,2,1)
hold on 
plot(c1_psi1)
plot(c2_psi1)
hold off
title("$\psi = 0.2$","Interpreter","latex")
xlabel("$t$","Interpreter","latex")
ylabel("$c$","Interpreter","latex")
subplot(2,2,2)
hold on
plot(c1_psi2)
plot(c2_psi2)
hold off
title("$\psi = 0.5$","Interpreter","latex")
xlabel("$t$","Interpreter","latex")
ylabel("$c$","Interpreter","latex")
subplot(2,2,3)
hold on
plot(c1_psi3)
plot(c2_psi3)
xlabel("$t$","Interpreter","latex")
ylabel("$c$","Interpreter","latex")
hold off
title("$\psi = 0.7$","Interpreter","latex")
subplot(2,2,4)
hold on
plot(c1_psi4)
plot(c2_psi4)
hold off
title("$\psi = 0.9$","Interpreter","latex")
xlabel("$t$","Interpreter","latex")
ylabel("$c$","Interpreter","latex")


fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('Agent 1','Agent 2');
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.5;
