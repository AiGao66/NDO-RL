clear all;
clc;
close all;
x1_0=0.2;x2_0=0.2;hat_x1_0=0.1;hat_x2_0=0.2;hat_epsilon1_0=0;hat_epsilon2_0=0;rho1_0=0.8;rho2_0=1;
hat_Wf1_0=0*ones(24,1);hat_Wf2_0=0*ones(32,1);
hat_Wc1_0=0.2*ones(24,1);hat_Wc2_0=0.2*ones(32,1);
hat_Wa1_0=0.4*ones(24,1);hat_Wa2_0=0.4*ones(32,1);
x0=[x1_0;x2_0;hat_x1_0;hat_x2_0;hat_epsilon1_0;hat_epsilon2_0;rho1_0;rho2_0;hat_Wf1_0;hat_Wf2_0;hat_Wc1_0;hat_Wc2_0;hat_Wa1_0;hat_Wa2_0];
T=[0,20];
[t,x]=ode45('model_NDORL',T,x0,[]);

L=length(t);
for i=1:L
  [dx,u(i),x1(i),x2(i),hat_x1(i),hat_x2(i),yr(i),z1(i),hat_alpha1(i),z2(i),d1(i),d2(i),hat_d1(i),hat_d2(i),hat_epsilon1(i),hat_epsilon2(i),rho1(i),rho2(i)]=control_NDORL(t(i),x(i,:)');
end
%% 跟踪误差
figure;
plot(t,x1,t,yr,'--r','linewidth',2);
hold on;
legend('$\xi_{1}$','$y_{r}$','interpreter','latex');
xlabel('Time(sec)');
%% 观测误差x1,x2
figure;
subplot(2,1,1);
plot(t,x1-hat_x1,'-r','linewidth',2);
ylim([-0.5 0.5])
hold on;
legend('$\tilde{\xi}_{1}$','interpreter','latex');
xlabel('Time(sec)');
subplot(2,1,2);
plot(t,x2-hat_x2,'-b','linewidth',2);
ylim([-5 5])
hold on;
legend('$\tilde{\xi}_{2}$','interpreter','latex');
xlabel('Time(sec)');
%% hat_d1,hat_d2
figure;
subplot(2,1,1);
plot(t,d1-hat_d1,'-r','linewidth',2);
ylim([-1 1])
hold on;
legend('$\tilde{d}_{1}$','interpreter','latex');
xlabel('Time(sec)');
subplot(2,1,2);
plot(t,d2-hat_d2,'-r','linewidth',2);
ylim([-2 2])
hold on;
legend('$\tilde{d}_{2}$','interpreter','latex');
xlabel('Time(sec)');
%% rho1,rho2
figure;
plot(t,rho1,t,rho2,'-r','linewidth',2);
hold on;
legend('$\rho_{1}$','$\rho_{2}$','interpreter','latex');
xlabel('Time(sec)');
%% hat_Wf1,hat_Wf2
figure;
subplot(2,1,1);
hat_Wf1=sqrt(sum(x(:,9:32).^2,2));
plot(t,hat_Wf1,'-r','linewidth',2);
hold on;
legend('$\Vert\hat{\Theta}_{f1}\Vert$','interpreter','latex');
xlabel('Time(sec)');
subplot(2,1,2);
hat_Wf2=sqrt(sum(x(:,33:64).^2,2));
plot(t,hat_Wf2,'-b','linewidth',2);
hold on;
legend('$\Vert\hat{\Theta}_{f2}\Vert$','interpreter','latex');
xlabel('Time(sec)');
%% hat_Wc1,hat_Wc2
figure;
subplot(2,1,1);
hat_Wc1=sqrt(sum(x(:,65:88).^2,2));
plot(t,hat_Wc1,'-r','linewidth',2);
hold on;
legend('$\Vert\hat{\Theta}_{c1}\Vert$','interpreter','latex');
xlabel('Time(sec)');
subplot(2,1,2);
hat_Wc2=sqrt(sum(x(:,89:120).^2,2));
plot(t,hat_Wc2,'-b','linewidth',2);
hold on;
legend('$\Vert\hat{\Theta}_{c2}\Vert$','interpreter','latex');
xlabel('Time(sec)');
%% hat_Wa1,hat_Wa2
figure;
subplot(2,1,1)
hat_Wa1=sqrt(sum(x(:,121:144).^2,2));
plot(t,hat_Wa1,'-r','linewidth',2);
hold on;
legend('$\Vert\hat{\Theta}_{a1}\Vert$','interpreter','latex');
xlabel('Time(sec)');
subplot(2,1,2)
hat_Wa2=sqrt(sum(x(:,144:176).^2,2));
plot(t,hat_Wa2,'-b','linewidth',2);
hold on;
legend('$\Vert\hat{\Theta}_{a2}\Vert$','interpreter','latex');
xlabel('Time(sec)');
%% u
figure;
plot(t,u,'b','linewidth',2);
hold on;
xlabel('Time(sec)')
legend('$u$','interpreter','latex');
 %% 性能指标函数h2
figure;
subplot(2,1,1)
h1=z1.*z1+hat_alpha1.*hat_alpha1;
plot(t,h1,'linewidth',2);
hold on;
xlabel('Time(sec)')
legend('$h_{1}=z_{1}^2+\alpha_{1}^{2}$','interpreter','latex');
subplot(2,1,2)
h2=z2.*z2+u.*u;
plot(t,h2,'linewidth',2);
xlabel('Time(sec)')
legend('$h_{2}=z_{2}^2+u^{2}$','interpreter','latex');

 %% rho和u
figure;
subplot(2,1,1)
plot(t,rho1,t,rho2,'-r','linewidth',2);
hold on;
legend('$\rho_{1}$','$\rho_{2}$','interpreter','latex');
xlabel('Time(sec)');
subplot(2,1,2)
plot(t,u,'b','linewidth',2);
hold on;
xlabel('Time(sec)')
legend('$u$','interpreter','latex');
