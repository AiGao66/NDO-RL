function [dx,u,x1,x2,hat_x1,hat_x2,yr,z1,hat_alpha1,z2,d1,d2,hat_d1,hat_d2,hat_epsilon1,hat_epsilon2,rho1,rho2]=control_NDORL(t,x)
%% 参数
c1=10;c2=18;l1=15;l2=20;gamma_c1=15;gamma_c2=15;gamma_a1=16;gamma_a2=16;k1=10;k2=20;gamma1=2;gamma2=5;
Gamma1=0.2;Gamma2=0.5;theta1=20;theta2=30;p=99/101;q=101/99;
%% 变量
x1=x(1);x2=x(2);hat_x1=x(3);hat_x2=x(4);hat_epsilon1=x(5);hat_epsilon2=x(6);rho1=x(7);rho2=x(8);
% hat_Wf1=x(9:32);hat_Wf2=x(33:64);
% hat_Wc1=x(65:88);hat_Wc2=x(89:120);
% hat_Wa1=x(121:144);hat_Wa2=x(145:176);
hat_Wf1=x(9:32);hat_Wf2=x(33:64);
hat_Wc1=x(65:88);hat_Wc2=x(89:120);
hat_Wa1=x(121:144);hat_Wa2=x(145:176);
%% 参考信号
f2=0.5*sin(x1^2)-x2;
d1=0.001;
% d1=0.006*t;
d2=0.04*sin(0.04*t)+0.02;
yr=0.2*sin(t);
yr_d=0.2*cos(t);
%yr_dd=-1.8*cos(0.6*t);
%% step1
z1=x1-yr;
rule=linspace(-8,8,24);
S_f1=calS([x1,yr]',rule);
S_J1=calS([x1,yr]',rule);
hat_d1=hat_epsilon1+(rho1+gamma1)*hat_x1;
hat_alpha1=-c1*z1^(2*p-1)-l1*z1^(2*q-1)-hat_Wf1'*S_f1-(hat_Wa1'*S_J1)/2-hat_d1-11*z1/4;
hat_Wf1_d=Gamma1*S_f1*z1-theta1*hat_Wf1;
hat_Wc1_d=-gamma_c1*(S_J1'*S_J1)*hat_Wc1;
hat_Wa1_d=-(S_J1'*S_J1)*(gamma_a1*(hat_Wa1-hat_Wc1)+gamma_c1*hat_Wc1);

%% step2
z2=x2-hat_alpha1;
rule=linspace(-8,8,32);
S_f2=calS([x1,x2]',rule);
S_J2=calS([x1,x2,yr,yr_d]',rule);
hat_d2=hat_epsilon2+(rho2+gamma2)*hat_x2;
u=-c2*z2^(2*p-1)-l2*z2^(2*q-1)-hat_Wf2'*S_f2-(hat_Wa2'*S_J2)/2-hat_d2-9*z2/4;
hat_Wf2_d=Gamma2*S_f2*z2-theta2*hat_Wf2;
hat_Wc2_d=-gamma_c2*(S_J2'*S_J2)*hat_Wc2;
hat_Wa2_d=-(S_J2'*S_J2)*(gamma_a2*(hat_Wa2-hat_Wc2)+gamma_c2*hat_Wc2);
%% System model
x1_d=x2+d1;
x2_d=u+f2+d2;
%% Observer equation
hat_x1_d=hat_x2+hat_Wf1'*S_f1+k1*(x1-hat_x1)+hat_d1;
% hat_x1_d=hat_x2+k1*(x1-hat_x1)+hat_d1;
hat_x2_d=u+hat_Wf2'*S_f2+k2*(x1-hat_x1)+hat_d2;
% hat_x2_d=u+k2*(x1-hat_x1)+hat_d2;
rho1_d=-rho1^2-2*gamma1*rho1;
rho2_d=-rho2^2-2*gamma2*rho2;
hat_epsilon1_d=-rho1_d*hat_x1-(rho1+gamma1)*(hat_epsilon1+hat_x2+(rho1+gamma1)*hat_x1);
hat_epsilon2_d=-rho2_d*hat_x2-(rho2+gamma2)*(hat_epsilon2+u+(rho2+gamma2)*hat_x2);
%% All the dynamic equation
dx=[x1_d;x2_d;hat_x1_d;hat_x2_d;hat_epsilon1_d;hat_epsilon2_d;rho1_d;rho2_d;hat_Wf1_d;hat_Wf2_d;hat_Wc1_d;hat_Wc2_d;hat_Wa1_d;hat_Wa2_d];
%% NN
function S_Z=calS(Z,Rule)
t_i=repmat(Z,1,length(Rule))-repmat(Rule,length(Z),1); 
t=exp(-0.5*sum(t_i.*t_i,1))';
S_Z=t;
