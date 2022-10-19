%% Exercise 9.7 (same as Example17.8 in NR - used to get the transfer function coefficients)
% script <a href="matlab:Example_17_8">Example_17_8</a>
% Compute the ODE coefficients in the horizontal dynamics of a three-story building.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Exampe 17.8.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap17">Chapter 17</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.

clear; clf; format compact; close all;
m=1000; p1=15407; p3=10507; p5=5607; d1=1000; d3=1000; d5=1000; ell=5;
syms L1 L2 L3 L4 L5 L6 L7 L8 c k D
kbar1=vpa(k-2*p1/ell+sqrt(2)*d1/ell) 
kbar2=vpa(k-2*p3/ell+sqrt(2)*d3/ell)
kbar3=vpa(k-2*p5/ell+sqrt(2)*d5/ell)
L1=m*D^2+2*c*D+(kbar1+kbar2);
L2=c*D+kbar1;
L3=c*D+kbar2;
L4=m*D^2+2*c*D+(kbar2+kbar3);
L5=c*D+kbar2;
L6=c*D+kbar3;
L7=m*D^2+c*D+kbar3;
L8=c*D+kbar3;
a   =coeffs(expand(L1*L4*L7-L5*L3*L7-L8*L1*L6)/m^3,D); a   =a   (end:-1:1);
b   =coeffs(expand(L1*L4-L5*L3)/m^3,D);                b   =b   (end:-1:1);
bbar=coeffs(expand(L8*L5*L2)/m^3,D);                   bbar=bbar(end:-1:1); disp(' ')

% Substitute in numerical values for k & c, and plot the impulse response of the structure
for w=1%:3
  disp(sprintf('case %d:',w)); 
  switch w
    case 1, k=10000, c=10,   T=500; disp('This case is highly oscillatory.')
    %case 2, k=10000, c=1000, T=25;  disp('This case is still oscillatory, but not as bad.')
    %case 3, k=1000,  c=10,   T=10;  disp('This case is unstable!')
  end
  a_num=double(eval(a)), b_num=double(eval(b)), bbar_num=double(eval(bbar)),
  p=roots((a_num))
  for i=1:6
    d(i)=bbar_num(4)*p(i)^3+bbar_num(3)*p(i)^2+bbar_num(2)*p(i)+bbar_num(1);
    for j=1:6, if (j~=i), d(i)=d(i)/(p(i)-p(j)); end, end
  end
  d, h=T/5000; t=[0:h:T];
  f=(d(1)*exp(p(1)*t)+d(2)*exp(p(2)*t)+d(3)*exp(p(3)*t)+d(4)*exp(p(4)*t)+d(5)*exp(p(5)*t)+d(6)*exp(p(6)*t));
  plot(t,real(f)), print('-depsc',sprintf('mae143a_hw2_plot%d.eps',w))
  disp(' ')%, if w<3, pause, end
end

% end script Example_17_8


num_a = bbar_num;
num_b = b_num;
den = a_num;

figure (1)
bode(tf(num_a,den))
figure (2)
bode(tf(num_b,den))

roots(den)
%% Exercise 9.8
% b)
clear; clc; close all;

syms M K C m k c s

L6 = expand((K + C*s)*(m*s^2 + c*s + k)) %numerator
L5 = expand((M*s^2 + C*s + K)*(m*s^2 + c*s + k) - (k + c*s)^2 + (k + c*s)*(m*s^2 + c*s + k)) %denominator

% c)
L6_c = subs(L6,m,0) %numerator
L5_c = subs(L5,m,0) %denominator

% d)
L6_d = subs(L6_c,[M K C],[10 10 0]); %numerator
L5_d = subs(L5_c,[M K C],[10 10 0]); %denominator
poles = solve(L5_d,s)
zeros = solve(L6_d,s)

num = flip(coeffs(L6_d,s)); %flip bc expecting coefs in descending order
den = flip(coeffs(L5_d,s));
TF = RR_tf(num,den);

figure (1)
bode(tf(double(TF.num.poly),double(TF.den.poly)))

%e
figure (2)
impulse(tf(double(TF.num.poly),double(TF.den.poly)))
xlim([0 10])
figure (3)
step(tf(double(TF.num.poly),double(TF.den.poly)))
xlim([1 10])

%f
L6_f = subs(L6_c,[M K C],[10 10 14.14]); %numerator
L5_f = subs(L5_c,[M K C],[10 10 14.14]); %denominator
poles_f = solve(L5_f,s)
zeros_F = solve(L6_f,s)

num_f = flip(coeffs(L6_f,s)); %flip bc expecting coefs in descending order
den_f = flip(coeffs(L5_f,s));
TF_f = RR_tf(num_f,den_f);

figure (4)
bode(tf(double(TF_f.num.poly),double(TF_f.den.poly)))

figure (5)
impulse(tf(double(TF_f.num.poly),double(TF_f.den.poly)))
figure (6)
step(tf(double(TF_f.num.poly),double(TF_f.den.poly)))

%h (incomplete, my thought was to resolve the transfer function, multiply
%by 1/s and recompute the poles, but as we discussed in class, this is
%impossible for a 5th order system.
case1 = [-.130 + 1.002i;-.130 - 1.002i;-.145+.873i;-.145-.873i];
case2 = [-.257+.890i;-.257-.890i;-.268+.628i;-.268-.628i];

[TF_h_num, TF_h_den] = zp2tf(1,case1,1)
TF_h_den = [TF_h_den 0] %multiplies the 1/s term

test = RR_tf(TF_h_num,TF_h_den)

poles_h = solve(poly2sym(TF_h_den),x)
%% Exercise 9.12
clear; clc; close all;
num1 = [1];
num2 = [1 -1];
num3 = [1 -1];
num4 = [1 101 100];
num5 = [1 0];
num6 = [1];
den1 = [1 1];
den2 = [1];
den3 = [1 1];
den4 = [1 0];
den5 = [1 101 100];
den6 = [1 1 1];

figure (1)
bode(tf(num1,den1))
figure (2)
bode(tf(num2,den2))
figure (3)
bode(tf(num3,den3))
figure (4)
bode(tf(num4,den4))
figure (5)
bode(tf(num5,den5))
figure (6)
bode(tf(num6,den6))

%create blank plots for hand sketching
figure (7)
bode(tf(num1,den1),'w')
figure (8)
bode(tf(num2,den2),'w')
figure (9)
bode(tf(num3,den3),'w')
figure (10)
bode(tf(num4,den4),'w')
figure (11)
bode(tf(num5,den5),'w')
figure (12)
bode(tf(num6,den6),'w')