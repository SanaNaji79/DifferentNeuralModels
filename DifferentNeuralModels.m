
%% README
%in this code the different nueral models are investigated and some fundamental concepts are shown they can be
%listed as below:
% 1) hodjkin-huxley model
% 2) izhikevich model
% 3) noisy output model
% 4) bifurcation investigation
% 5) phase plane investigation
% 6) Simple two-neuron network modeling

%% problem1 (hodjkin-huxley model)

dt = 0.01; % Simulation time step
Duration = 20000; % Simulation length 
T = ceil(Duration/dt); 
vRest1 = -60 ;
t = (1:T) * dt; % Simulation time points in ms 
Cm = 1; % Membrane capacitance in micro Farads 
gNa = 120; % in Siemens, maximum conductivity of Na+ Channel 
gK = 36; % in Siemens, maximum conductivity of K+ Channel 
gl = 0.3; % in Siemens, conductivity of leak Channel 
ENa = 55; % in mv, Na+ nernst potential 
EK = -72; % in mv, K+ nernst potential 
El = -49.4; % in mv, nernst potential for leak channel 
vRest = -60; % in mv, resting potential 
V = vRest1 * ones(1,T); % Vector of output voltage 
I = zeros(1,T); % in uA, external stimulus (external current) 
 % for example: I(1:10000) = 2; % an input current pulse


a = linspace(-80 , 120 , 20000) ;
u = vRest - a;
alpha_n = (.1 * u + 1)./(exp(1 + .1 * u) - 1) / 10;
beta_n = .125 * exp(u/80);
alpha_m = (u+25) ./ (exp(2.5+.1*u)-1)/10;
beta_m = 4*exp(u/18);
alpha_h = .07 * exp(u/20);
beta_h = 1 ./ (1+exp(3 + .1*u));

% part one (time constant plots)
ta_n = 1./(alpha_n + beta_n) ;
ta_m = 1./(alpha_m + beta_m) ;
ta_h = 1./(alpha_h + beta_h) ;
plot(a , ta_h , a , ta_m , a , ta_n , 'LineWidth' , 1.5) ;
grid on ;
title('Time Constants') ;
legend('\tau_h' , '\tau_m' , '\tau_n') ;
legend('Location','NorthEast') ;

n_inf = alpha_n./(alpha_n + beta_n) ;
m_inf = alpha_m./(alpha_m + beta_m) ;
h_inf = alpha_h./(alpha_h + beta_h) ;
plot(a , h_inf , a , m_inf , a , n_inf , 'LineWidth' , 1.5) ;
grid on ;
title('Steady-State Values') ;
legend('h_\infty' , 'm_\infty' , 'n_\infty') ;
legend('Location','NorthEast') ;

% part 2 (the membrane's voltage)

% parameters declaration
n = n_inf(vRest1+80)*ones(1,T);
m = m_inf(vRest1+80)*ones(1,T);
h = h_inf(vRest1+80)*ones(1,T);
% I = [40*ones(1 , 5000) , zeros(1 , T-5000)] ;
% I = 5.80005*[zeros(1 , 10000) , ones(1 , T-10000) ] ;
%I = 2.171*ones(1 , T) ;
%I(1000000:1010000) = 10 ;
% k = 10000 ;
% I(k:k+122) = 6 ;
% d = 200 ;
% for i = 0:1000
%    I(10000+i*d+1:10000+i*d+d) = i/10 ;
% end
for i = 1:T 
    I(i) = 10*sin(5*t(i)) ;
end
a = T;
for i = 1:a-1
    V(i+1) = V(i) - (dt/Cm)*(gl*(V(i)-El)+gK*((n(i)).^4)*(V(i)-EK)+gNa*(m(i).^3)*h(i)*(V(i)-ENa)-I(i));
    n(i+1) = n(i) + (dt/ta_n(100*ceil(V(i)+80)))*(-n(i)+n_inf(100*ceil(V(i)+80))) ;
    m(i+1) = m(i) + (dt/ta_m(100*ceil(V(i)+80)))*(-m(i)+m_inf(100*ceil(V(i)+80))) ;
    h(i+1) = h(i) + (dt/ta_h(100*ceil(V(i)+80)))*(-h(i)+h_inf(100*ceil(V(i)+80))) ;
end
plot( t(1:a/10),V(1:a/10)) ;
xlabel('time(ms)') ;
ylabel('voltage(mV)') ;

%% problem2 (izhikevich model)
% part 3
clc 
clear
dt = 0.01; % Simulation time step
Duration = 100; % Simulation length 
T = ceil(Duration/dt); 
t = (1:T) * dt; % Simulation time points in ms 
a = 0.02 ;
b = 0.20 ;
c = -50 ;
d = 2 ;
h = 15 ;
v0 = c ;
u0 = -14 ;
I = h*[zeros(1 , 10) , ones(1 , T-10)] ;
V = v0*ones(1 , T) ;
u = u0*ones(1 , T) ;

for i = 1:T-1
    V(i+1) = V(i) + dt*(0.04*(V(i).^2)+5*V(i)+140-u(i)+I(i)) ;
    u(i+1) = u(i) + dt*a*(b*V(i)-u(i)) ;

    if(V(i) >= 30)
        V(i+1) = c ;
        u(i+1) = u(i) + d ;
   end 

end

plot(t , V) ;
xlabel('time(ms)') ;
ylabel('voltage(mV)') ;
%% problem3 (noisy output model)
clc 
clear

beta =  0.5 ;
gama = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] ;
syms F1(x) F2(x) F3(x) F4(x) F5(x) F6(x) F7(x) F8(x) F9(x)
syms f1(x) f2(x) f3(x) f4(x) f5(x) f6(x) f7(x) f8(x) f9(x)
A = [F1(x) , F2(x) , F3(x) , F4(x) , F5(x) , F6(x) , F7(x) , F8(x) , F9(x)] ;
B = [f1(x) , f2(x) , f3(x) , f4(x) , f5(x) , f6(x) , f7(x) , f8(x) , f9(x)] ;
for i = 1:9
    A(i) = beta*(1 + tanh(gama(i)*x)) ;
 
end
for i = 1:9
    hold on ;
    fplot(A(i) , [-50 , 50]) ;
end
title('CDF')
legend('\gamma=0.1' , '\gamma=0.2' , '\gamma=0.3' , '\gamma=0.4' , '\gamma=0.5' , '\gamma=0.6' , '\gamma=0.7' , '\gamma=0.8' , '\gamma=0.9') ;
hold off ;
for i = 1:9
    B(i) = beta*gama(i)*(4/(exp(gama(i)*x)+exp(-gama(i)*x))^2) ; ;
 
end
for i = 1:9
    hold on ;
    fplot(B(i) , [-50 , 50]) ;
end
title('PDF')
legend('\gamma=0.1' , '\gamma=0.2' , '\gamma=0.3' , '\gamma=0.4' , '\gamma=0.5' , '\gamma=0.6' , '\gamma=0.7' , '\gamma=0.8' , '\gamma=0.9') ;
hold off ;

% part 2
RI = 30 ;
Vth = -45 ;
Vrest = -70 ;
Vpeak = 30 ;
ta_m = 2 ;
dt = 0.001; % Simulation time step
Duration = 20; % Simulation length 
T = ceil(Duration/dt); 
t = (1:T) * dt; % Simulation time points in ms 

v = Vrest*ones(1 , T) ;
% for i = 1:T-1
%     if((v(i) >= Vth)&&(v(i) < Vpeak))
%         v(i+1) = Vpeak ;
%     elseif(v(i) == Vpeak)
%         v(i+1) = Vrest ;
%     else
%         v(i+1) = v(i) + (dt/ta_m)*(-v(i) + RI) ;
%     end  
% end
% plot(t , v) ;

f = 0 ;
k =60 ;
F = zeros(1 , k+1) ;
current = (0:60)*0.5 ;
for j = 0:k
    for i = 1:T-1
    if((v(i) >= Vth)&&(v(i) < Vpeak))
        v(i+1) = Vpeak ;
        f = f+1 ;
    elseif(v(i) == Vpeak)
        v(i+1) = Vrest ;
    else
        v(i+1) = v(i) + (dt/ta_m)*(-v(i) + j/2) ;
    end  
    end
    F(j+1) = 10*f ;
    f = 0 ;
end

plot(current , F) ;
xlabel('current') ;
ylabel('freguency') ;
title('firing rate') ;

%% problem5 (bifurcation investigation)
clc 
clear

% a = -9 ;
% [t , x] = meshgrid(-5:0.2:5 , -5:0.2:5) ;
% s = x.*x + a ;
% l = sqrt(1 + s.^2) ;
% quiver(t , x , 1./l , s./l , 0.45) ;
% 
syms y1(k) y2(k)

y1(k) = sqrt(-k) ;
y2(k) = -sqrt(-k) ;

fplot(k , y1 , '--') ;
hold on ;
plot(0 , 0 , '.') ;
hold on ;
fplot(k , y2) ;
xlabel ( 'parameter a') ;
ylabel ( 'fixed points') ;




%% problem6 (phase plane)

clc 
clear 
figure 
% initial parameters 
noip = 15; 
interval = 60; 
Mee = 1.25; 
Mei = -1; 
Mie = 1; 
Mii = 0; 
Ye = -10; 
Yi = 10; 
Te = 0.01; 
Ti = 0.05; 
%f = @(t,Y) [(-Y(1)+Mee*Y(1)+Mei*Y(2)-Ye)/Te;(-Y(2)+Mie*Y(1)+Mii*Y(2)-Yi)/Ti]; 
f = @(t,Y) [(-Y(1)+afunc(Mee*Y(1)+Mei*Y(2)-Ye))/Te;(-Y(2)+afunc(Mie*Y(1)+Mii*Y(2)-Yi))/Ti]; 
y1 = linspace(-interval,interval,20); 
y2 = linspace(-interval,interval,20); 
% creates two matrices one for all the x-values on the grid, and one for 
% all the y-values on the grid. Note that x and y are matrices of the same 
% size and shape, in this case 20 rows and 20 columns 
[x,y] = meshgrid(y1,y2); 
u = zeros(size(x)); 
V = zeros(size(x)); 
% we can use a single loop over each element to compute the derivatives at 
% each point (y1, y2) 
t=0; % we want the derivatives at each point at t=0, i.e. the starting time 
for i = 1:numel(x) 
Yprime = f(t,[x(i); y(i)]); 
u(i) = Yprime(1); 
V(i) = Yprime(2); 
end 
quiver(x,y,u,V,'r'); 
xlabel('V_E') 
ylabel('V_I') 
% axis tight equal; 
hold on 
for i = 1:noip 
[ts,ys] = ode45(f,[0,50],[rand()*interval*((-1)^floor(rand()*interval)); ... 
rand()*interval*((-1)^floor(rand()*interval))]); 
plot(ys(:,1),ys(:,2),'b') 
plot(ys(1,1),ys(1,2),'bo') % starting point 
plot(ys(end,1),ys(end,2),'ks') % ending point 
xlim([-interval interval]); 
ylim([-interval interval]); 
end 
syms t 
fplot((t+Ye-Mee*t)/Mei,'m','LineWidth',2); 
fplot((Yi-Mie*t)/(Mii-1),'g','LineWidth',2); 
title("phase plane") 
hold('off') 
%% simple two-neurons modeling
%% problem2
clear 
clc
ta_m = 20 ;
E_L = -70 ;
IR = 25 ;
r_m = 100 ;
V_th = -54 ;
ta_peak = 10 ;
E_s = -80 ;
v0 = -80 ;
% time 
dt = 0.01 ; % Simulation time step
Duration = 200 ; % Simulation length 
T = ceil(Duration/dt) ; 
t = (1:T) * dt; % Simulation time points in ms
f = 50 ; % the pulse frequency
% 
K = (((t(1))/ta_peak)*(exp(1 - (t(1))/ta_peak)))*ones(1 , T) ;
v = (v0)*ones(1 , T) ;
% forloop
for i = 1:T-1
    K(i+1) = (t(i+1)/ta_peak)*(exp(1 - t(i+1)/ta_peak)) ;
end
a = 0 ;
j = 0 ;
while f*(j) < T
      a = a + K(f*(j) + 1) ;
      j = j + 1 ;
end
g0 = a ;
g = (g0)*ones(1 , T) ;
for i = 1:T-1
    a = 0 ;
    j = 0 ;
    while f*(j) + mod(i+1 , f) + 1 < T
       a = a + K(f*(j) + mod(i+1 , f) + 1) ;
       j = j + 1 ;
    end
    g(i+1) = a ;
    v(i+1) = v(i) +(dt/ta_m)*(IR + -(v(i)-E_L+(g(i)*(v(i)-E_s)*r_m/1000))) ;
    if v(i) > V_th
        v(i+1) = 0 ;
    end
    if v(i) == 0
        v(i+1) = -80 ;
    end
end
figure ;
b = T ;
plot(t(1:b) , v(1:b)) ;
xlabel('time(ms)') ;
ylabel('potential(mV)') ;
title('neuron potential') ;
% figure ;
% c = T ;
% plot(t , K) ;
% xlabel('time(ms)') ;
% title('K(t)') ;
% figure ;
% plot(t(1:c) , g(1:c)) ;
% xlabel('time(ms)') ;
% title('g(t)') ;

%% problem2 part 2
clear 
clc
ta_m = 20 ;
E_L = -70 ;
IR = 25 ;
r_m = 100 ;
V_th = -54 ;
ta_peak = 10 ;
E_s1 = 0 ;
E_s2 = -80 ;
v_rest = -80 ;
v0_1 = -80 ;
v0_2 = -60 ;
% time 
dt = 0.01 ; % Simulation time step
Duration = 600 ; % Simulation length 
T = ceil(Duration/dt) ; 
t = (1:T) * dt; % Simulation time points in ms

K = (((t(1))/ta_peak)*(exp(1 - (t(1))/ta_peak)))*ones(1 , T) ;
v1 = (v0_1)*ones(1 , T) ;
v2 = (v0_2)*ones(1 , T) ;
g1 = ones(1 , T) ;
g2 = ones(1 , T) ;
spike1 = [] ;
spike2 = [] ;

for i = 1:T-1
    v1(i+1) = v1(i) +(dt/ta_m)*(IR + -(v1(i)-E_L + (g1(i)*(v1(i)- E_s1)*r_m/1000))) ;
    v2(i+1) = v2(i) +(dt/ta_m)*(IR + -(v2(i)-E_L + (g2(i)*(v2(i)- E_s2)*r_m/1000))) ;
    if v1(i) > V_th
        v1(i+1) = 0 ;
        spike1 = [spike1 , i+1] ;
    end
    if v1(i) == 0
        v1(i+1) = v_rest ;
    end
    if v2(i) > V_th
        v2(i+1) = 0 ;
        spike2 = [spike2 , i+1] ;
    end
    if v2(i) == 0
        v2(i+1) = v_rest ;
    end
    a = 0 ;
    for j = spike2
       a = a + K(1 + (i + 1 - j)) ; 
    end
    g1(i+1) = a ;
    
    a = 0 ;
    for j = spike1
       a = a + K(1 + (i + 1 - j)) ; 
    end
    g2(i+1) = a ;
end

figure ;
b = T ;
plot(t(1:b) , v1(1:b) , t(1:b) , v2(1:b)) ;
xlabel('time(ms)') ;
ylabel('potential(mV)') ;
title('neuron potential') ;
%plot(t(1:b) , g1(1:b)) ;


function x = afunc(y)
if(y>0)
    x = 1;
else
    x = 0 ;
end
end