clear all
close all
clc

% Global constants

global cw ca A K mpunkt_quer F_quer tc alpha lrampe r0

% Basic Design Data of Rocket

m1=2.
m0=80
cw=0.4
ca=0.
alpha=0.
lrampe=10.
D=0.2
A=D^2*pi/4

% Genereric Constants

r0=6371000
rho=1.225
K=3.9658e14
g0=9.81

% Motor Data

m8=26.
Is=1700.
tc=15.
mpunkt_quer=m8/tc
Itot=m8*Is
F_quer=Itot/tc

% Calculated Rocket Data

mc=m0-m8
R=m0/mc

% Numerical solution
gamma=80*pi/180

[T,Y]=ode15s(@Rocket_2DOF,[0 101],[0 r0 m0 gamma 0]);

% Calculate accelerations as derivatives from velocities
    deltav=diff(Y(:,1));
    deltat=diff(T);
    lt=length(T)
    for j=1 : lt-1
        acc(j)=deltav(j)/deltat(j);
    end
    acc(lt)=acc(lt-1)
    if length (acc) > lt
        la=length(acc);
        acc(lt+1:la)=[];
    end
% plot results   
figure (1)
plot (T,acc,'b-')
ylabel ('Acceleration [m/s^2]')
xlabel ('Time [s]')
figure (2)
plot (T,Y(:,1),'r-')
ylabel ('Velocity [m/s]')
xlabel ('Time [s]')
figure (3)
plot (T,Y(:,2)-r0,'b-')
ylabel ('Altitude [m]')
xlabel ('Time [s]')
figure (4)
plot (T,Y(:,3))
ylabel ('Mass [kg]')
xlabel ('Time [s]')
figure (5)
plot (T,Y(:,4)*180/pi)
ylabel ('gamma [°]')
xlabel ('Time [s]')


