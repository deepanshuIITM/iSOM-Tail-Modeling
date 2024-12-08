% This function does Symbolic Calculation for the amplitude of the system 
% with 1,2 and 3-dof (and it's variants). This function can be used for
% calculating amplitude of the system in tuned mass damper system
function [y, Data] = TMD_3DOF(n,absorber,damp2,damp3,spr2,spr3,xx)

%% INPUTS:
% n:        no of response you want to generate
% absorber: 1 for 2-Dof and 2 for 3-Dof
% spr1:     it takes logical values 0 & 1; if spr1= 0; no spring in
%           system; if spr1 = 1; spring is present in system
% spr2:     it takes logical values 0 & 1; if spr2= 0; no spring in
%           absorber 1; if spr2 = 1; spring is present in absorber 1
% spr3:     it takes logical values 0 & 1; if spr3= 0; no spring in
%           absorber 2; if spr2 = 1; spring is present in absorber 2
% damp1:    it takes logical values 0 & 1; if damp1= 0; no damper in
%           system; if damp2 = 1; damper is present in system
% damp2:    it takes logical values 0 & 1; if damp2= 0; no damper in
%           absorber 1; if damp2 = 1; damper is present in absorber 1
% damp3:    it takes logical values 0 & 1; if damp3= 0; no damper in
%           absorber 2; if damp3 = 1; damper is present in absorber 2
% xx:       Input in terms of r1, r2 and/ r3 (This is optional. Use only if you want to calculate
%           response on given input values)

%% OUTPUT: 
% Y_func:        Amplitude of the System in terms of MATLAB Function. Use this
% function and input values of J (damping ratio), R (Mass Ratio) & random 
% variables (r1,r2,r3) 

%% Symbolic calculation
syms M w r1 r2 r3 J R real;
w1 = r1*w; w2 = r2*w; w3 = r3*w;
a=1;b=1;c=1; % Change value of a,b,c to adjust mass and damping
m1 = M; m2=R*M; m3=a*R*M;
k1= m1*w1^2; k2= m2*w2^2; k3= m3*w3^2;
c1 = 2*J*m1*w1; c2 = 2*J*m2*w2*b; c3 = 2*J*m3*w3*c;        
% c3 = 0; k3 = 0;m3=0;c2=0;r3=0;

%% Conditions
if absorber ==1
    m3 = 0;
end
if damp2 ==0
    c2 = 0; 
end
if damp3 ==0
    c3 = 0; 
end
if spr2 ==0
    k2 = 0; 
end
if spr3 ==0
    k3 = 0; 
end
    
%% Matrix
A = [(k1+k2-m1*w^2+(c1+c2)*w*1i),      -(k2+c2*w*1i),                   0;
         -(k2+c2*w*1i),          (k2+k3-m2*w^2+w*(c2+c3)*1i),     -(k3+c3*w*1i);
                0,                     -(k3+c3*w*1i),          (k3-m3*w^2+c3*w*1i)];
f = [1;0;0];

if absorber ==1
    A = A(1:2,1:2); f = f(1:2);
end
X = A\f;
sol = X(1)*m1*w1^2;

%% Symbol to function
J=0.01; R=0.025; %0.025
y_func = matlabFunction(sol);

if absorber ==1
    mu = 1; sig = 0.025; % Change the value of mu and sigma here 0.025
    if (nargin == 6)
        r1=normrnd(mu,sig,n,1); r2=normrnd(mu,sig,n,1);
    elseif (nargin == 7)
        r1 = xx(:,1); r2 = xx(:,2);
    end 
    y = abs(y_func(J,R,r1,r2))/27;
    Data = [r1 r2 y];
else
    J=0.01; R=0.025;
    mu = 1; sig = 0.05;  % Change the value of mu and sigma here (0.05)
    if (nargin == 6)
        r1=normrnd(mu,sig,n,1); r2=normrnd(mu,sig,n,1);r3=normrnd(mu,sig,n,1);
    elseif (nargin == 7)
        r1 = xx(:,1); r2 = xx(:,2); r3= xx(:,3);
    end 
    y = abs(y_func(J,R,r1,r2,r3))/27;
    Data = [r1 r2 r3 y]; 
end
end