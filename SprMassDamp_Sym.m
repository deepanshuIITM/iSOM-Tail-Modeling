% This function does Symbolic Calculation for the amplitude ratio of the system 
% with 1,2 and 3-dof (and it's variants). This function can be used for
% calculating amplitude ratio of the system in tuned mass damper problem

function [Y_func] = SprMassDamp_Sym(dof,spr1,damp1,spr2,damp2,spr3,damp3)

%% INPUTS:
% dof:      degrees of freedom of spring mass system(1,2,3)
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
% xx:       Input data in terms of r1,r2,r3

%% OUTPUT: 
% Y_func:   Amplitude ratio of the System in terms of MATLAB Function. 
%           Use this function by inputting values of J (damping ratio), R 
%           (Mass Ratio) & random variables (r1,r2,r3) 

%% Symbolic calculation
syms M w r1 r2 r3 J R real;
w1 = r1*w; w2 = r2*w; w3 = r3*w;
a=1;b=1;c=1; % Change value of a,b,c to adjust mass and damping
m1 = M; m2=R*M; m3=a*R*M;
k1= m1*w1^2; k2= m2*w2^2; k3= m3*w3^2;
c1 = 2*J*m1*w1; c2 = 2*J*m2*w2*b; c3 = 2*J*m3*w3*c;        

%% Conditions
if dof ==1
    m2=0;m3=0;damp2=0;damp3=0;spr2=0; spr3=0; 
    if nargin~=3
        disp('Invalid Input values! Check again...');Y_func = 'not valid'; return
    end
end
if dof==2
    m3=0; damp3=0;spr3=0; 
    if nargin~=5
        disp('Invalid Input values! Check again...');Y_func = 'not valid'; return
    end
end

if (spr1 ==0 ||spr1 <damp1 ||spr2 <damp2 ||spr3 <damp3)
    disp('Invalid Input values! Check again...');Y_func = 'not valid'; return
end

if damp1 ==0; c1 = 0; end
if damp2 ==0; c2 = 0; end
if damp3 ==0; c3 = 0; end
if spr2 ==0; k2 = 0; end
if spr3 ==0; k3 = 0; end

%% Matrix
A = [(k1+k2-m1*w^2+(c1+c2)*w*1i),      -(k2+c2*w*1i),                   0;
         -(k2+c2*w*1i),          (k2+k3-m2*w^2+w*(c2+c3)*1i),     -(k3+c3*w*1i);
                0,                     -(k3+c3*w*1i),          (k3-m3*w^2+c3*w*1i)];
f = [1;0;0];

if dof ==1
    A = A(1,1); f = f(1);
end
if dof ==2
    A = A(1:2,1:2); f = f(1:2);
end

eig_A = matlabFunction(eig(simplify(A/(m1*w1^2))));
%% Solution
X = A\f;
Y_func = matlabFunction(abs(simplify(X(1)*m1*w1^2)));

if (dof ==1 && nargin(Y_func)>3)
    disp('Invalid Input values! Check again...');
    Y_func = 'not valid';
end
if (dof ==2 && nargin(Y_func)>4)
    disp('Invalid Input values! Check again...');
    Y_func = 'not valid';
end
if (dof ==3 && nargin(Y_func)>5)
    disp('Invalid Input values! Check again...');
     Y_func = 'not valid';
end

end