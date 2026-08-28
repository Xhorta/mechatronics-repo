%six bar linkage
%Static Equilibrium 

clc; 
clear; 

%define joints 

A=[7 4 0];
B=[5 16 0];
C=[25 25 0];
D=[23 10 0];
E=[18 35 0];
F=[43 32 0];
G=[45 17 0];

%Define the lengths of the links 
lAB = norm(B-A);
lBC = norm(C-B);
lBE = norm(E-B);
lCD = norm(D-C);
lEF = norm(F-E);
LFG = norm(G- F);

%Weight of Each link 
WAB = [0 -1 0];
WBEC =[0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];

%Center of Mass for each link

S1=(A+B)/2;
S2=(B+E+C)/3; % not precise 
S4=(E+F)/2;
S3=(C+D)/2;
S5=(F+G)/2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA= [FAx FAy 0];
ForceB= [FBx FBy 0];
ForceC= [FCx FCy 0];
ForceD= [FDx FDy 0];
ForceE= [FEx FEy 0];
ForceF= [FFx FFy 0]; 
ForceG= [FGx FGy 0];
InputTorque = [0 0 Tin];


%Applied Force 
AppliedForce = [50 0 0];

%Static Equilibrium Conditions for Link AB 

% Sum of Forces= 0 
%Fa + Fb + WAB = 0
eqn1= ForceA + ForceB + WAB == 0; 

%Sum of moments = 0 with respect to CoM of link AB
% S1A x FA + S1B x FB + InputTorque= 0 

eqn2= cross(A-S1,ForceA)+ cross(B-S1,ForceB)+ InputTorque == 0;

% Link BEC 
%Sum of forces = 0 
% -Fb + Fc + Fe + WBEC = 0 
eqn3= -ForceB + ForceC + ForceE + WBEC ==0; 

%sum of moments = 0 
% Sum of moments = 0 with respect to CoM of link BEC
% S2B x -FB + S2C x FC + S2E x FE  = 0 
eqn4 = cross(S2-B, -ForceB) + cross(S2-C, ForceC) + cross(S2-E, ForceE) == 0;


%link  CD  
%Sum of forces = 0 
% -Fc + Fd + WCD = 0 
eqn5 = -ForceC + ForceD + WCD == 0; 

%sum of moments=0
% Link CD
% Sum of moments = 0 with respect to CoM of link CD
% S3C x -FC + S3D x FD  = 0
eqn6 = cross(S3-C, -ForceC) + cross(S3-D, ForceD) == 0;

%link  EF   
%Sum of forces = 0 
% -EF + FF + WFG = 0 
eqn7 = -ForceE + ForceF + WEF == 0; 

% Sum of moments = 0 with respect to CoM of link FG
% S4E x -FE + S4F x FF = 0
eqn8= cross(E-S4, -ForceE)+ cross(F-S4,ForceF) == 0; 



% Link FG 
% Sum of Force  = 0
eqn9 = -ForceF + ForceG + WFG + AppliedForce== 0;

%Sum  of Moments with respect to CoM of link FG
% S5F x -FF + S5G x FG = 0
eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == 0;

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution=solve(eqnMatrix,[FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);



