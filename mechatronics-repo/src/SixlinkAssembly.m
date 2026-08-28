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


Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);
Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);
Input_Torque = double(StaticSolution.Tin);


%Display the results of the forces and input torque
disp('Force A:');
disp([Force_Ax , Force_Ay]);
disp('Force B:');
disp([Force_Bx , Force_By]);
disp('Force C:');
disp([Force_Cx, Force_Cy]);
disp('Force D:');
disp([Force_Dx , Force_Dy]);
disp('Force E:');
disp([Force_Ex, Force_Ey]);
disp('Force F:');
disp([Force_Fx, Force_Fy]);
disp('Force G:');
disp([Force_Gx ,Force_Gy]);
disp('Input Torque:'); 
disp(Input_Torque);

%Angular Velocity Calc
%Loop ABCDA

 syms wBEC wCD wEF wFG
omega_AB=[0 0 1];
omega_BEC=[0 0 wBEC];%assuming its in the counter clockwise position
omega_CD= [ 0 0 wCD];


eqn11= cross(omega_AB, B-A)+ cross(omega_BEC, C-B)+cross(omega_CD, D-C)==0;


%---------------------------------------------------------------------------------------
%solving unknown variables using the looping method 
%using the ABCDA, DCEFGD loops we were able to solve for the angular
%velociteis 
loop1Solution= solve(eqn11,[wBEC wCD]);

%Extract angular velocities from loop solution
angularVelocity_BEC = double(loop1Solution.wBEC);
angularVelocity_CD = double(loop1Solution.wCD);

%second loop 
%DCEFGD

omegaBEC= [ 0 0 angularVelocity_BEC];
omegaCD=[0 0 angularVelocity_CD];
omega_FG=[0 0 wFG];
omega_EF= [ 0 0 wEF];


eqn12= cross(omegaCD, C-D)+cross(omegaBEC,E-C) + cross(omega_EF,F-E)+cross(omega_FG,G-F)==0;
loop2Solution=solve(eqn12,[wEF wFG]);

%Extract angular velocities from loop solution
angularVelocity_EF = double(loop2Solution.wEF);
angularVelocity_FG = double(loop2Solution.wFG);

% this finds the angular acceleration through using 
% loop 1 ABCDA 


syms aBEC aCD aFG aEF 

alpha_AB= [0 0 0]; % assume there is no acceleration for AB
alpha_BEC=[0 0 aBEC];
alpha_CD=[0 0 aCD];

a_B_A= cross(alpha_AB, B-A)+ cross(omega_AB, cross(omega_AB, B-A));
a_C_B= cross(alpha_BEC, C-B)+ cross(omegaBEC, cross(omegaBEC, C-B));
a_D_C= cross(alpha_CD, D-C)+ cross(omegaCD, cross(omegaCD, D-C));

eqn13=a_D_C+ a_C_B+ a_B_A ==0; 

loop1AccSolution= solve(eqn13,[aBEC aCD]);
alphaBEC= double(loop1AccSolution.aBEC);
alphaCD= double(loop1AccSolution.aCD);

%loop 2
%DCEFGD

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];

alpha_FG=[0 0 aFG];
alpha_EF=[0 0 aEF];

angVel_FG= [0 0 angularVelocity_FG];
angVel_EF= [0 0 angularVelocity_EF];



%a_C_D + a_E_C+ a_F_E+ a_G_F= 0
a_C_D= cross(alphaCD_vector, C-D)+ cross(omegaCD, cross(omegaCD, C-D));

a_E_C= cross(alphaBEC_vector, E-C)+ cross(omegaBEC, cross(omegaBEC, E-C));

a_F_E= cross(alpha_EF, F-E)+ cross(angVel_EF, cross(angVel_EF, F-E));

a_G_F= cross(alpha_FG, G-F)+ cross(angVel_FG, cross(angVel_FG, G-F));

eqn14= a_E_C+ a_C_D+ a_G_F + a_F_E ==0;

loop2AccSolution= solve(eqn14,[aEF aFG]);
alphaEF= double(loop2AccSolution.aEF);
alphaFG= double(loop2AccSolution.aFG);

vB_A = cross(omega_AB,B-A);
vC_B = cross(omega_BEC, C-B);



v_E_B= cross(omegaBEC,E-B);

VE_A = v_E_B +vB_A;

vC_D= cross(omegaCD, C-D);


vF_G= cross(angVel_FG, F-G);


V_S4_F= cross(angVel_EF,S4-F);
vS4_G= V_S4_F + vF_G; 

vS1_A = cross(omega_AB, S1-A);

vS2_A = cross(omegaBEC, S2-B)+ vB_A;

vS3_D = cross(omegaCD, S3-D);

vS5_G = cross(angVel_FG, S5-G);



%Velocity at joint 
%VS1_A= cross(omega_AB, S1-A);
%VS2_A= cross(omegaBEC,S2-B)+vB_A;
%VS2_D= cross(omegaCD,S3-D);
%VS5_G= cross(angVel_FG,S5-G);

aB_A = cross(alpha_AB, B-A)+ cross(omega_AB,cross(omega_AB, B-A));

aC_D = cross(alphaCD_vector, C-D)+ cross(omegaCD,cross(omegaCD, C-D));

aE_D = cross(alphaBEC_vector, E-C)+ cross(omegaBEC,cross(omegaBEC, E-C))+ aC_D;

aF_G = cross([0 0 alphaFG], F-G)+ cross(angVel_FG,cross(angVel_FG, F-G));



aS1_A= cross(alpha_AB,S1- A) + cross(omega_AB, cross(omega_AB,S1-A));
aS2_A= cross(alphaBEC_vector,S2-B)+ cross(omegaBEC,cross(omegaBEC,S2-B)) + aB_A;
aS3_D= cross(alphaCD_vector,S3-D)+ cross(omegaCD,cross(omegaCD,S3-D));
aS4_G= cross([ 0 0 alphaEF], S4-F)+ cross(angVel_EF,cross(angVel_EF,S4-F)) + aF_G;
aS5_G= cross([ 0 0 alphaFG], S5-G)+ cross(angVel_FG,cross(angVel_FG,S5-G));


%Newtons Second law Implementation 

MassAB= 1;
MassBEC= 1; 
MassCD = 1;
MassEF = 1; 
MassFG = 1;

%mass moment of interia
%obtain through solid works through mass properties 
J_AB = 1;
J_BEC = 1; 
J_CD = 1; 
J_EF = 1; 
J_FG = 1; 

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

NForceA = [NFAx NFAy 0];
NForceB= [NFBx NFBy 0];
NForceC= [NFCx NFCy 0];
NForceD= [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF= [ NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NForceTorque= [0 0 NTin];

%Eqn for link AB 
%Sun of forces 
eqn15= NForceA+ NForceB+ WAB == MassAB *aS1_A;
eqn16 = cross(A-S1,NForceA)+ cross(B-S1,NForceB)+ NForceTorque== J_AB* alpha_AB;


%eqn for link BEC
%Sum of forces 
eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBEC * aS2_A;

eqn18= cross(B-S2, -NForceB)+ cross(C-S2,NForceC)+ cross(E-S2, NForceE) == J_BEC * alphaBEC_vector;


%Eqn for link CD 
%sum of forces 

eqn19= - NForceC+ NForceD + WCD == MassCD *aS3_D;

eqn20 = cross(C-S3, -NForceC )+ cross(D-S3, NForceD )== J_CD * alphaCD_vector;

%eqn for link Ef

eqn21= -NForceE + NForceF+ WEF == MassEF * aS4_G;


eqn22 = cross(E-S4, -NForceE )+ cross(F-S4, NForceF)== J_EF *[0 0 alphaEF];


eqn23= -NForceF + NForceG + WFG + AppliedForce == MassFG *aS5_G;


eqn24= cross(F-S5, -NForceF)+ cross(G-S5,NForceG) == J_FG * [ 0 0 alphaFG];

NeqMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19,eqn20, eqn21, eqn22 , eqn23 ,eqn24];

DynamicSolution= solve(NeqMatrix, [NFAx, NFAy, NFBx, NFBy, NFCx,NFCy,NFDx,NFDy,NFEx,NFEy,NFFx, NFFy, NFGx, NFGy, NTin]);

%extracts force 
NForce_Ay= double(DynamicSolution.NFAy)
NForce_Ax= double(DynamicSolution.NFAx)
NForce_Cx= double(DynamicSolution.NFCx)
NForce_Cy= double(DynamicSolution.NFCy)
NForce_Bx= double(DynamicSolution.NFBx)
NForce_By= double(DynamicSolution.NFBy)
NForce_Dx= double(DynamicSolution.NFDx)
NForce_Dy= double(DynamicSolution.NFDy)
NForce_Ex= double(DynamicSolution.NFEx)
NForce_Ey= double(DynamicSolution.NFEy)
NForce_Fx= double(DynamicSolution.NFFx)
NForce_Fy= double(DynamicSolution.NFFy)
NForce_Gx= double(DynamicSolution.NFGx)
NForce_Gy= double(DynamicSolution.NFGy)

% displays velocities at joints 

disp(' ');
disp('______ VELOCITIES AT JOINTS _______');

disp('Velocity of Joint B:');
disp(vB_A);

disp('Velocity of Joint C:');
disp(vC_B);

disp('Velocity of Joint D:');
disp(vC_D);

disp('Velocity of Joint E:');
disp(VE_A);

disp('Velocity of Joint F:');
disp(vF_G);


% displays velocities at center of mass for Links AB, BEC,CD,EF,FG

disp(' ');
disp('________VELOCITIES AT CENTER OF MASS _________');

disp('Velocity of S1 (Link AB):');
disp(vS1_A);

disp('Velocity of S2 (Link BEC):');
disp(vS2_A);

disp('Velocity of S3 (Link CD):');
disp(vS3_D);

disp('Velocity of S4 (Link EF):');
disp(vS4_G);

disp('Velocity of S5 (Link FG):');
disp(vS5_G);

%displays acceleration at joints
disp(' ');
disp('_________ACCELERATIONS AT JOINTS _________');

disp('Acceleration of Joint B:');
disp(aB_A);

disp('Acceleration of Joint C:');
disp(aC_D);


disp('Acceleration of Joint E:');
disp(aE_D);

disp('Acceleration of Joint F:');
disp(aF_G);



%Displays acceleration at center of mass for Links AB, BEC,CD,EF,FG
disp(' ');
disp('===== ACCELERATIONS AT CENTER OF MASS _______');

disp('Acceleration of S1 (Link AB):');
disp(aS1_A);

disp('Acceleration of S2 (Link BEC):');
disp(aS2_A);

disp('Acceleration of S3 (Link CD):');
disp(aS3_D);

disp('Acceleration of S4 (Link EF):');
disp(aS4_G);


disp('Acceleration of S5 (Link FG):');
disp(aS5_G);