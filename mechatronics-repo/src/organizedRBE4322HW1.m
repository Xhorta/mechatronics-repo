%% ME/RBE 4322 - Homework 1
% Six-Bar Linkage Analysis
%used ChatGPT to reorganize/clean/debug/code assistance my code 
% Position, Velocity, Acceleration, Static Equilibrium, and Dynamic Force Analysis
% Full 0-360 degree input crank sweep in 1-degree increments

% Loop 1: A-B-C-D-A
% Loop 2: D-E-F-G-D through the rigid CDE ternary body

% Assumption: one full revolution of AB corresponds to one completed part
% Assumption: the gripper extends 1.843 m from F along the G->F direction

clear;
clc;
close all;

%% ASSUMPTIONS 
%% ============================================================

parts = 12500;                          %parts in 9 hours
operationTime = 9*3600;                % 9hours in seconds
cyclesPerSecond = parts/operationTime;
omegaInput = 2*pi*cyclesPerSecond;     % rad/s
alphaInput = 0;                       % constant input speed

g = 9.81;                              % m/s^2
lbToKg = 0.45359237;                   % conversion from lb to kg

% Link masses/s moments of inertia from on shape, converted from lb to kg
% Parts were modeled on the Front plane
% Lyy values were converted from lb*in^2 to kg*m^2

mAB  = 11.982  * lbToKg;               % kg
mBC  = 44.907  * lbToKg;               % kg
mCDE = 123.634 * lbToKg;               % kg
mEF  = 37.456  * lbToKg;               % kg
mFG  = 81.159  * lbToKg; % kg

IAB  = 0.1585100;  % kg*m^2
IBC  = 3.6210122;    % kg*m^2
ICDE = 32.5765628;    % kg*m^2
IEF  = 2.1145146;      % kg*m^2
IFG  = 21.1582145;   % kg*m^2

% Artifact mass corresponding to the specified 200 N weight
massArtifact = 200/g; % kg

fprintf('Required input angular velocity = %.4f rad/s\n',omegaInput);
fprintf('Link masses converted to SI units (kg).\n');

%% 
%% INITIAL JOINT COORDINATES
%% ============================================================
A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.20 0];

% 2-D for the sweep calculations
Axy = A(1:2);
Bxy = B(1:2);
Cxy = C(1:2);
Dxy = D(1:2);
Exy = E(1:2);
Fxy = F(1:2);
Gxy = G(1:2);

%% GRIPPER / ARTIFACT LOCATION
%% ============================================================
gripperDistance = 1.843;              % m from F
UV_GF = (F-G)/norm(F-G);             % unit vector from G toward F
H = F + gripperDistance*UV_GF;       % load location
% Artifact weight, downward
AppliedForce = [0 -200 0];             % N
%% LINK LENGTHS
%% ============================================================
lAB = norm(B-A);
lBC = norm(C-B);
lCD = norm(D-C);
lDE = norm(E-D);
lEF = norm(F-E);
lFG = norm(G-F);
lAD = norm(D-A);
lDG = norm(G-D);

fprintf('\n--- Link Lengths ---\n');
fprintf('AB = %.4f m\n',lAB);
fprintf('BC = %.4f m\n',lBC);
fprintf('CD = %.4f m\n',lCD);
fprintf('DE = %.4f m\n',lDE);
fprintf('EF = %.4f m\n',lEF);
fprintf('FG = %.4f m\n',lFG);

%% ============================================================
%% LINK WEIGHTS
%% ============================================================
% Weights use SI masses from Onshape.
WAB  = [0 -mAB*g  0];
WBC  = [0 -mBC*g  0];
WCDE = [0 -mCDE*g 0];
WEF  = [0 -mEF*g  0];
WFG  = [0 -mFG*g  0];
%% ============================================================
%% INITIAL CENTER-OF-MASS LOCATIONS
%% ============================================================
% Straight links are symmetric, so their COM is at the midpoint.
S1 = (A+B)/2;
S2 = (B+C)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

% CDE COM from Onshape.
% COM:  X = -2.273 in, Z = 40.183 in
% Joint D:X = -7.029 in, Z = -7.884 in
% Relative D->COM = [4.756, 48.067] in = [0.1208024, 1.2209018] m.
% Assumption: Onshape +X maps to MATLAB +x and Onshape +Z maps to MATLAB +y.
rD_S3_initial = [0.1208024, 1.2209018];  % m
S3xy = Dxy + rD_S3_initial;
S3 = [S3xy 0];


%% ============================================================
%% POSITION ANALYSIS - FULL 0 TO 360 DEGREE SWEEP
%% ============================================================

% Fixed angular offset between D->C and D->E on rigid ternary body C-D-E
theta_DC = atan2(Cxy(2)-Dxy(2),Cxy(1)-Dxy(1));
theta_DE = atan2(Exy(2)-Dxy(2),Exy(1)-Dxy(1));
offset = theta_DE-theta_DC;
theta_CDE_initial = theta_DC;

% Initial input angle A->B
theta2_AB = atan2(Bxy(2)-Axy(2),Bxy(1)-Axy(1));

% 0 through 360 degrees, 1 degree increment 
N = 361;
theta2 = theta2_AB + linspace(0,2*pi,N);
inputAngle = rad2deg(theta2-theta2_AB);

% Preallocate position arrays
A_pos = repmat(Axy,N,1); %grounded
B_pos = NaN(N,2);
C_pos = NaN(N,2);
D_pos = repmat(Dxy,N,1); %grounded
E_pos = NaN(N,2);
F_pos = NaN(N,2);
G_pos = repmat(Gxy,N,1); %grounded
H_pos = NaN(N,2);
% Preserve assembly 
C_prev = Cxy;
F_prev = Fxy;

for i = 1:N

    % Joint B is driven by the input crank AB
    Bi = Axy + lAB*[cos(theta2(i)),sin(theta2(i))];
    B_pos(i,:) = Bi;
    % Joint C: intersection of circle(B,lBC) and circle(D,lCD)

    Ci = circleIntersect(Bi,lBC,Dxy,lCD,C_prev);
    if isempty(Ci)
        warning('Loop ABCD cannot be solved at input angle %.1f deg',inputAngle(i));
        continue
    end
    C_pos(i,:) = Ci;
    C_prev = Ci;
    % Joint E rotates rigidly with C-D-E about grounded pivot D
    theta_DC_i = atan2(Ci(2)-Dxy(2),Ci(1)-Dxy(1));
    theta_DE_i = theta_DC_i + offset;

    Ei = Dxy + lDE*[cos(theta_DE_i),sin(theta_DE_i)];
    E_pos(i,:) = Ei;

    % Joint F: intersection of circle(E,lEF) and circle(G,lFG)
    Fi = circleIntersect(Ei,lEF,Gxy,lFG,F_prev);

    if isempty(Fi)
        warning('Loop DEFG cannot be solved at input angle %.1f deg',inputAngle(i));
        continue
    end

    F_pos(i,:) = Fi;
    F_prev = Fi;

    % Artifact location extends beyond F along G->F
    unit_GF_i = (Fi-Gxy)/norm(Fi-Gxy);
    H_pos(i,:) = Fi + gripperDistance*unit_GF_i;
end

%% ============================================================
%% VELOCITY ANALYSIS - ALL 361 POSITIONS
%% ============================================================

omega_AB = omegaInput;

% Angular velocities
omega_BC_all  = NaN(N,1);
omega_CDE_all = NaN(N,1);
omega_EF_all  = NaN(N,1);
omega_FG_all  = NaN(N,1);

% Joint velocities
VB_all = NaN(N,2);
VC_all = NaN(N,2);
VE_all = NaN(N,2);
VF_all = NaN(N,2);
VH_all = NaN(N,2);

% Mass-center velocities
VS1_all = NaN(N,2);
VS2_all = NaN(N,2);
VS3_all = NaN(N,2);
VS4_all = NaN(N,2);
VS5_all = NaN(N,2);

for i = 1:N

    Ai = A_pos(i,:);
    Bi = B_pos(i,:);
    Ci = C_pos(i,:);
    Di = D_pos(i,:);
    Ei = E_pos(i,:);
    Fi = F_pos(i,:);
    Gi = G_pos(i,:);
    Hi = H_pos(i,:);

    if any(isnan([Bi Ci Ei Fi Hi]))
        continue
    end

    % Loop 1: A-B-C-D-A
    rAB = Bi-Ai;
    rBC = Ci-Bi;
    rDC = Ci-Di;

    % Velocity of B
    VBi = crossZ(omega_AB,rAB);
    VB_all(i,:) = VBi;

    % V_B + omega_BC x rBC = omega_CDE x rDC
    M1 = [-rBC(2),  rDC(2); rBC(1), -rDC(1)];

    b1 = [-VBi(1); -VBi(2)];

    if rcond(M1) < 1e-10
        warning('Velocity Loop 1 is near singular at %.1f deg',inputAngle(i));
        continue
    end

    omegaSolution1 = M1\b1;
    omega_BC_all(i)  = omegaSolution1(1);
    omega_CDE_all(i) = omegaSolution1(2);
    % Velocity of C
    VCi = crossZ(omega_CDE_all(i),rDC);
    VC_all(i,:) = VCi;
    % Loop 2: D-E-F-G-D
    rDE = Ei-Di;
    rEF = Fi-Ei;
    rGF = Fi-Gi;
    % Velocity of E on rigid CDE body
    VEi = crossZ(omega_CDE_all(i),rDE);
    VE_all(i,:) = VEi;
    
    % V_E + omega_EF x rEF = omega_FG x rGF
    M2 = [-rEF(2),  rGF(2); rEF(1), -rGF(1)];

    b2 = [-VEi(1);-VEi(2)];

    if rcond(M2) < 1e-10
        warning('Velocity Loop 2 is near singular at %.1f deg',inputAngle(i));
        continue
    end

    omegaSolution2 = M2\b2;

    omega_EF_all(i) = omegaSolution2(1);
    omega_FG_all(i) = omegaSolution2(2);

    % Velocity of F
    VFi = crossZ(omega_FG_all(i),rGF);
    VF_all(i,:) = VFi;

    % Velocity of artifact point H extended pass F of link FG about G
    rGH = Hi-Gi;
    VHi = crossZ(omega_FG_all(i),rGH);
    VH_all(i,:) = VHi;

    %% Mass-center velocities
    S1i = (Ai+Bi)/2;
    S2i = (Bi+Ci)/2;
    theta_CDE_i = atan2(Ci(2)-Di(2),Ci(1)-Di(1));
    deltaTheta_CDE = theta_CDE_i-theta_CDE_initial;
    rD_S3_i = rotate2D(rD_S3_initial,deltaTheta_CDE);
    S3i = Di+rD_S3_i;
    S4i = (Ei+Fi)/2;
    S5i = (Fi+Gi)/2;

    VS1_all(i,:) = crossZ(omega_AB,S1i-Ai);
    VS2_all(i,:) = VBi + crossZ(omega_BC_all(i),S2i-Bi);
    VS3_all(i,:) = crossZ(omega_CDE_all(i),S3i-Di);
    VS4_all(i,:) = VEi + crossZ(omega_EF_all(i),S4i-Ei);
    VS5_all(i,:) = crossZ(omega_FG_all(i),S5i-Gi);

end

%% ============================================================
%% ACCELERATION ANALYSIS - ALL 361 POSITIONS
%% ============================================================

alpha_AB = alphaInput;

% Angular accelerations
alpha_BC_all  = NaN(N,1);
alpha_CDE_all = NaN(N,1);
alpha_EF_all  = NaN(N,1);
alpha_FG_all  = NaN(N,1);

% Joint accelerations
aB_all = NaN(N,2);
aC_all = NaN(N,2);
aE_all = NaN(N,2);
aF_all = NaN(N,2);
aH_all = NaN(N,2);

% Mass-center accelerations
aS1_all= NaN(N,2);
aS2_all= NaN(N,2);
aS3_all= NaN(N,2);
aS4_all= NaN(N,2);
aS5_all= NaN(N,2);

for i = 1:N
    Ai = A_pos(i,:);
    Bi = B_pos(i,:);
    Ci = C_pos(i,:);
    Di = D_pos(i,:);
    Ei = E_pos(i,:);
    Fi = F_pos(i,:);
    Gi = G_pos(i,:);
    Hi = H_pos(i,:);

    if any(isnan([Bi Ci Ei Fi Hi]))
        continue
    end

    wBC  = omega_BC_all(i);
    wCDE = omega_CDE_all(i);
    wEF  = omega_EF_all(i);
    wFG  = omega_FG_all(i);

    if any(isnan([wBC wCDE wEF wFG]))
        continue
    end

    %% Loop 1
    rAB = Bi-Ai;
    rBC = Ci-Bi;
    rDC = Ci-Di;
    % a_B = alpha_AB x rAB - omega_AB^2*rAB
    aBi = crossZ(alpha_AB,rAB) - omega_AB^2*rAB;
    aB_all(i,:) = aBi;

    M1 = [-rBC(2),  rDC(2); rBC(1), -rDC(1)];
    % aB + alphaBC x rBC - wBC^2*rBC = alphaCDE x rDC - wCDE^2*rDC
    RHS1 = -aBi + wBC^2*rBC - wCDE^2*rDC;

    if rcond(M1) < 1e-10
        warning('Acceleration Loop 1 is near singular at %.1f deg',inputAngle(i));
        continue
    end

    alphaSolution1 = M1\RHS1(:);

    alpha_BC_all(i) = alphaSolution1(1);
    alpha_CDE_all(i)= alphaSolution1(2);

    % Acceleration of C
    aCi = crossZ(alpha_CDE_all(i),rDC) - wCDE^2*rDC;
    aC_all(i,:) = aCi;

    %% Loop 2
    rDE = Ei-Di;
    rEF = Fi-Ei;
    rGF = Fi-Gi;

    % Acceleration of E
    aEi = crossZ(alpha_CDE_all(i),rDE) - wCDE^2*rDE;
    aE_all(i,:) = aEi;

    M2 = [-rEF(2),  rGF(2);
           rEF(1), -rGF(1)];

    RHS2 = -aEi + wEF^2*rEF - wFG^2*rGF;

    if rcond(M2) < 1e-10
        warning('Acceleration Loop 2 is near singular at %.1f deg',inputAngle(i));
        continue
    end

    alphaSolution2 = M2\RHS2(:);

    alpha_EF_all(i) = alphaSolution2(1);
    alpha_FG_all(i) = alphaSolution2(2);

    % Acceleration of F
    aFi = crossZ(alpha_FG_all(i),rGF) - wFG^2*rGF;
    aF_all(i,:) = aFi;

    % Acceleration of artifact point H
    rGH = Hi-Gi;
    aHi = crossZ(alpha_FG_all(i),rGH) - wFG^2*rGH;
    aH_all(i,:) = aHi;

    % Mass-center accelerations
    S1i = (Ai+Bi)/2;
    S2i = (Bi+Ci)/2;
    theta_CDE_i = atan2(Ci(2)-Di(2),Ci(1)-Di(1));
    deltaTheta_CDE = theta_CDE_i-theta_CDE_initial;
    rD_S3_i = rotate2D(rD_S3_initial,deltaTheta_CDE);
    S3i = Di+rD_S3_i;
    S4i = (Ei+Fi)/2;
    S5i = (Fi+Gi)/2;

    % Compute mass-center accelerations: tangential (alpha x r) and centripetal (-omega^2 r)
    aS1_all(i,:)= crossZ(alpha_AB,S1i-Ai) - omega_AB^2*(S1i-Ai);

    aS2_all(i,:)= aBi + crossZ(alpha_BC_all(i),S2i-Bi) - wBC^2*(S2i-Bi);

    aS3_all(i,:)= crossZ(alpha_CDE_all(i),S3i-Di) - wCDE^2*(S3i-Di);

    aS4_all(i,:)= aEi+ crossZ(alpha_EF_all(i),S4i-Ei)- wEF^2*(S4i-Ei);

    aS5_all(i,:)= crossZ(alpha_FG_all(i),S5i-Gi)- wFG^2*(S5i-Gi);
end

%% ============================================================
%% STATIC EQUILIBRIUM - ALL 361 POSITIONS
%% ============================================================
% Unknown vector for each position:
% x = [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]^T

FA_static = NaN(N,2);
FB_static = NaN(N,2);
FC_static = NaN(N,2);
FD_static = NaN(N,2);
FE_static = NaN(N,2);
FF_static = NaN(N,2);
FG_static = NaN(N,2);
Tin_static = NaN(N,1);

% For r x F = rx*Fy - ry*Fx
momentCoeff = @(r)[-r(2),r(1)];

for i = 1:N

    Ai = A_pos(i,:);
    Bi = B_pos(i,:);
    Ci = C_pos(i,:);
    Di = D_pos(i,:);
    Ei = E_pos(i,:);
    Fi = F_pos(i,:);
    Gi = G_pos(i,:);
    Hi = H_pos(i,:);

    if any(isnan([Bi Ci Ei Fi Hi]))
        continue
    end

    % Current mass-center locations
    S1i=(Ai+Bi)/2;
    %center of mass BC
    S2i=(Bi+Ci)/2;
    theta_CDE_i = atan2(Ci(2)-Di(2),Ci(1)-Di(1));
    deltaTheta_CDE = theta_CDE_i-theta_CDE_initial;
    rD_S3_i = rotate2D(rD_S3_initial,deltaTheta_CDE);
    %center of mass DE
    S3i=Di+rD_S3_i;
    %center of Mass EF
    S4i=(Ei+Fi)/2;
    %center of Mass FG
    S5i=(Fi+Gi)/2;

    % 15 equations, 15 unknowns
    A_static = zeros(15,15);
    B_static = zeros(15,1);

    %% Link AB: FA + FB + WAB = 0
    A_static(1,1)= 1;        % FAx
    A_static(1,3)= 1;        % FBx
    B_static(1) =-WAB(1);

    A_static(2,2)= 1;        % FAy
    A_static(2,4)= 1;        % FBy
    B_static(2)= -WAB(2);

    % Moment about S1
    A_static(3,1:2)= momentCoeff(Ai-S1i); 
    A_static(3,3:4)=momentCoeff(Bi-S1i);
    A_static(3,15)=1;       % Tin
    B_static(3) = 0;

    %% Link BC: -FB + FC + WBC = 0
    A_static(4,3)= -1;  %-FBx
    A_static(4,5)= 1;   %FCx
    B_static(4)= -WBC(1);

    A_static(5,4)=-1;  %-FBy
    A_static(5,6)= 1;  %FCy
    B_static(5)= -WBC(2);

    % Moment about S2
    A_static(6,3:4)= -momentCoeff(Bi-S2i);
    A_static(6,5:6)= momentCoeff(Ci-S2i);
    B_static(6)= 0;

    %% Link CDE: -FC + FD - FE + WCDE = 0
    A_static(7,5)=-1; %-FCx
    A_static(7,7)=1;    %FDx
    A_static(7,9)=-1; %-FEx
    B_static(7)=-WCDE(1);

    A_static(8,6) = -1; % -FCy
    A_static(8,8) = 1;  %FDy
    A_static(8,10)=-1;  %-FEy
    B_static(8)=-WCDE(2);

    % Moment about S3
    A_static(9,5:6)=-momentCoeff(Ci-S3i);
    A_static(9,7:8)= momentCoeff(Di-S3i);
    A_static(9,9:10)=-momentCoeff(Ei-S3i);
    B_static(9)= 0;

    %% Link EF: FE - FF + WEF = 0
    A_static(10,9)=  1;  %FEx
    A_static(10,11)= -1; %-FFx
    B_static(10)= -WEF(1);

    A_static(11,10)=1;  %FEx
    A_static(11,12)=-1;  %-FFy
    B_static(11) = -WEF(2);

    % Moment about S4
    A_static(12,9:10)= momentCoeff(Ei-S4i);
    A_static(12,11:12)= -momentCoeff(Fi-S4i);
    B_static(12)= 0;

    %% Link FG: FF - FG + WFG + AppliedForce = 0
    A_static(13,11)= 1;  %FFx
    A_static(13,13)=-1;   %-FGx
    B_static(13)= -(WFG(1)+AppliedForce(1));

    A_static(14,12)=  1; %FFy
    A_static(14,14) = -1; %-FGy
    B_static(14) = -(WFG(2)+AppliedForce(2));

    % Moment about S5
    A_static(15,11:12) =  momentCoeff(Fi-S5i);
    A_static(15,13:14) = -momentCoeff(Gi-S5i);

    artifactMoment = (Hi(1)-S5i(1))*AppliedForce(2) - (Hi(2)-S5i(2))*AppliedForce(1);

    B_static(15) = -artifactMoment;

    %% Solve this position
    if rcond(A_static) < 1e-12
        warning('Static eq matrix is near singular at %.1f deg',inputAngle(i));
        continue
    end

    staticSolution = A_static\B_static;

    % Store results
    FA_static(i,:) =staticSolution(1:2).';
    FB_static(i,:)=staticSolution(3:4).';
    FC_static(i,:)=staticSolution(5:6).';
    FD_static(i,:)=staticSolution(7:8).';
    FE_static(i,:)=staticSolution(9:10).';
    FF_static(i,:)=staticSolution(11:12).';
    FG_static(i,:)=staticSolution(13:14).';
    Tin_static(i) =staticSolution(15).';
end

%% ============================================================
%% DYNAMIC FORCE / TORQUE ANALYSIS - NEWTON'S SECOND LAW
%% ============================================================
% Unknown vector for each position:
% [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]^T
% For each moving link:
% sum(F) = m*a_S
% sum(M about mass center) = I*alpha
% The 200 N artifact is treated as a point mass rigidly carried at H.
% force on the gripper is W_artifact - mass_artifact*a_H.

FA_dynamic = NaN(N,2);
FB_dynamic = NaN(N,2);
FC_dynamic = NaN(N,2);
FD_dynamic = NaN(N,2);
FE_dynamic = NaN(N,2);
FF_dynamic = NaN(N,2);
FG_dynamic = NaN(N,2);
Tin_dynamic = NaN(N,1);
dynamic_ArtifactForce = NaN(N,2);

for i = 1:N
    Ai = A_pos(i,:);
    Bi = B_pos(i,:);
    Ci = C_pos(i,:);
    Di = D_pos(i,:);
    Ei = E_pos(i,:);
    Fi = F_pos(i,:);
    Gi = G_pos(i,:);
    Hi = H_pos(i,:);

    if any(isnan([Bi Ci Ei Fi Hi])) ||any(isnan([aS1_all(i,:) aS2_all(i,:) aS3_all(i,:) aS4_all(i,:) aS5_all(i,:) aH_all(i,:)])) || ...
       any(isnan([alpha_BC_all(i) alpha_CDE_all(i) alpha_EF_all(i) alpha_FG_all(i)]))
        continue
    end

    % Current centers of mass
    S1i = (Ai+Bi)/2;
    S2i = (Bi+Ci)/2;
    theta_CDE_i = atan2(Ci(2)-Di(2),Ci(1)-Di(1));
    deltaTheta_CDE = theta_CDE_i-theta_CDE_initial;
    rD_S3i = rotate2D(rD_S3_initial,deltaTheta_CDE);
    S3i = Di+rD_S3i;
    S4i = (Ei+Fi)/2;
    S5i = (Fi+Gi)/2;

    % Dynamic force exerted by the artifact on the gripper/link FG.
    % Artifact equation: F_gripper_artifact + W = m*a_H
    % Therefore artifact_on_gripper = W - m*a_H.
    F_artifact = AppliedForce(1:2) - massArtifact*aH_all(i,:);
    dynamic_ArtifactForce(i,:) = F_artifact;

    A_dynamic = zeros(15,15);
    B_dynamic = zeros(15,1);

    %% Link AB
    A_dynamic(1,1)=1;  
    A_dynamic(1,3)=1;
    B_dynamic(1)=mAB*aS1_all(i,1)-WAB(1); 

    A_dynamic(2,2)=1; 
    A_dynamic(2,4)=1;
    B_dynamic(2)=mAB*aS1_all(i,2)-WAB(2);

    A_dynamic(3,1:2)=momentCoeff(Ai-S1i);
    A_dynamic(3,3:4)=momentCoeff(Bi-S1i);
    A_dynamic(3,15)=1;
    B_dynamic(3)=IAB*alpha_AB;

    %% Link BC
    A_dynamic(4,3)=-1; 
    A_dynamic(4,5)=1;
    B_dynamic(4)=mBC*aS2_all(i,1)-WBC(1);

    A_dynamic(5,4)=-1; A_dynamic(5,6)=1;
    B_dynamic(5)=mBC*aS2_all(i,2)-WBC(2);

    A_dynamic(6,3:4)=-momentCoeff(Bi-S2i);
    A_dynamic(6,5:6)= momentCoeff(Ci-S2i);
    B_dynamic(6)=IBC*alpha_BC_all(i);

    %% Link CDE
    A_dynamic(7,5)=-1; 
    A_dynamic(7,7)=1;
    A_dynamic(7,9)=-1;
    B_dynamic(7)=mCDE*aS3_all(i,1)-WCDE(1);

    A_dynamic(8,6)=-1; 
    A_dynamic(8,8)=1;
    A_dynamic(8,10)=-1;
    B_dynamic(8)=mCDE*aS3_all(i,2)-WCDE(2);

    A_dynamic(9,5:6)=-momentCoeff(Ci-S3i);
    A_dynamic(9,7:8)= momentCoeff(Di-S3i);
    A_dynamic(9,9:10)=-momentCoeff(Ei-S3i);
    B_dynamic(9)=ICDE*alpha_CDE_all(i);

    %% Link EF
    A_dynamic(10,9)=1; 
    A_dynamic(10,11)=-1;

    B_dynamic(10)=mEF*aS4_all(i,1)-WEF(1);

    A_dynamic(11,10)=1; A_dynamic(11,12)=-1;
    B_dynamic(11)=mEF*aS4_all(i,2)-WEF(2);

    A_dynamic(12,9:10)= momentCoeff(Ei-S4i);
    A_dynamic(12,11:12)=-momentCoeff(Fi-S4i);
    B_dynamic(12)=IEF*alpha_EF_all(i);

    %% Link FG
    % FF - FG + WFG + F_artifact = mFG*aS5
    A_dynamic(13,11)=1; 
    A_dynamic(13,13)=-1;
    B_dynamic(13)=mFG*aS5_all(i,1)-WFG(1)-F_artifact(1);

    A_dynamic(14,12)=1; 
    A_dynamic(14,14)=-1;
    B_dynamic(14)=mFG*aS5_all(i,2)-WFG(2)-F_artifact(2);

    % Moments about S5
    A_dynamic(15,11:12)= momentCoeff(Fi-S5i);
    A_dynamic(15,13:14)=-momentCoeff(Gi-S5i);
    artifactMomentDynamic=(Hi(1)-S5i(1))*F_artifact(2) - (Hi(2)-S5i(2))*F_artifact(1);
    B_dynamic(15)=IFG*alpha_FG_all(i)-artifactMomentDynamic;

    if rcond(A_dynamic) < 1e-12
        warning('Dynamic matrix is near singular at %.1f deg',inputAngle(i));
        continue
    end

    dynamicSolution=A_dynamic\B_dynamic;
    
    FA_dynamic(i,:)=dynamicSolution(1:2).';
    FB_dynamic(i,:)=dynamicSolution(3:4).';
    FC_dynamic(i,:)=dynamicSolution(5:6).';
    FD_dynamic(i,:)=dynamicSolution(7:8).';
    FE_dynamic(i,:)=dynamicSolution(9:10).';
    FF_dynamic(i,:)=dynamicSolution(11:12).';
    FG_dynamic(i,:)=dynamicSolution(13:14).';
    Tin_dynamic(i)=dynamicSolution(15);
end

%% FIRST-POSITION RESULTS
%% ============================================================
fprintf('FIRST POSITION RESULTS\n');
fprintf('============================================================\n');

fprintf('\n--- Joint Coordinates ---\n');
fprintf('A = (%.4f, %.4f) m\n',A_pos(1,:));
fprintf('B = (%.4f, %.4f) m\n',B_pos(1,:));
fprintf('C = (%.4f, %.4f) m\n',C_pos(1,:));
fprintf('D = (%.4f, %.4f) m\n',D_pos(1,:));
fprintf('E = (%.4f, %.4f) m\n',E_pos(1,:));
fprintf('F = (%.4f, %.4f) m\n',F_pos(1,:));
fprintf('G = (%.4f, %.4f) m\n',G_pos(1,:));
fprintf('H = (%.4f, %.4f) m  [artifact location]\n',H_pos(1,:));

fprintf('\n--- Angular Velocities ---\n');
fprintf('omega AB  = %.4f rad/s\n',omega_AB);
fprintf('omega BC  = %.4f rad/s\n',omega_BC_all(1));
fprintf('omega CDE = %.4f rad/s\n',omega_CDE_all(1));
fprintf('omega EF  = %.4f rad/s\n',omega_EF_all(1));
fprintf('omega FG  = %.4f rad/s\n',omega_FG_all(1));

fprintf('\n--- Joint Velocities ---\n');
fprintf('VB = (%.4f, %.4f) m/s\n',VB_all(1,:));
fprintf('VC = (%.4f, %.4f) m/s\n',VC_all(1,:));
fprintf('VE = (%.4f, %.4f) m/s\n',VE_all(1,:));
fprintf('VF = (%.4f, %.4f) m/s\n',VF_all(1,:));

fprintf('\n--- Angular Accelerations ---\n');
fprintf('alpha BC  = %.4f rad/s^2\n',alpha_BC_all(1));
fprintf('alpha CDE = %.4f rad/s^2\n',alpha_CDE_all(1));
fprintf('alpha EF  = %.4f rad/s^2\n',alpha_EF_all(1));
fprintf('alpha FG  = %.4f rad/s^2\n',alpha_FG_all(1));

fprintf('\n--- Joint Accelerations ---\n');
fprintf('aB = (%.4f, %.4f) m/s^2\n',aB_all(1,:));
fprintf('aC = (%.4f, %.4f) m/s^2\n',aC_all(1,:));
fprintf('aE = (%.4f, %.4f) m/s^2\n',aE_all(1,:));
fprintf('aF = (%.4f, %.4f) m/s^2\n',aF_all(1,:));

fprintf('\n--- Mass-Center Accelerations ---\n');
fprintf('aS1 = (%.4f, %.4f) m/s^2\n',aS1_all(1,:));
fprintf('aS2 = (%.4f, %.4f) m/s^2\n',aS2_all(1,:));
fprintf('aS3 = (%.4f, %.4f) m/s^2\n',aS3_all(1,:));
fprintf('aS4 = (%.4f, %.4f) m/s^2\n',aS4_all(1,:));
fprintf('aS5 = (%.4f, %.4f) m/s^2\n',aS5_all(1,:));

fprintf('\n--- Static Joint Forces ---\n');
fprintf('FA = (%.4f, %.4f) N\n',FA_static(1,:));
fprintf('FB = (%.4f, %.4f) N\n',FB_static(1,:));
fprintf('FC = (%.4f, %.4f) N\n',FC_static(1,:));
fprintf('FD = (%.4f, %.4f) N\n',FD_static(1,:));
fprintf('FE = (%.4f, %.4f) N\n',FE_static(1,:));
fprintf('FF = (%.4f, %.4f) N\n',FF_static(1,:));
fprintf('FG = (%.4f, %.4f) N\n',FG_static(1,:));
fprintf('Input Torque = %.4f N*m\n',Tin_static(1));

fprintf('\n--- Artifact Kinematics ---\n');
fprintf('VH = (%.4f, %.4f) m/s\n',VH_all(1,:));
fprintf('aH = (%.4f, %.4f) m/s^2\n',aH_all(1,:));

fprintf('\n--- Dynamic Joint Forces / Torque ---\n');
fprintf('FA = (%.4f, %.4f) N\n',FA_dynamic(1,:));
fprintf('FB = (%.4f, %.4f) N\n',FB_dynamic(1,:));
fprintf('FC = (%.4f, %.4f) N\n',FC_dynamic(1,:));
fprintf('FD = (%.4f, %.4f) N\n',FD_dynamic(1,:));
fprintf('FE = (%.4f, %.4f) N\n',FE_dynamic(1,:));
fprintf('FF = (%.4f, %.4f) N\n',FF_dynamic(1,:));
fprintf('FG = (%.4f, %.4f) N\n',FG_dynamic(1,:));
fprintf('Dynamic Input Torque = %.4f N*m\n',Tin_dynamic(1));

%% PLOTS
%% ============================================================

% ---------- Kinematic outline: first position ----------
figure('Name','Kinematic Outline - Inital Position');
hold on;
axis equal;
grid on;

plot([A_pos(1,1) B_pos(1,1)],[A_pos(1,2) B_pos(1,2)],'-o','LineWidth',2,'DisplayName','AB');
plot([B_pos(1,1) C_pos(1,1)],[B_pos(1,2) C_pos(1,2)],'-o','LineWidth',2,'DisplayName','BC');
plot([C_pos(1,1) D_pos(1,1) E_pos(1,1)],[C_pos(1,2) D_pos(1,2) E_pos(1,2)],'-o','LineWidth',2,'DisplayName','CDE');
plot([E_pos(1,1) F_pos(1,1)],[E_pos(1,2) F_pos(1,2)],'-o','LineWidth',2,'DisplayName','EF');
plot([G_pos(1,1) F_pos(1,1)],[G_pos(1,2) F_pos(1,2)],'-o','LineWidth',2,'DisplayName','GF');
plot([F_pos(1,1) H_pos(1,1)],[F_pos(1,2) H_pos(1,2)],'--','LineWidth',1.5,'DisplayName','Gripper extension');

plot(A_pos(1,1),A_pos(1,2),'ks','MarkerFaceColor','k','DisplayName','Ground A');
plot(D_pos(1,1),D_pos(1,2),'ks','MarkerFaceColor','k','DisplayName','Ground D');
plot(G_pos(1,1),G_pos(1,2),'ks','MarkerFaceColor','k','DisplayName','Ground G');

joint = {'A','B','C','D','E','F','G','H'};
pointsFirst = [A_pos(1,:);B_pos(1,:);C_pos(1,:);D_pos(1,:);E_pos(1,:);F_pos(1,:);G_pos(1,:);H_pos(1,:)];

for k = 1:size(pointsFirst,1)
    text(pointsFirst(k,1)+0.03,pointsFirst(k,2)+0.03,joint{k});
end


xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Six-Bar Linkage - inital Position');
legend('Location','best');

% ---------- Joint position traces ----------
figure('Name','Joint Position Traces');
hold on;
axis equal;
grid on;
plot(B_pos(:,1),B_pos(:,2),'LineWidth',1.5,'DisplayName','B');
plot(C_pos(:,1),C_pos(:,2),'LineWidth',1.5,'DisplayName','C');
plot(E_pos(:,1),E_pos(:,2),'LineWidth',1.5,'DisplayName','E');
plot(F_pos(:,1),F_pos(:,2),'LineWidth',1.5,'DisplayName','F');
plot(H_pos(:,1),H_pos(:,2),'LineWidth',1.5,'DisplayName','Artifact H');
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Joint / Artifact Position ');
legend('Location','best');

% ---------- Angular velocities ----------
figure('Name','Angular Velocities');
plot(inputAngle,omega_BC_all,'LineWidth',1.5);
hold on;
plot(inputAngle,omega_CDE_all,'LineWidth',1.5);
plot(inputAngle,omega_EF_all,'LineWidth',1.5);
plot(inputAngle,omega_FG_all,'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocities vs Input Crank Rotation');
legend('BC','CDE','EF','FG','Location','best');

% ---------- Joint speed magnitudes ----------
figure('Name','Joint Velocity Magnitudes');
plot(inputAngle,vecnorm(VB_all,2,2),'LineWidth',1.5);
hold on;
plot(inputAngle,vecnorm(VC_all,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(VE_all,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(VF_all,2,2),'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Joint Speed (m/s)');
title('Joint Velocity Magnitudes vs Input Crank Rotation');
legend('V_B','V_C','V_E','V_F','Location','best');

% ---------- Angular accelerations ----------
figure('Name','Angular Accelerations');
plot(inputAngle,alpha_BC_all,'LineWidth',1.5);
hold on;
plot(inputAngle,alpha_CDE_all,'LineWidth',1.5);
plot(inputAngle,alpha_EF_all,'LineWidth',1.5);
plot(inputAngle,alpha_FG_all,'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Angular Acceleration (rad/s^2)');
title('Angular Accelerations vs Input Crank Rotation');
legend('BC','CDE','EF','FG','Location','best');

% ---------- Joint acceleration magnitudes ----------
figure('Name','Joint Acceleration Magnitudes');
plot(inputAngle,vecnorm(aB_all,2,2),'LineWidth',1.5);
hold on;
plot(inputAngle,vecnorm(aC_all,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(aE_all,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(aF_all,2,2),'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Joint Acceleration (m/s^2)');
title('Joint Acceleration Magnitudes vs Input Crank Rotation');
legend('a_B','a_C','a_E','a_F','Location','best');

% ---------- Static input torque ----------
figure('Name','Static Input Torque');
plot(inputAngle,Tin_static,'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Static Input Torque (N*m)');
title('Static Input Torque vs Input Crank Rotation');

% ---------- Static joint force magnitudes ----------
figure('Name','Static Joint Force Magnitudes');
plot(inputAngle,vecnorm(FA_static,2,2),'LineWidth',1.5);
hold on;
plot(inputAngle,vecnorm(FB_static,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FC_static,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FD_static,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FE_static,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FF_static,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FG_static,2,2),'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Static Joint Force Magnitude (N)');
title('Static Joint Force Magnitudes vs Input Crank Rotation');
legend('A','B','C','D','E','F','G','Location','best');


% ---------- Dynamic input torque ----------
figure('Name','Dynamic Input Torque');
plot(inputAngle,Tin_dynamic,'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Dynamic Input Torque (N*m)');
title('Dynamic Input Torque vs Input Crank Rotation');

% ---------- Dynamic joint force magnitudes ----------
figure('Name','Dynamic Joint Force Magnitudes');
plot(inputAngle,vecnorm(FA_dynamic,2,2),'LineWidth',1.5);
hold on;
plot(inputAngle,vecnorm(FB_dynamic,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FC_dynamic,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FD_dynamic,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FE_dynamic,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FF_dynamic,2,2),'LineWidth',1.5);
plot(inputAngle,vecnorm(FG_dynamic,2,2),'LineWidth',1.5);
grid on;
xlabel('Input Crank Rotation (deg)');
ylabel('Dynamic Joint Force Magnitude (N)');
title('Dynamic Joint Force Magnitudes vs Input Crank Rotation');
legend('A','B','C','D','E','F','G','Location','best');

%% STATIC MAXIMUM / MINIMUM TORQUE SUMMARY
%% ============================================================

Torque = ~isnan(Tin_static);

if any(Torque)
    validIndices = find(Torque);
    [maxTorqueValue,localMaxIndex] = max(Tin_static(Torque));
    [minTorqueValue,localMinIndex] = min(Tin_static(Torque));
    [maxAbsTorqueValue,localAbsIndex] = max(abs(Tin_static(Torque)));

    maxTorqueIndex = validIndices(localMaxIndex);
    minTorqueIndex = validIndices(localMinIndex);
    maxAbsTorqueIndex = validIndices(localAbsIndex);

    fprintf('\n--- Static Input Torque Extremes ---\n');
    fprintf('Maximum signed torque = %.4f N*m at %.1f deg\n', maxTorqueValue,inputAngle(maxTorqueIndex));
    fprintf('Minimum signed torque = %.4f N*m at %.1f deg\n', minTorqueValue,inputAngle(minTorqueIndex));
    fprintf('Maximum absolute torque = %.4f N*m at %.1f deg\n', maxAbsTorqueValue,inputAngle(maxAbsTorqueIndex));
end

%% DYNAMIC MAXIMUM / MINIMUM TORQUE SUMMARY
DynamicTorqueValid = ~isnan(Tin_dynamic);
if any(DynamicTorqueValid)
    validIndices=find(DynamicTorqueValid);
    [maxDynamicTorqueValue,localMaxIndex]=max(Tin_dynamic(DynamicTorqueValid));
    [minDynamicTorqueValue,localMinIndex]=min(Tin_dynamic(DynamicTorqueValid));
    [maxAbsDynamicTorqueValue,localAbsIndex]=max(abs(Tin_dynamic(DynamicTorqueValid)));

    maxDynamicTorqueIndex=validIndices(localMaxIndex);
    minDynamicTorqueIndex=validIndices(localMinIndex);
    maxAbsDynamicTorqueIndex=validIndices(localAbsIndex);

    fprintf('\n--- Dynamic Input Torque Extremes ---\n');
    fprintf('Maximum signed torque = %.4f N*m at %.1f deg\n',maxDynamicTorqueValue,inputAngle(maxDynamicTorqueIndex));
    fprintf('Minimum signed torque = %.4f N*m at %.1f deg\n',minDynamicTorqueValue,inputAngle(minDynamicTorqueIndex));
    fprintf('Maximum absolute torque = %.4f N*m at %.1f deg\n',maxAbsDynamicTorqueValue,inputAngle(maxAbsDynamicTorqueIndex));
end

% Largest static and dynamic joint-force magnitudes for discussion
StaticForceMag=[vecnorm(FA_static,2,2),vecnorm(FB_static,2,2),vecnorm(FC_static,2,2),vecnorm(FD_static,2,2),vecnorm(FE_static,2,2),vecnorm(FF_static,2,2),vecnorm(FG_static,2,2)];

DynamicForceMag=[vecnorm(FA_dynamic,2,2),vecnorm(FB_dynamic,2,2),vecnorm(FC_dynamic,2,2),vecnorm(FD_dynamic,2,2),vecnorm(FE_dynamic,2,2),vecnorm(FF_dynamic,2,2),vecnorm(FG_dynamic,2,2)];
jointNames={'A','B','C','D','E','F','G'};

fprintf('\n--- Maximum Joint Force Magnitudes ---\n');
for j=1:7
    fprintf('Joint %s: static max = %.3f N, dynamic max = %.3f N\n',...
        jointNames{j},max(StaticForceMag(:,j),[],'omitnan'),max(DynamicForceMag(:,j),[],'omitnan'));
end

%% SAVE RESULTS
%% ============================================================

save('HW1_AllResults.mat', 'theta2','theta2_AB','inputAngle','omegaInput','alphaInput','A_pos','B_pos','C_pos','D_pos','E_pos','F_pos','G_pos','H_pos', ...
    'omega_BC_all','omega_CDE_all','omega_EF_all','omega_FG_all', 'VB_all','VC_all','VE_all','VF_all','VH_all', 'VS1_all','VS2_all','VS3_all','VS4_all','VS5_all', ...
    'alpha_BC_all','alpha_CDE_all','alpha_EF_all','alpha_FG_all','aB_all','aC_all','aE_all','aF_all','aH_all','aS1_all','aS2_all','aS3_all','aS4_all','aS5_all', ...
    'FA_static','FB_static','FC_static','FD_static','FE_static','FF_static','FG_static','Tin_static', ...
    'FA_dynamic','FB_dynamic','FC_dynamic','FD_dynamic','FE_dynamic','FF_dynamic','FG_dynamic','Tin_dynamic','dynamic_ArtifactForce', ...
    'mAB','mBC','mCDE','mEF','mFG','massArtifact','IAB','IBC','ICDE','IEF','IFG','rD_S3_initial');

fprintf('\nAnalysis complete. Results saved to HW1_AllResults.mat\n');

%% LOCAL FUNCTIONS
%% ============================================================

function P = circleIntersect(C1,r1,C2,r2,prevGuess)
    % uses the circle intersection method of two circles and chooses the point closest
    d = norm(C2-C1);
    if d > (r1+r2) || d < abs(r1-r2) || d == 0
        P = [];
        return
    end

    a = (r1^2-r2^2+d^2)/(2*d);
    h2 = r1^2-a^2;

    %protects against tiny negative roundoff values.
    if h2 <-1e-10
        P = [];
        return
    end
    h2 = max(h2,0);
    h = sqrt(h2);

    mid = C1 + a*(C2-C1)/d;
    perp = [-(C2(2)-C1(2)),(C2(1)-C1(1))]/d;

    P1 = mid + h*perp;
    P2 = mid - h*perp;

    if norm(P1-prevGuess) <= norm(P2-prevGuess)
        P = P1;
    else
        P = P2;
    end
end

function v = crossZ(omega,r)
    v = omega*[-r(2),r(1)];
end

function rRot = rotate2D(r,theta)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rRot = (R*r(:)).';
end