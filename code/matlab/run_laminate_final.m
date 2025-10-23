
%% 1. Material properties
Ef = 85.5;  Em = 4.0;      % GPa
Gf = 26.72; Gm = 1.44;     % GPa
vf = 0.22;  vm = 0.39;
Vf = 0.6;   Vm = 1 - Vf;

E1  = Ef*Vf + Em*Vm;
E2  = (Ef*Em)/(Ef*Vm + Em*Vf);
G12 = (Gf*Gm)/(Gf*Vm + Gm*Vf);
nu12 = vf*Vf + vm*Vm;
nu21 = nu12 * E2 / E1;

Q11 = E1/(1 - nu12*nu21);
Q22 = E2/(1 - nu12*nu21);
Q12 = nu12*E2/(1 - nu12*nu21);
Q66 = G12;
Qred = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];  % kN/mm²

%% 2. Stack definitions
plyThk = 0.4;
stacks(1).name = '[+60/0/−60] Anti-symmetric';
stacks(1).angles = [ 60  0 -60];
stacks(2).name = '[0/+60/−60] Asymmetric';
stacks(2).angles = [  0 60 -60];

%% 3. Loop through both stacks
for s = 1:2
    ang = stacks(s).angles;
    n   = numel(ang);
    z   = linspace(-n*plyThk/2, n*plyThk/2, n+1);

    [A,B,D] = buildABD(Qred, ang, z);
    ABD = [A B; B D];

    stacks(s).A = A; stacks(s).B = B; stacks(s).D = D;

    % Apply in-plane & bending loads
    N = [1; 0; 0];     % kN/mm
    M = [10; 0; 0];    % kN·mm/mm
    loadVec = [N; M];
    x = ABD \ loadVec;
    eps0 = x(1:3); kappa = x(4:6);
    
    % Stress and strain through thickness
    zEval = linspace(z(1), z(end), 200);
    strainZ = zeros(3,numel(zEval));
    stressZ = zeros(3,numel(zEval));

    for j = 1:numel(zEval)
        zj = zEval(j);
        strain = eps0 + zj * kappa;
        ply = find(z(1:end-1) <= zj & z(2:end) >= zj, 1);
        Qbar = QbarCS(Qred, ang(ply));
        stress = Qbar * strain;
        strainZ(:,j) = strain;
        stressZ(:,j) = stress;
    end

    % Store for plotting
    stacks(s).zEval = zEval;
    stacks(s).strainZ = strainZ;
    stacks(s).stressZ = stressZ;
end

%% 4. Plotting
for s = 1:2
    zEval = stacks(s).zEval;
    strainZ = stacks(s).strainZ;
    stressZ = stacks(s).stressZ;

    figure('Name', ['Strain - ' stacks(s).name]);
    plot(strainZ(1,:),zEval,'r', strainZ(2,:),zEval,'g', strainZ(3,:),zEval,'b');
    xlabel('Strain'); ylabel('z (mm)'); grid on;
    legend('\epsilon_x','\epsilon_y','\gamma_{xy}');
    title(['Mid-plane Strain Distribution: ' stacks(s).name]);

    figure('Name', ['Stress - ' stacks(s).name]);
    plot(stressZ(1,:),zEval,'r', stressZ(2,:),zEval,'g', stressZ(3,:),zEval,'b');
    xlabel('Stress (kN/mm²)'); ylabel('z (mm)'); grid on;
    legend('\sigma_x','\sigma_y','\tau_{xy}');
    title(['Stress Distribution: ' stacks(s).name]);
end

%% 5. Export CSV
for s = 1:2
    writematrix([stacks(s).zEval.' stacks(s).strainZ.'], ['strain_' sName(stacks(s).name) '.csv']);
    writematrix([stacks(s).zEval.' stacks(s).stressZ.'], ['stress_' sName(stacks(s).name) '.csv']);
end

function Qb = QbarCS(Q,theta)
c = cosd(theta); s = sind(theta); c2=c^2; s2=s^2; c4=c2^2; s4=s2^2;
Q11=Q(1,1); Q22=Q(2,2); Q12=Q(1,2); Q66=Q(3,3);
Qb11 = Q11*c4 + 2*(Q12+2*Q66)*c2*s2 + Q22*s4;
Qb22 = Q11*s4 + 2*(Q12+2*Q66)*c2*s2 + Q22*c4;
Qb12 = (Q11+Q22-4*Q66)*c2*s2 + Q12*(c4+s4);
Qb66 = (Q11+Q22-2*Q12-2*Q66)*c2*s2 + Q66*(c4+s4);
Qb16 = (Q11-Q12-2*Q66)*c^3*s - (Q22-Q12-2*Q66)*c*s^3;
Qb26 = (Q11-Q12-2*Q66)*c*s^3 - (Q22-Q12-2*Q66)*c^3*s;
Qb   = [Qb11 Qb12 Qb16; Qb12 Qb22 Qb26; Qb16 Qb26 Qb66];
end

function [A,B,D] = buildABD(Qred,angles,z)
A = zeros(3); B = zeros(3); D = zeros(3);
for k = 1:numel(angles)
    Qb = QbarCS(Qred, angles(k));
    zk = z(k+1); zk1 = z(k);
    A = A + Qb*(zk - zk1);
    B = B + 0.5*Qb*(zk^2 - zk1^2);
    D = D + (1/3)*Qb*(zk^3 - zk1^3);
end
end

function tag = sName(name)
tag = lower(regexprep(name,'[^a-zA-Z0-9]','_'));
end
