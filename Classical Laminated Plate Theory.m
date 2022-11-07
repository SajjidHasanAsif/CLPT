clc
clear all

% Input data 
% Geometry
% Material properties
% Constitutive relations in the laminae coordinate
% Stacking sequence of the laminae
% Constitutive relations in the laminae using Transformation matrix
% Calculations of the unknown field variables
% Calculations of stresses, strains & displacements

% Geometry of the laminate
Lx=input('Length of the plate in X-direction, Lx=');
Ly=input('Length of the plate in Y-direction, Ly=');
H=input('Length of the plate in Y-direction, H=');
N=input('Number of plies in the composite plate, N=');
Npt=input('Number of points in Z-direction of the plate, Npt=');
a=Lx;
b=Ly;

% Thickness of each laminae
sum=0;
for i=1:N+1
    if (i==1)
        fprintf('Please enter the thickness of the laminae %i layer', i-1)
        thk(i)=0;
        sum=sum+thk(i);
        z_lam(i)=-H/2+sum(i)
    else
        thk(i)=input('of plate:=');
        h1(i-1)=thk(i)/2;
        sum=sum+thk(i);
        z_lam(i)=-H/2+sum
    end
end
%h=h(1)';
z_l=z_lam'

% Material properties of the laminae
E1=25;
E2=1;
G12=0.5;
G13=0.5;
G23=0.2;
n12=0.25;
n21=n12*E2/E1;

% Constitutive relations in the laminae coordinate
Q11=E1/(1-n12*n21);
Q12=E2*n12/(1-n12*n21);
Q22=E2/(1-n12*n21);
Q66=G12;
Q44=G23;
Q55=G13;

% Specific coordinate points for computation
XU1=a/2;
XU2=a/2;
XU3=a/2;
YU1=b/2;
YU2=b/2;
YU3=b/2;

XS1=a/2;
XS2=a/2;
XS3=a/2;
YS1=b/2;
YS2=b/2;
YS3=b/2;

XT13=0;
XT23=a/2;
XT12=0;
YT13=b/2;
YT23=0;
YT12=0;

% Load case:
% Sinusoidal load
q_not=1;
q_mn=q_not;

% Ply stack sequence of the laminated plate
sum_Aij=0;
sum_Bij=0;
sum_Dij=0;
m_m=1;
n_n=1;
alpha=m_m.*pi/a;
beta=n_n.*pi/b;
for k=1:N
    i=k;
    k_k=i+1;
    fprintf('Enter the orientation of %i layer',i);
    thta=input('of fibers:=')
    th=pi*thta/180;
    m=cos(th);
    n=sin(th);
    
    % Constitutive relations in the laminae using Transformation matrix
    Q11_bar=Q11*m^4+Q22*n^4+2*(Q12+2*Q66)*m^2*n^2;
    Q12_bar=(Q11+Q22-4*Q66)*m^2*n^2+Q12*(m^4+n^4);
    Q22_bar=Q11*n^4+Q22*m^4+2*(Q12+2*Q66)*m^2*n^2;
    Q16_bar=(Q11-Q22-2*Q66)*m^3*n-(Q22-Q12-2*Q66)*m*n^3;
    Q26_bar=(Q11-Q22-2*Q66)*m*n^3-(Q22-Q12-2*Q66)*m^3*n;
    Q66_bar=(Q11+Q22-2*Q12-2*Q66)*m^2*n^2+Q66*(m^4+n^4);
    Q44_bar=Q44*m^2+Q55*n^2;
    Q45_bar=(Q55-Q44)*m*n;
    Q55_bar=Q55*m^2+Q44*n^2;
    
    Q_mat1=[Q11_bar Q12_bar Q16_bar;...
            Q12_bar Q22_bar Q26_bar;...
            Q16_bar Q26_bar Q66_bar];
    Q_mat2=[Q44_bar Q45_bar ; Q45_bar Q55_bar];
    
    Aij=Q_mat1.*(z_l(k+1)-z_l(k));
    sum_Aij=sum_Aij+Aij;
    Bij=0.5*Q_mat1.*(z_l(k+1)^2-z_l(k)^2);
    sum_Bij=sum_Bij+Bij;
    Dij=(1/3).*Q_mat1.*(z_l(k+1)^3-z_l(k)^3);
    sum_Dij=sum_Dij+Dij
end

A_mat=sum_Aij;
B_mat=sum_Bij;
D_mat=sum_Dij;

A11=A_mat(1,1);
A12=A_mat(1,2);
A16=A_mat(1,3);
A22=A_mat(2,2);
A26=A_mat(2,3);
A66=A_mat(3,3);

B11=B_mat(1,1);
B12=B_mat(1,2);
B16=B_mat(1,3);
B22=B_mat(2,2);
B26=B_mat(2,3);
B66=B_mat(3,3);

D11=D_mat(1,1);
D12=D_mat(1,2);
D16=D_mat(1,3);
D22=D_mat(2,2);
D26=D_mat(2,3);
D66=D_mat(3,3);

% SOLUTION
C11_hat=A11*alpha^2+A66*beta^2;
C12_hat=(A12+A66)*alpha*beta;
C13_hat=-(B11*alpha^2+(B12+2*B66)*beta^2).*alpha;
C22_hat=A66*alpha^2+A22*beta^2;
C23_hat=-(B22*beta^2+(B12+2*B66)*alpha^2)*beta;
C33_hat=(D11*alpha^4+2*(D12+2*D66)*beta^2*alpha^2+D22*beta^4);

C_hat=[C11_hat C12_hat C13_hat;...
       C12_hat C22_hat C23_hat;...
       C13_hat C13_hat C33_hat];
   
F_vec=[0  0  q_not]';
W_1=q_not./C33_hat;
disp_vec=inv(C_hat)*F_vec;

U_mn=disp_vec(1);
V_mn=disp_vec(2);
W_mn=disp_vec(3);

U_0=U_mn*cos(alpha*XS1)*sin(beta*YS1);
V_0=V_mn*sin(alpha*XS2)*cos(beta*YS2);
W_0=W_mn*sin(alpha*XS3)*sin(beta*YS3);

% Anti-Symmetric cross-ply laminate conditions
Q16_bar=0;
Q26_bar=0;
Q45_bar=0;
A16=0;
A26=0;
A45=0;
B16=0;
B26=0;
D16=0;
D26=0;

l1=0;
for pp=1:N
    jj=pp;
    kk=pp;
    fprintf('Enter the orientation of %i layer',pp);
    thta=input('of fibers:=')
    th=pi*thta/180;
    
    z_cord=[linspace(z_l(pp),z_l(pp+1),Npt)]';
    m=cos(th);
    n=sin(th);
    
    Q11_bar=Q11*m^4+Q22*n^4+2*(Q12+2*Q66)*m^2*n^2;
    Q12_bar=(Q11+Q22-4*Q66)*m^2*n^2+Q12*(m^4+n^4);
    Q22_bar=Q11*n^4+Q22*m^4+2*(Q12+2*Q66)*m^2*n^2;
    Q66_bar=(Q11+Q22-2*Q12-2*Q66)*m^2*n^2+Q66*(m^4+n^4);
    Q44_bar=Q44*m^2+Q55*n^2;
    Q55_bar=Q55*m^2+Q44*n^2;
    
    Q_mat1=[Q11_bar Q12_bar Q16_bar;...
            Q12_bar Q22_bar Q26_bar;...
            Q16_bar Q26_bar Q66_bar];
    Q_mat2=[Q44_bar Q45_bar ; Q45_bar Q55_bar];
    
    Ai_mn=(alpha^2.*Q11_bar+beta^2.*Q66_bar).*U_mn+alpha*beta.*(Q12_bar+Q66_bar);
    Ci_mn=alpha*beta.*(Q12_bar+Q66_bar).*U_mn+(alpha^2.*Q66_bar+beta^2.*Q22_bar);
    Bi_mn=-(alpha^3.*Q11_bar+alpha*beta^2.*(Q12_bar+2*Q66_bar)).*W_mn;
    Di_mn=-(beta^3.*Q22_bar+alpha^2*beta.*(Q12_bar+2*Q66_bar)).*W_mn;
    Si_mn=alpha.*Bi_mn+beta.*Di_mn;
    
    % In-plane Displacement
    X_mn=0;
    Y_mn=0;
    
    A_const=(X_mn+z_cord.*U_mn);
    B_const=(Y_mn+z_cord.*V_mn);
    
    u=A_const.*cos(alpha.*XU1).*sin(beta.*YU1);
    v=B_const.*sin(alpha.*XU2).*cos(beta.*YU2);
    w=W_mn.*sin(alpha.*XU3).*sin(beta.*YU3);
    
    % In-plane Strain
    eps_xx=(U_mn.*(-alpha)+z_cord.*alpha^2.*W_mn).*(alpha.*XS1).*sin(beta.*YS1);
    eps_yy=(V_mn.*(-beta)+z_cord.*beta^2.*W_mn).*(alpha.*XS2).*sin(beta.*YS2);
    gam_xy=(U_mn.*beta+V_mn.*alpha+z_cord.*(-2.*alpha*beta.*W_mn)).*cos(alpha.*XT12);
    
    % In-plane Stress
    sigma_xx=Q11_bar.*eps_xx+Q12_bar.*eps_yy;
    sigma_yy=Q12_bar.*eps_xx+Q22_bar.*eps_yy;
    tau_xy=Q66_bar.*gam_xy;
    
    % Out-of-plane stress using the Equilibrium stress
    if pp==1
    tau_xz=((z_cord-z_l(pp)).*Ai_mn+0.5.*(z_cord.^2-z_l(pp).^2).*Bi_mn.*cos(alpha.*XT13).*sin(beta.*YT13));
    tau_yz=((z_cord-z_l(pp)).*Ci_mn+0.5.*(z_cord.^2-z_l(pp).^2).*Di_mn.*sin(alpha.*XT23).*cos(beta.*YT23));
    else
    tau_xz=((z_cord-z_l(pp)).*Ai_mn+0.5.*(z_cord.^2-z_l(pp).^2).*Bi_mn.*cos(alpha.*XT13).*sin(beta.*YT13)+tau_xz(Npt));
    tau_yz=((z_cord-z_l(pp)).*Ci_mn+0.5.*(z_cord.^2-z_l(pp).^2).*Di_mn.*sin(alpha.*XT23).*cos(beta.*YT23)+tau_yz(Npt));
    end
    
    final_z_cord(l1*(Npt)+1:(Npt)*kk,1)=z_cord;
    final_tau_xz(l1*(Npt)+1:(Npt)*pp,1)=tau_xz;
    final_tau_yz(l1*(Npt)+1:(Npt)*pp,1)=tau_yz;
    final_sigma_xx(l1*(Npt)+1:(Npt)*kk,1)=sigma_xx;
    final_sigma_yy(l1*(Npt)+1:(Npt)*kk,1)=sigma_yy;
    final_tau_xy(l1*(Npt)+1:(Npt)*kk,1)=tau_xy;
    
    final_Strain_XX(l1*(Npt)+1:(Npt)*kk,1)=eps_xx;
    final_Strain_YY(l1*(Npt)+1:(Npt)*kk,1)=eps_yy;
    final_GAM_XY(l1*(Npt)+1:(Npt)*kk,1)=gam_xy;
    final_U1(l1*(Npt)+1:(Npt)*kk,1)=u;
    final_V1(l1*(Npt)+1:(Npt)*kk,1)=v;
    final_W1(l1*(Npt)+1:(Npt)*kk,1)=w;
    
    l1=l1+1;
end

tau_xz=final_tau_xz;
tau_yz=final_tau_yz;
tau_xz_B=((H)/(q_not*a)).*tau_xz;
tau_yz_B=((H)/(q_not*a)).*tau_yz;

Sigma_xx=final_sigma_xx;
Sigma_xy=final_tau_xy;
Sigma_yy=final_sigma_yy;

u_disp=final_U1;
v_disp=final_V1;
w_disp=final_W1;

SIG_XXf=Sigma_xx;
TAU_XYf=Sigma_xy;
SIG_YYf=Sigma_yy;

u_b=u_disp*E2*H^3/(q_not*H*a^3);
v_b=v_disp*E2*H^3/(q_not*H*a^3);
w_B=((10*E2*(H)^3)/(q_not*H*a^4))*w_disp;

zcord=final_z_cord;
z_bar=zcord/(H);

TAU_XY_B=(((H)^2)/(q_not*a^2)).*TAU_XYf;
SIG_XX_B=(((H)^2)/(q_not*a^2)).*SIG_XXf;
SIG_YY_B=(((H)^2)/(q_not*a^2)).*SIG_YYf;

MATRIX=[SIG_XX_B SIG_YY_B tau_xz_B tau_yz_B TAU_XY_B z_bar u_b v_b w_B];

vpa(MATRIX);
A=MATRIX;

figure (1)
plot (SIG_XX_B, z_bar,'k','LineWidth',1.5);
title('CLPT');
xlabel('\sigma_{XX}');
ylabel ('$\bar{Z}$');
grid on

figure (2)
plot (SIG_YY_B, z_bar,'r','LineWidth',1.5);
title('CLPT');
xlabel('\sigma_{YY}');
ylabel ('$\bar{Z}$');
grid on

figure (3)
plot (tau_xz_B, z_bar,'g','LineWidth',1.5);
title('CLPT');
xlabel('\sigma_{XZ}');
ylabel ('$\bar{Z}$');
grid on

figure (4)
plot (tau_yz_B, z_bar,'c','LineWidth',1.5);
title('CLPT');
xlabel('\sigma_{YZ}');
ylabel ('$\bar{Z}$');
grid on

figure (5)
plot (TAU_XY_B, z_bar,'m','LineWidth',1.5);
title('CLPT');
xlabel('\sigma_{XY}');
ylabel ('$\bar{Z}$');
grid on

figure (6)
plot (u_b, z_bar,'b','LineWidth',1.5);
title('CLPT');
xlabel('U');
ylabel ('$\bar{Z}$');
grid on

figure (7)
plot (v_b, z_bar,'b','LineWidth',1.5);
title('CLPT');
xlabel('V');
ylabel ('$\bar{Z}$');
grid on

figure (8)
plot (w_B, z_bar,'b','LineWidth',1.5);
title('CLPT');
xlabel('W');
ylabel ('$\bar{Z}$');
grid on