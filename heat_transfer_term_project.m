%------------Heat Transfer Term Project-------------------------------
%--------?—å?å¯?-104303542--?³æ?è«?-104303540---------------------------
clear
clc
close all
%------------Dimensions-----------------------------------------------
Lx=8;                        %length(mm)
Ly=4;                        %width(mm)
dx=2;                        %delta x(mm)
dy=2;                        %delta y(mm)
dx1=1;                       %delta x1(mm)
dy1=1;                       %delta y1(mm)
dx2=2^-1;                    %delta x2(mm)
dy2=2^-1;                    %delta y2(mm)
dx3=2^-2;                    %delta x3(mm)
dy3=2^-2;                    %delta y3(mm)
dx4=2^-3;                    %delta x4(mm)
dy4=2^-3;                    %delta y4(mm)
Nx=Lx/dx;                    %number of nodes in x direction
Ny=Ly/dy;                    %number of nodes in y direction
Nx1=Lx/dx1;                  %number of nodes in x direction(1)
Ny1=Ly/dy1;                  %number of nodes in y direction(1)
Nx2=Lx/dx2;                  %number of nodes in x direction(0.5)
Ny2=Ly/dy2;                  %number of nodes in y direction(0.5)
Nx3=Lx/dx3;                  %number of nodes in x direction(0.25)
Ny3=Ly/dy3;                  %number of nodes in y direction(0.25)
Nx4=Lx/dx4;                  %number of nodes in x direction(0.125)
Ny4=Ly/dy4;                  %number of nodes in y direction(0.125)
%------------Properties-----------------------------------------------
k=10;                        %thermal conductivity(W/(m*K))                      
h=600;                       %convection coefficient(W/(m^2*K))
Bi=(h*dx*10^-3)/k;           %Biot number
Bi1=(h*dx1*10^-3)/k;         %Biot number(1mm)
Bi2=(h*dx2*10^-3)/k;         %Biot number(0.5mm)
Bi3=(h*dx3*10^-3)/k;         %Biot number(0.25mm)
Bi4=(h*dx4*10^-3)/k;         %Biot number(0.125mm)
%------------Boundary Conditions&Ambient Conditions-------------------
Tb=50;                       %temperature of the base (C)
Tinf=25;                     %temperature of the air (C)
%------------Solve Finite Difference Equations------------------------
T=Tinf*ones(Nx+1,Ny+1);      %initialize the temperature matrix
T1=Tinf*ones(Nx1+1,Ny1+1);
T2=Tinf*ones(Nx2+1,Ny2+1);
T3=Tinf*ones(Nx3+1,Ny3+1);
T4=Tinf*ones(Nx4+1,Ny4+1);

for j=1:Ny+1                 %initialize the temperature at base
    T(1,j)=Tb;
end
for o=1:Ny1+1                %initialize the temperature at base(1mm)
    T1(1,o)=Tb;
end
for u=1:Ny2+1                %initialize the temperature at base(0.5mm)
    T2(1,u)=Tb;
end
for b=1:Ny3+1                %initialize the temperature at base(0.25mm)
    T3(1,b)=Tb;
end
for e=1:Ny4+1                %initialize the temperature at base(0.125mm)
    T4(1,e)=Tb;
end
Ttol=10^-5;                  %tolerance temperature difference to stop
nmax=10^9;                   %max time to do iteration

for n=1:nmax
    Told=T;
    for j=2:Ny
        for i=2:(Nx+1)
            if i==(Nx+1)
                %case 3 node at a plane surface with convection
                T(Nx+1,j)=(1/(2*(Bi+2)))*(2*T(Nx,j)+T(Nx+1,j+1)+T(Nx+1,j-1)+2*Bi*Tinf);
            else
                %case 1 interior node
                T(i,j)=(1/4)*(T(i,j+1)+T(i,j-1)+T(i+1,j)+T(i-1,j));
            end
        end
    end
    
    for i=2:(Nx+1)
        if i==Nx+1
            %case 4 node at an external corner with convection
            T(Nx+1,1)=(1/(2*(Bi+1)))*(T(Nx+1,2)+T(Nx,1)+2*Bi*Tinf);
            T(Nx+1,Ny+1)=(1/(2*(Bi+1)))*(T(Nx+1,Ny)+T(Nx,Ny+1)+2*Bi*Tinf);
        else
            %case 3 node at a plane surface with convection
            T(i,Ny+1)=(1/(2*(Bi+2)))*(2*T(i,Ny)+T(i-1,Ny+1)+T(i+1,Ny+1)+2*Bi*Tinf);
            T(i,1)=(1/(2*(Bi+2)))*(2*T(i,2)+T(i-1,1)+T(i+1,1)+2*Bi*Tinf);
        end  
    end
    
    %calculate the max difference between Told and T
    maxdiff=max(max(abs(Told-T)));
    fprintf('Iter = %8.0f - max difference = %10.6f deg. C\n', n, maxdiff);
    if(maxdiff<Ttol)
        break
    end
end
fprintf('Number of iterations = \t %8.0f \n\n', n) % Print how many steps
for p=1:nmax
    Told1=T1;
    for o=2:Ny1
        for r=2:(Nx1+1)
            if r==(Nx1+1)
                %case 3 node at a plane surface with convection
                T1(Nx1+1,o)=(1/(2*(Bi1+2)))*(2*T1(Nx1,o)+T1(Nx1+1,o+1)+T1(Nx1+1,o-1)+2*Bi1*Tinf);
            else
                %case 1 interior node
                T1(r,o)=(1/4)*(T1(r,o+1)+T1(r,o-1)+T1(r+1,o)+T1(r-1,o));
            end
        end
    end
    
    for r=2:(Nx1+1)
        if r==Nx1+1
            %case 4 node at an external corner with convection
            T1(Nx1+1,1)=(1/(2*(Bi1+1)))*(T1(Nx1+1,2)+T1(Nx1,1)+2*Bi1*Tinf);
            T1(Nx1+1,Ny1+1)=(1/(2*(Bi1+1)))*(T1(Nx1+1,Ny1)+T1(Nx1,Ny1+1)+2*Bi1*Tinf);
        else
            %case 3 node at a plane surface with convection
            T1(r,Ny1+1)=(1/(2*(Bi1+2)))*(2*T1(r,Ny1)+T1(r-1,Ny1+1)+T1(r+1,Ny1+1)+2*Bi1*Tinf);
            T1(r,1)=(1/(2*(Bi1+2)))*(2*T1(r,2)+T1(r-1,1)+T1(r+1,1)+2*Bi1*Tinf);
        end  
    end
    
    %calculate the max difference between Told1 and T1
    maxdiff1=max(max(abs(Told1-T1)));
    fprintf('Iter = %8.0f - max difference = %10.6f deg. C\n', p, maxdiff1);
    if(maxdiff1<Ttol)
        break
    end
end
fprintf('Number of iterations = \t %8.0f \n\n', p) % Print how many steps
for v=1:nmax
    Told=T2;
    for u=2:Ny2
        for w=2:(Nx2+1)
            if w==(Nx2+1)
                %case 3 node at a plane surface with convection
                T2(Nx2+1,u)=(1/(2*(Bi2+2)))*(2*T2(Nx2,u)+T2(Nx2+1,u+1)+T2(Nx2+1,u-1)+2*Bi2*Tinf);
            else
                %case 1 interior node
                T2(w,u)=(1/4)*(T2(w,u+1)+T2(w,u-1)+T2(w+1,u)+T2(w-1,u));
            end
        end
    end
    
    for w=2:(Nx2+1)
        if w==Nx2+1
            %case 4 node at an external corner with convection
            T2(Nx2+1,1)=(1/(2*(Bi2+1)))*(T2(Nx2+1,2)+T2(Nx2,1)+2*Bi2*Tinf);
            T2(Nx2+1,Ny2+1)=(1/(2*(Bi2+1)))*(T2(Nx2+1,Ny2)+T2(Nx2,Ny2+1)+2*Bi2*Tinf);
        else
            %case 3 node at a plane surface with convection
            T2(w,Ny2+1)=(1/(2*(Bi2+2)))*(2*T2(w,Ny2)+T2(w-1,Ny2+1)+T2(w+1,Ny2+1)+2*Bi2*Tinf);
            T2(w,1)=(1/(2*(Bi2+2)))*(2*T2(w,2)+T2(w-1,1)+T2(w+1,1)+2*Bi2*Tinf);
        end  
    end
    
    %calculate the max difference between Told2 and T2
    maxdiff2=max(max(abs(Told-T2)));
    fprintf('Iter = %8.0f - max difference = %10.6f deg. C\n', v, maxdiff2);
    if(maxdiff2<Ttol)
        break
    end
end
fprintf('Number of iterations = \t %8.0f \n\n', v) % Print how many steps
for a=1:nmax
    Told=T3;
    for b=2:Ny3
        for c=2:(Nx3+1)
            if c==(Nx3+1)
                %case 3 node at a plane surface with convection
                T3(Nx3+1,b)=(1/(2*(Bi3+2)))*(2*T3(Nx3,b)+T3(Nx3+1,b+1)+T3(Nx3+1,b-1)+2*Bi3*Tinf);
            else
                %case 1 interior node
                T3(c,b)=(1/4)*(T3(c,b+1)+T3(c,b-1)+T3(c+1,b)+T3(c-1,b));
            end
        end
    end
    
    for c=2:(Nx3+1)
        if c==Nx3+1
            %case 4 node at an external corner with convection
            T3(Nx3+1,1)=(1/(2*(Bi3+1)))*(T3(Nx3+1,2)+T3(Nx3,1)+2*Bi3*Tinf);
            T3(Nx3+1,Ny3+1)=(1/(2*(Bi3+1)))*(T3(Nx3+1,Ny3)+T3(Nx3,Ny3+1)+2*Bi3*Tinf);
        else
            %case 3 node at a plane surface with convection
            T3(c,Ny3+1)=(1/(2*(Bi3+2)))*(2*T3(c,Ny3)+T3(c-1,Ny3+1)+T3(c+1,Ny3+1)+2*Bi3*Tinf);
            T3(c,1)=(1/(2*(Bi3+2)))*(2*T3(c,2)+T3(c-1,1)+T3(c+1,1)+2*Bi3*Tinf);
        end  
    end
    
    %calculate the max difference between Told3 and T3
    maxdiff3=max(max(abs(Told-T3)));
    fprintf('Iter = %8.0f - max difference = %10.6f deg. C\n', a, maxdiff3);
    if(maxdiff3<Ttol)
        break
    end
end
fprintf('Number of iterations = \t %8.0f \n\n', a) % Print how many steps
for f=1:nmax
    Told=T4;
    for e=2:Ny4
        for g=2:(Nx4+1)
            if g==(Nx4+1)
                %case 3 node at a plane surface with convection
                T4(Nx4+1,e)=(1/(2*(Bi4+2)))*(2*T4(Nx4,e)+T4(Nx4+1,e+1)+T4(Nx4+1,e-1)+2*Bi4*Tinf);
            else
                %case 1 interior node
                T4(g,e)=(1/4)*(T4(g,e+1)+T4(g,e-1)+T4(g+1,e)+T4(g-1,e));
            end
        end
    end
    
    for g=2:(Nx4+1)
        if g==Nx4+1
            %case 4 node at an external corner with convection
            T4(Nx4+1,1)=(1/(2*(Bi4+1)))*(T4(Nx4+1,2)+T4(Nx4,1)+2*Bi4*Tinf);
            T4(Nx4+1,Ny4+1)=(1/(2*(Bi4+1)))*(T4(Nx4+1,Ny4)+T4(Nx4,Ny4+1)+2*Bi4*Tinf);
        else
            %case 3 node at a plane surface with convection
            T4(g,Ny4+1)=(1/(2*(Bi4+2)))*(2*T4(g,Ny4)+T4(g-1,Ny4+1)+T4(g+1,Ny4+1)+2*Bi4*Tinf);
            T4(g,1)=(1/(2*(Bi4+2)))*(2*T4(g,2)+T4(g-1,1)+T4(g+1,1)+2*Bi4*Tinf);
        end  
    end
    
    %calculate the max difference between Told4 and T4
    maxdiff4=max(max(abs(Told-T4)));
    fprintf('Iter = %8.0f - max difference = %10.6f deg. C\n', f, maxdiff4);
    if(maxdiff4<Ttol)
        break
    end
end
fprintf('Number of iterations = \t %8.0f \n\n', f) % Print how many steps

if (n == nmax) % Warn if number of iterations exceeds maximum 
    fprintf('interation time=%f',kmax) 
    fprintf('\n') 
end
if (p == nmax) % Warn if number of iterations exceeds maximum 
    fprintf('interation time=%f',kmax) 
    fprintf('\n') 
end
if (v == nmax) % Warn if number of iterations exceeds maximum 
    fprintf('interation time=%f',kmax) 
    fprintf('\n') 
end
if (a == nmax) % Warn if number of iterations exceeds maximum 
    fprintf('interation time=%f',kmax) 
    fprintf('\n') 
end
if (f == nmax) % Warn if number of iterations exceeds maximum 
    fprintf('interation time=%f',kmax) 
    fprintf('\n') 
end

%--------calculate the heat rate from the base-----------------------------
qup=h*(dx*(10^-3)/2)*(T(1,Ny+1)-Tinf)+k*((dy*(10^-3))/2)*((T(1,Ny+1)-T(2,Ny+1))/(dx*10^-3));
qdown=h*(dx*(10^-3)/2)*(T(1,1)-Tinf)+k*((dy*(10^-3))/2)*((T(1,1)-T(2,1))/(dx*10^-3));
deltaTmid=0;
for j=2:Ny
    deltaTmid=sum(T(1,j)-T(2,j))+deltaTmid; 
end
qmid=k*(dy*10^-3)*(deltaTmid/(dx*10^-3));
q=qmid+qup+qdown;
fprintf('heat rate per unit length from the base q''=%f (W/m)',q);
fprintf('\n') 
%--------independence test-------------------------------------------------
A=T(5,2);
B=T1(9,3);
C=T2(17,5);
D=T3(33,9);
E=T4(65,17);
%--------figures-----------------------------------------------------------
figure(1)
colormap hot
subplot(2,3,1)
mesh(T);
subplot(2,3,2)
pcolor(T)
subplot(2,3,3)
contour(T)
subplot(2,3,4)
meshz(T)
subplot(2,3,5)
[X,Y] = meshgrid(1:Ny+1,1:Nx+1);
[Tx, Ty] = gradient(T,-1,-1);
contour(T)
hold on
quiver(Tx,Ty)
hold on
sx = 1:Ny+1;
sy = ones(size(sx));
streamline(X,Y,Tx,Ty,sx,sy)
hold off

figure(2)
colormap hot
subplot(2,3,1)
mesh(T1);
subplot(2,3,2)
pcolor(T1)
subplot(2,3,3)
contour(T1)
subplot(2,3,4)
meshz(T1)
subplot(2,3,5)
[X1,Y1] = meshgrid(1:Ny1+1,1:Nx1+1);
[Tx1, Ty1] = gradient(T1,-1,-1);
contour(T1)
hold on
quiver(Tx1,Ty1)
hold on
sx1 = 1:Ny1+1;
sy1 = ones(size(sx1));
streamline(X1,Y1,Tx1,Ty1,sx1,sy1)
hold off
subplot(2,3,6)
line([1 2 3 4 5]',[A B C D E]')