% 【单染料-耦合系数微调地图】

clear; clc;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
L=[512.5 555];  % 纯染料峰波长 [左 右]
E=He./L;


p=10;
Lp=450:p:850;
L0=450-p;
n=(400/p)+1;


Ls=550;
S=He./Ls;

M=zeros(n,n,3);


for Le=Lp
    for Ri=Lp
        disp([num2str(Le),num2str(Ri)])
        
        La=[Le,Ri];
        Ea=He./La;
        
        syms D1 D2;
        H=[
             S  , D1 , D2 ;
             D1 ,E(1),  0 ;
             D2 ,  0 ,E(2)
        ];
        I=eye(3);

        Eq1 = sym(det(H-(Ea(1)*I))==0);
        Eq2 = sym(det(H-(Ea(2)*I))==0);

        Sol = solve(Eq1,Eq2);
        D1a = double(abs(Sol.D1(1)));
        D2a = double(abs(Sol.D2(1)));
        
        A=[
        S   ,D1a ,D2a ;
        D1a ,E(1),  0 ;
        D2a,  0 ,E(2);
        ];
        
        cnl = (Le-L0)/p;
        cnr = (Ri-L0)/p;
        M(cnl,cnr,:) = eig(A);
        
        
    end
end
M = He ./ M;

x=zeros(n,n);
y=zeros(n,n);
for i=1:n
    for j=1:n
        x(i,j)=p*i+L0;
        y(i,j)=p*j+L0;
    end
end


























mesh(x,y,x); hold on;
mesh(x,y,y); hold on;
mesh(x,y,M(:,:,3)); hold on;
mesh(x,y,M(:,:,2)); hold on;

mesh(x,y,M(:,:,1));











