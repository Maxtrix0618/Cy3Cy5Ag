% 【双染料-银膜耦合 染料交叉项调试-变染料浓度比】
% 考虑染料峰间相互作用，即含非零dye-dye非对角元，生成理论值曲线

clear; clc;

He=1243.125; % 波长和能量换算常数 （nm->eV, E=hc/lamda） 

L=[523 561 621 665];	% 纯染料峰波长 [Cy3左 Cy3右 Cy5左 Cy5右]
E=He./L;                % 各波长对应能量

S=He/610;

D=[0.16139      0.26151      0.15529      0.35152];   % 耦合参数Delta_1~4

p=0.01;         % 步长
Lp=0.01:p:2;    % 横轴：u[Cy3]/u[Cy5]
n=2/p;          % 离散点组数目
C=11;           % 染料交叉项取值数
cp=0.1;         % 染料交叉项取值步长
Eg=zeros(C,n,5);    % C*n*5列的2D数组，记录在一个染料交叉值下第行数个Dye比值的5个矩阵特征值

% 染料交叉项
C12=0;
C13=0;
C14=0;
C23=0;
C24=0;
C34=0;
CC='24';

% 染料交叉项取值遍历
for c=1:C
    C24=(c-1)*cp;
    for x=1:n
        D(3)=D(1)/sqrt(x*p);
        D(4)=D(2)/sqrt(x*p);
        A=[
        S   ,D(1),D(2),D(3),D(4);
        D(1),E(1),C12 ,C13 ,C14 ;
        D(2),C12 ,E(2),C23 ,C24 ;
        D(3),C13 ,C23 ,E(3),C34 ;
        D(4),C14 ,C24 ,C34 ,E(4);
        ];
        Eg(c,x,:) = sort(eig(A),'descend');   % 求A的（5个）特征值（能量）
    end
end
Lac = He ./ Eg;    % 将特征能量转换到特征波长（数组）


% 作图
figure
set(gcf, 'Position', [0, 0, 700, 600]); 
colp=200/C;
cNc=colp*C;
for lc=1:C
    cnc=colp*lc;
    currColor = [255-cnc cnc^2/cNc sqrt(cnc*cNc)]/255;
    for li=1:5
        plot(Lp,Lac(lc,:,li),'-','color',currColor); hold on;
    end
end
xlabel('u[Cy3]/u[Cy5]');
ylabel('Polariton peak(nm)');
title(['Cy3-Cy5@Ag 染料交叉项:C',CC]);


% 作子图
figure
TLY = tiledlayout(2,3); % 窗口布局
set(gcf, 'Position', [0, 0, 1600, 1200]); 
colp=200/C;
cNc=colp*C;
for li=1:5
    nexttile
    for lc=1:C
        cnc=colp*lc;
        currColor = [255-cnc cnc^2/cNc sqrt(cnc*cNc)]/255;
        plot(Lp,Lac(lc,:,li),'-','color',currColor); hold on;
    end
    ylim([300 1000]);
    xlabel('u[Cy3]/u[Cy5]');
    ylabel('Polariton peak(nm)');
    title(['峰:',num2str(li),' | C',CC]);
end
nexttile
for lc=1:C
    cnc=colp*lc;
    currColor = [255-cnc cnc^2/cNc sqrt(cnc*cNc)]/255;
    plot([0,1],[(lc-1)*cp,(lc-1)*cp],'-','color',currColor); hold on;
end
ylabel(['C',CC]);
title(['C',CC,' - 色彩对照']);

% 输出数据
data=zeros(n,C*5+1);
data(:,1)=Lp;
for lc=1:C
    for li=1:5
        data(:,(li-1)*C+lc+1)=Lac(lc,:,li);
%         data(:,(lc-1)*5+li+1)=Lac(lc,:,li);
    end
end
save(['Dye_C',CC,'.txt'], 'data', '-ascii');
