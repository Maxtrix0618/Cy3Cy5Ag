% 【调整简并耦合参数D拟合地图 单染料-银膜耦合】
% 对同一批次的每个样品分别单向调整D，制作实验-对照地图

clear; clc;

I=eye(3);
dG = [4 157 107]/255;
TLY = tiledlayout(2,4); % 窗口布局

He=1239.841;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）
L=[498 538];	% 纯染料峰波长 [左 右]
E=He./L;

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
Exp=Exp';
N=length(Exp(:,1));

p=0.02;         % 调整步长
cU=0.6;         % 调整上限
n=(cU/p)+1;
CPR=0:p:cU;

CalS=zeros(N,n,3);
for EI=1:N
    disp(['>>开始组',num2str(EI)]);
    S=He/Exp(EI,1);
    
    for ni=0:(n-1)
        D=p*ni;
        A=[
            S   , D  , D  ;
            D   ,E(1), 0  ;
            D   , 0  ,E(2);
          ];
        Cal=He ./ sort(eig(A),'descend');
        CalS(EI,ni+1,:) = Cal;
        
    end
    disp(['!组',num2str(EI),'计算完成']);
end

% 作图
CB = struct('C', {{'k','r','b'}});
for EI=1:N
    nexttile
    PpR=[450,750];
    DeltaR=[0,cU];
    axis([DeltaR(1) DeltaR(2), PpR(1) PpR(2)]);
    for i=1:2
        plot(DeltaR,[L(i),L(i)], ':','color',dG,'linewidth',1.8);  hold on    % 染料双峰准线
    end
    plot(DeltaR,[Exp(EI,1),Exp(EI,1)], 'm-. ');  hold on    % 银膜准线

    for i=1:3
        plot(CPR,CalS(EI,:,i),strcat(CB.C{i},'o ')); hold on;             % 计算值离散点
        plot(DeltaR,[Exp(EI,i+1),Exp(EI,i+1)],strcat(CB.C{i},'- '),'linewidth',1.2);  hold on    % 实验值准线
    end

    title(['样品',num2str(EI),' | ',num2str(Exp(EI,1)),'nm'])
    xlabel('D(eV)');
    ylabel('Polariton peak(nm)');
end
