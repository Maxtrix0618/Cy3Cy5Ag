% 【单染料-银膜耦合 全组测试（耦合参数简并）】

clear; clc;
darkGreen = [4 157 107]/255; lw = 1.5;
p = 1;
Lp=450:p:850;
n=(400/p)+1;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
L=[512.5 555];  % 纯染料峰波长 [左 右]

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
Exp=Exp';

TLY = tiledlayout(2,4); % 窗口布局
disp('NSF    Dup    Dmid    Ddn    D_ave        左-实验:计算    中-实验:计算    右-实验:计算  ');

for EI=1:length(Exp(:,1))

    La=Exp(EI,:);

    E=He./L;
    S=He./La(1);
    Ea=He./La(2:4);

    syms D;
    H=[
        S  , D  ,  D ;
        D  ,E(1),  0 ;
        D  , 0  ,E(2)
    ];
    I=eye(3);

    Ds=zeros(3);
    Da = 0;
    for i=1:3
        Sol = solve(sym(det(H-(Ea(i)*I))==0));
        D_ans = double(abs(Sol(1)));

        Ds(i) = D_ans;
        Da = Da + D_ans;
    end
    Da = Da / 3;
    
    Eg=zeros(n,3);
    Ep=zeros(n,3);
    for i=1:n
        S=He/Lp(i);
        A=[
            S   ,Da  ,Da  ;
            Da  ,E(1),  0 ;
            Da  ,  0 ,E(2);
        ];
        Eg(i,:) = eig(A);
        
        for j=1:3
            B=[
                S    ,Ds(j),Ds(j);
                Ds(j),E(1) ,  0  ;
                Ds(j),  0  ,E(2) ;
            ];
            EigB = eig(B);
            Ep(i,j) = EigB(4-j);
        end
        
    end
    Lac = He ./ Eg;
    Lpc = He ./ Ep;
    
    % 打印计算数据
    cL=Exp(EI,1);
    cLN=cL-450+1;
    BK='  '; BC=':';
    disp([num2str(cL),BK,num2str(Ds(1)),BK,num2str(Ds(2)),BK,num2str(Ds(3)),BK,num2str(Da),BK,BK,num2str(Exp(EI,2)),BC,num2str(Lpc(cLN,1)),BK,num2str(Exp(EI,3)),BC,num2str(Lpc(cLN,2)),BK,num2str(Exp(EI,4)),BC,num2str(Lpc(cLN,3))]);

    % 作图
    nexttile;
    axis([450,850 450,850]);
    for i=1:2
        plot([450,850],[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
    end
    plot(Lp,Lp, 'b-. ');
    plot(Lp,Lac(:,1),'r- ',Lp,Lac(:,2),'r- ',Lp,Lac(:,3),'r- ')     % 简并耦合系数的计算结果
    plot(Lp,Lpc(:,1),'r: ',Lp,Lpc(:,2),'r: ',Lp,Lpc(:,3),'r: ')     % 三个耦合系数分别计算结果
    
    % 实验数据
    for i=2:4
        plot(Exp(:,1),Exp(:,i),'k+ '); hold on;
    end

    title(['三峰计算',num2str(EI),'-耦合系数按实验数据组Ag',num2str(Exp(EI,1)),'nm给出'])
    xlabel('Bare plasmon peak(nm)');    ylabel('Polariton peak(nm)');

end






