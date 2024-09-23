% 【双染料-银膜耦合 全组测试（耦合参数简并）】

clear; clc;
darkGreen = [4 157 107]/255; lw = 1.5;
I=eye(5);
p = 1;
Lp=450:p:850;
n=(400/p)+1;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
L=[525 565 605 660];

fid=fopen('Exp5.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[6,inf]);
Exp=Exp';

TLY = tiledlayout(3,3); % 窗口布局
disp('NSF    DL2    DL1    DMi    DR1    DR2    D_ave');

for EI=1:length(Exp(:,1))

    La=Exp(EI,:);

    E=He./L;
    S=He./La(1);
    Ea=He./La(2:6);

    syms D;
    H=[
         S  , D  , D  , D  , D  ;
         D  ,E(1),  0 ,  0 ,  0 ;
         D  ,  0 ,E(2),  0 ,  0 ;
         D  ,  0 ,  0 ,E(3),  0 ;
         D  ,  0 ,  0 ,  0 ,E(4);
    ];

    Ds=zeros(5);
    Da = 0;
    for i=1:5
        Sol = solve(sym(det(H-(Ea(i)*I))==0));
        D_ans = double(abs(Sol(1)));

        Ds(i) = D_ans;
        Da = Da + D_ans;
    end
    Da = Da / 5;
    
    Eg=zeros(n,5);
    Ep=zeros(n,5);
    for i=1:n
        S=He/Lp(i);
        A=[
             S  , Da , Da , Da , Da ;
             Da ,E(1),  0 ,  0 ,  0 ;
             Da ,  0 ,E(2),  0 ,  0 ;
             Da ,  0 ,  0 ,E(3),  0 ;
             Da ,  0 ,  0 ,  0 ,E(4);
        ];
        Eg(i,:) = eig(A);
        for j=1:5
            B=[
                 S   ,Ds(j),Ds(j),Ds(j),Ds(j);
                 Ds(j),E(1),  0  ,  0 ,   0  ;
                 Ds(j),  0 ,E(2) ,  0 ,   0  ;
                 Ds(j),  0 ,  0  ,E(3) ,  0  ;
                 Ds(j),  0 ,  0  ,  0  ,E(4) ;
            ];
            EigB = eig(B);
            Ep(i,j) = EigB(6-j);
        end

    end
    Lac = He ./ Eg;
    Lpc = He ./ Ep;
    
    % 打印计算数据
    cL=Exp(EI,1);
    cLN=cL-450+1;
    BK='  '; BC=':';
    disp([num2str(cL),BK,num2str(Ds(1)),BK,num2str(Ds(2)),BK,num2str(Ds(3)),BK,num2str(Ds(4)),BK,num2str(Ds(5)),BK,num2str(Da)]);

    % 作图
    nexttile;
    axis([450,850 450,850]);
    for i=1:4
        plot([450,850],[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
    end
    plot(Lp,Lp, 'b-. ');
    plot(Lp,Lac(:,1),'r- ',Lp,Lac(:,2),'r- ',Lp,Lac(:,3),'r- ',Lp,Lac(:,4),'r- ',Lp,Lac(:,5),'r- ')     % 简并耦合系数的计算结果
    plot(Lp,Lpc(:,1),'r: ',Lp,Lpc(:,2),'r: ',Lp,Lpc(:,3),'r: ',Lp,Lpc(:,4),'r: ',Lp,Lpc(:,5),'r: ')     % 三个耦合系数分别计算结果
    
    % 打印实验数据
    for i=2:6
        plot(Exp(:,1),Exp(:,i),'k+ '); hold on;
    end

    title(['耦合系数D=',num2str(Da),'按数据组',num2str(EI),':Ag',num2str(Exp(EI,1)),'nm'])
    xlabel('Bare plasmon peak(nm)');    ylabel('Polariton peak(nm)');

end






