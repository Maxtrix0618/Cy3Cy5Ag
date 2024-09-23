% 【单染料-银膜耦合 全组测试】
% 对一个批次的所有组实验数据分别计算耦合系数D1、D2，然后根据耦合系数做出理论计算的极化峰图线，每组数据做出一张图线。

clear; clc;

I=eye(3);
darkGreen = [4 157 107]/255; lw = 1.5;
p = 1;
Lp=450:p:850;
n=(400/p)+1;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
L=[512.5 555];	% 纯染料峰波长 [左 右]
E=He./L;

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
N=length(Exp(1,:));

TLY = tiledlayout(2,4); % 窗口布局：2*4=8张图
Da=[0,0];

disp('NSF    D1    D2');

for Exp_I=1:N

    La=Exp(:,Exp_I);

    S=He./La(1);
    Ea=He./La(2:4);

    syms D1 D2;
    H=[
         S  , D1 , D2 ;
         D1 ,E(1),  0 ;
         D2 ,  0 ,E(2)
    ];

    Eq1 = sym(det(H-(Ea(1)*I))==0);
    Eq2 = sym(det(H-(Ea(3)*I))==0);

    Sol = solve(Eq1,Eq2);
    D=[double(abs(Sol.D1(1))),double(abs(Sol.D2(1)))];
    Da = Da + D;

    disp([num2str(Exp(1,Exp_I)),'  ',num2str(D(1)),'  ',num2str(D(2))]);
    
    Eg=zeros(n,3);
    for i=1:n
        S=He/Lp(i);
        A=[
        S   ,D(1),D(2);
        D(1),E(1),  0 ;
        D(2),  0 ,E(2);
        ];
        Eg(i,:) = eig(A);
    end
    Lac = He ./ Eg;

    % 作图
    nexttile;
    axis([450,850 450,850]);
    for i=1:2
        plot([450,850],[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
    end
    plot(Lp,Lp, 'b-. ');
    plot(Lp,Lac(:,1),'r- ',Lp,Lac(:,2),'r- ',Lp,Lac(:,3),'r- ')
    
    % 实验数据
    for i=2:4
        plot(Exp(1,:),Exp(i,:),'k+ '); hold on;
    end

    title(['三峰计算-耦合系数按实验数据组：Ag',num2str(Exp(1,Exp_I)),'nm'])
    xlabel('Bare plasmon peak(nm)');
    ylabel('Polariton peak(nm)');

end

Da = Da ./ N;
disp(['均值','  ',num2str(Da(1)),'  ',num2str(Da(2))]);

Eg=zeros(n,3);
for i=1:n
    S=He/Lp(i);
    A=[
    S   ,Da(1),Da(2);
    Da(1),E(1),  0  ;
    Da(2),  0  ,E(2);
    ];
    Eg(i,:) = eig(A);
end
Lac = He ./ Eg;
nexttile;
axis([450,850 450,850]);
for i=1:2
    plot([450,850],[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
end
plot(Lp,Lp, 'b-. ');
plot(Lp,Lac(:,1),'r- ',Lp,Lac(:,2),'r- ',Lp,Lac(:,3),'r- ')
for i=2:4
    plot(Exp(1,:),Exp(i,:),'k+ '); hold on;
end

title('三峰计算-平均耦合系数')
xlabel('Bare plasmon peak(nm)');
ylabel('Polariton peak(nm)');



