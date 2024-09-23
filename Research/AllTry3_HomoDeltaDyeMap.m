% 【单染料-染料微调地图（简并耦合系数-系数由同批次的每组数据分别计算给出）】

clear; clc;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
Exp=Exp';
n=length(Exp(:,1));
I=eye(3);

O=[512.5 555];	% 初始纯染料峰波长 [左 右]
p=2;  R=20;     % 微调步长与单向步数
N=2*R+1;
Lo=[O(1)-p*R, O(1)+p*R];
Ro=[O(2)-p*R, O(2)+p*R];
LeRange = Lo(1) : p : Lo(2);
RiRange = Ro(1) : p : Ro(2);

ErrMap=zeros(N,N,3);

for Le=LeRange
    for Ri=RiRange
        E=He./[Le,Ri];
        
        % 求出同批次的每组实验数据的耦合系数
        DSS=Exp;
        for EI=1:n
            La=Exp(EI,:);
            S=He./La(1);
            Ea=He./La(2:4);
            syms D;
            H=[
                S  , D  ,  D ;
                D  ,E(1),  0 ;
                D  , 0  ,E(2)
            ];
            for i=1:3
                Sol = solve(sym(det(H-(Ea(i)*I))==0));
                D_ans = double(abs(Sol(1)));
                DSS(EI,i+1) = D_ans;
            end
        end
        
        % 根据所得耦合系数计算理论峰线
        Eg=zeros(n,4);
        for line=1:3
            for i=1:n
                S=He/Exp(i,1);
                Da=DSS(n,line+1);
                A=[
                    S  ,Da  ,Da  ;
                    Da ,E(1),  0 ;
                    Da ,  0 ,E(2);
                ];
                EigA = eig(A);
                Eg(i,line+1) = EigA(4-line);
            end
        end
        Lac = He ./ Eg;
        Lac(:,1)=Exp(:,1);

        % 给出三条计算峰线与实验数据的标准差
        le = (Le-Lo(1))/p +1;
        ri = (Ri-Ro(1))/p +1;
        for line=1:3
            for i=1:n
                T = (Lac(i,line+1)-Exp(i,line+1))^2;
                ErrMap(le,ri,line) = ErrMap(le,ri,line) + T;
            end
        end
        
        disp([le,ri]);
    end
end

ErrMap = ErrMap./n;

ErrMap = sqrt(ErrMap);


x=zeros(N,N);
y=zeros(N,N);
for i=1:N
    for j=1:N
        x(i,j)=p*i+Lo(1);
        y(i,j)=p*j+Ro(1);
    end
end

Mid=floor((N+1)/2);
for line=1:3
    surf(x,y,ErrMap(:,:,line)); hold on;
    disp([num2str(line),'  ',num2str(ErrMap(Mid,Mid,line))])
end
plot3([O(1),O(1)],[O(2),O(2)],[0,40],'-.','Color',[0.4,0.4,0.4]);

title('微调染料峰：极化峰实验与计算的标准差')
xlabel('染料左峰');
ylabel('染料右峰');



