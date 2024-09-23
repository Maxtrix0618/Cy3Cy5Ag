%% 固定总浓度，改变Cy3占比,且C24是变量，为[Cy3]*[Cy5]的函数
%%
function All5_DyeCross_C24goesC3xC5
clear; clc;

He=1243.125; % 波长和能量换算常数 （nm->eV, E=hc/lamda） 

L=[523 561 621 665];	% 纯染料峰波长 [Cy3左 Cy3右 Cy5左 Cy5右]
E=He./L;                % 各波长对应能量

S=He/610;

KD1=@(VK1) VK1; VK1=0;	% D1变化调比
KD3=@(VK3) VK3; VK3=0;	% D3变化调比

D=[0 0 0 0];    % 耦合参数Delta_1~4

p=0.01;         % 步长
Lp=0:p:1;       % 横轴：u[Cy3]/(u[Cy3]+u[Cy5])
n=1/p+1;        % 离散点组数目
Kc=@(VC) VC; VC=0.5;

% 染料交叉项
C12=0; C13=0; C14=0; C23=0 ;C24=0; C34=0;

F = figure(1);
set(gcf, 'Position', [800, 200, 700, 700]); 
LineSet = struct('lineStyle', {{'--', '-'}});
ColorSet = [[0 0 0];[255 0 0];[0 0 255];[24 157 137];[192 0 237]]/255;

H1 = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.12 0.75 0.04],...
    'min' , 0,...
    'max' , 1,...
    'value' , VK1,...
    'callback', @(s, e) @sliderCallback);
H3 = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.07 0.75 0.04],...
    'min' , 0,...
    'max' , 1,...
    'value' , VK3,...
    'callback', @(s, e) @sliderCallback);
HC = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.02 0.75 0.04],...
    'min' , 0,...
    'max' , 20,...
    'value' , VC,...
    'callback', @(s, e) @sliderCallback);

addlistener(H1, 'Value', 'PostSet', @sliderCallback);
addlistener(H3, 'Value', 'PostSet', @sliderCallback);
addlistener(HC, 'Value', 'PostSet', @sliderCallback);

T1 = annotation('textbox',[0.05 0.12 0.8 0.04],'String',['KD1= ', num2str(sprintf('%.2f',VK1))],'EdgeColor','none');
T3 = annotation('textbox',[0.05 0.07 0.8 0.04],'String',['KD3= ', num2str(sprintf('%.2f',VK3))],'EdgeColor','none');
TC = annotation('textbox',[0.05 0.02 0.8 0.04],'String',['Kc= ', num2str(sprintf('%.3f',VC))],'EdgeColor','none');

% 实验数据
fid=fopen('Exp5U.txt','r');
Exp=fscanf(fid,'%f',[6,inf]);
Exp=Exp';
for i=2:6
    for j=1:length(Exp(:,1))
        ev = Exp(j,i);
        if ev >= 0
            plot(Exp(j,1),ev,'+','color',ColorSet(i-1,:)); hold on;
        end
    end
end

P = [];
for lcc=1:2
    for lii=1:5
        P(lcc,lii) = plot(0:0,0:0); hold on;
    end
end
set(gca, 'position', [0.1 0.25 0.85 0.7]);

function sliderCallback(~,~)
    for lc=1:2
        for li=1:5
            delete(P(lc,li));
        end
    end
    delete(T1); delete(T3); delete(TC);
    
    KD1 = get(H1, 'value');
    KD3 = get(H3, 'value');
    Kc = get(HC, 'value');
    Kcs=[0 Kc];
    
    T1 = annotation('textbox',[0.05 0.12 0.8 0.04],'String',['KD1= ', num2str(sprintf('%.2f',KD1))],'EdgeColor','none');
    T3 = annotation('textbox',[0.05 0.07 0.8 0.04],'String',['KD3= ', num2str(sprintf('%.2f',KD3))],'EdgeColor','none');
    TC = annotation('textbox',[0.05 0.02 0.8 0.04],'String',['Kc= ', num2str(sprintf('%.3f',Kc))],'EdgeColor','none');
    
    Eg=zeros(2,n,5);
    for ci=1:2
        for xi=1:n
            x=xi-1;
            u3=x*p;
            u5=1-u3;
            D(1)=KD1*sqrt(u3);
            D(2)=KD1*sqrt(u3);
            D(3)=KD3*sqrt(u5);
            D(4)=KD3*sqrt(u5);
            C24=Kcs(ci)*(u3*u5);
            A=[
                S   ,D(1),D(2),D(3),D(4);
                D(1),E(1),C12 ,C13 ,C14 ;
                D(2),C12 ,E(2),C23 ,C24 ;
                D(3),C13 ,C23 ,E(3),C34 ;
                D(4),C14 ,C24 ,C34 ,E(4);
            ];
            Eg(ci,xi,:) = sort(eig(A),'descend');   % 求A的（5个）特征值（能量）
        end
    end
    Lac = He ./ Eg;     % 将特征能量转换到特征波长（数组）
    
    for lc=1:2
        for li=1:5
            P(lc,li) = plot(Lp,Lac(lc,:,li),LineSet.lineStyle{lc},'color',ColorSet(li,:)); hold on;
        end
    end
    xlabel('u[Cy3] / (u[Cy3]+u[Cy5])');
    ylabel('Polariton peak(nm)');
    title('Cy3-Cy5@Ag | C24=Kc[Cy3][Cy5]');
    
end

end

