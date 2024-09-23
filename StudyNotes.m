%% MATLAB 学习笔记

%% 命令
%%
%'clc' 清空命令行窗口
%'clear' 清除所有变量
%'edit' 新建编辑器脚本
% [选中脚本内容]+F9 执行所选内容
% <Ctrl>+C 中断当前运算
% <Ctrl>+R/T 快速注释/取消注释

% 'ver' 查看版本信息


%% 数据结构
%%
%矩阵
E = zeros(10,5,3)
E(:,:,1) = rand(10,5) 

%元胞数组
A = cell(1,6)
A{2} = eye(3)   %单位矩阵
A{5} = magic(5) %幻方矩阵

%结构体
books = struct('name', {{'Java', 'Python'}}, 'price', [24, 16, 32])
books.name
books.name(1)
books.name{1}
books.price(3)

%矩阵运算
A = [1 2 3 4; 5 6 7 8]  %不同行用分号'；'隔开，同一行中不同矩阵元用空格' '或逗号','隔开
a = A(2,3)  %第2行第3个元素
al = A(2,:) %第2行所有元素
A1 = A(2:7) %逐列读取，从第2个元素到第7个元素，输出为一行
A2 = A(:)   %逐列读取，所有元素，输出为一列
A3 = A'     %转置

B = 1:2:9               %从1到9，步长为2，逐个取值赋给矩阵元
B1 = repmat(B, 3, 2)    %将矩阵B竖向扩展到3倍，横向扩展到2倍

C = [1, 2; 3, 4]
C1 = inv(C)     %逆矩阵
C2 = C * C1     %矩阵相乘

D = zeros(2,4)
D1 = ones(2,4)

E1 = [1 3 5; 7 9 11]
E2 = [0 2 4; 6 8 10]
E3 = E1 * E2'
E4 = E1' * E2
E5 = E1 .* E2   %'.*'表示同位矩阵元相乘

F = rand(5, 8)      %生成5行8列的矩阵，矩阵元均为[0~1]的伪随机数
F1 = rand(5, 8, 'double')   %同上，指定矩阵元为double数据类型
F2 = randn(5, 8)    %生成矩阵元为标准正态分布的伪随机数5x8矩阵
f = randi([2, 64])    %生成[2~64]的伪随机数
F3 = randi([2, 64], 5, 8)   %生成矩阵元为[2~64]的伪随机数5x8矩阵
f1 = find(F3 > 20)      %给出F3中所有大于20的矩阵元的索引（逐列读取）
[m, n] = find(F3 > 20)  %同上，但分成行索引和列索引分别输出


%% 程序结构
%%
%遍历循环(for)
s = 0;
for i = 1:5
    s = s+i^2;
end
s

%条件循环(while)
s = 0;
n = 0;
while s <= 100
    s = (s+1)*2;
    n = n+1;
end
s, n

%条件判断(if...else...)
a = 3;
b = 1;
if a > b
    boo = 'true';
else
    boo = 'false';
end
boo

%条件分支(switch...case...)
a = 2;
str = '无';
switch a
    case 1
        str = '一';
    case 2
        str = '二';
    case 3
        str = '三';
end
str


%% 绘图
%%
% 二维绘图
% Example 1
x = 0 : 0.01 : 2*pi;    %横轴x取值范围[0~2pi]，间隔0.01
y = sin(x);             %纵轴y对于横轴x的函数，正弦函数sin(x)
figure      %建立幕布
plot(x, y)  %绘制图像
title('y=sin(x) 图像')
xlabel('x')
ylabel('y')
xlim([0, 2*pi]) %划定横轴x的显示范围，[0~2pi]

% Example 2
x = 0 : 0.01 : 20;
y1 = 200*exp(-0.05*x) .* sin(x);
y2 = 0.8*exp(-0.5*x) .* sin(10*x);
figure
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');  %AX：轴组，H1：图线1，H2：图线2
title('Multiply Decay Rates')
xlabel('Time(\musec)')                              %设置横轴名字
set(get(AX(1), 'Ylabel'), 'String', 'Slow Decay')   %设置纵轴名字
set(get(AX(2), 'Ylabel'), 'String', 'Fast Decay')
set(H1, 'LineStyle','--')   %设置图线线型
set(H2, 'LineStyle',':')

% 颜色：红r，绿g，蓝b，黄y，粉m，青c，白w，黑k.
% 线型：实线'-'，虚线'--'，点线'：'，点横线'-.'.
% 点型：点·，十字+；圆圈o，星号*，叉号x，正方形s，菱形d，上三角^，下三角v，左三角<，右三角>，五角星p，六角星h.


% 三维绘图
% 曲线（参数方程）
t = 0 : pi/50 : 10*pi;      %自变量轴t取值范围[0~10pi]，间隔pi/50
plot3(sin(t), cos(t), t)    %绘图，前两轴x,y的值是第三轴z的变量t的函数，分别为sin和cos函数
xlabel('sin(t)')
ylabel('cos(t)')
zlabel('t')
axis square     %正方轴（三轴等长显示）。改为'equal'-等长刻度，'tight'-紧凑，'auto'-默认（矩形），'on/off'-开启和关闭轴显示
box off         %关闭上和右边框
grid on         %添加网格


% 曲面（多元函数）
% Eample 1
[x, y, z] = peaks(30);
mesh(x, y, z)   %绘图，无阴影。可改为'surf(x, y, z)'函数以使用阴影风格

% Eample 2
[x, y] = meshgrid(-4 : 0.1 : 4, -4 : 0.1 : 4);
z = cos(x) .* sin(y);
mesh(x, y, z);


%% 文件控制
%% （需先将左侧'当前文件夹目录'定位到文件所在文件夹）

% 写出数据到文件
dataW = rand(6,2);
fileW = fopen('test.txt', 'w');     % 打开文件；'r'-只读， 'w'-只写， 'a'-追加
fprintf(fileW, '%f %f\n', dataW);   % 写入数据
fclose(fileW);                      % 关闭文件

% 从文件读入数据
fileR = fopen('test.txt','r');
dataR = fscanf(fileR, '%d%f', [4,inf]);  % 扫描file内容并读入4行inf（无穷）列的矩阵，每个数据先后以整数'%d'和小数'%f'形式逐列读入矩阵，然后赋值给data
fclose(fileR);


% 直接打印文件
type test.txt



%% END











