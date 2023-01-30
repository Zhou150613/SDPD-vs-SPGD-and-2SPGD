%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：基于经典SPGD算法的Number_DM单元点源自适应光学仿真程序
% 功能：利用各种性能指标作为优化目标函数进行闭环校正，利用峰值Strehl ratio作为评价函数。
% 利用Atm_Turbulence函数生成符合Kolmogorov大气湍流功率谱的畸变波前；
% 服从泊松分布的微小量作为扰动电压；
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% load ('A_effect32.mat');%变形镜 zhou
% load ('A_effect32.mat');
% load ('A_effect52.mat');
% load ('A_effect80.mat');
% load ('A_effect112.mat');
% load ('A_effect156.mat');
% gama=0.2;

% load('A_effect32.mat');
% load('Aberration_WaveByDM32.mat')
% Number_Iteration=800;%每帧畸变相位屏迭代Number_Iteration次
% gama=2.5;


load('A_effect156.mat');
% load('Aberration_WaveByDM156.mat');
load('Aberration_Wave1028.mat');
Number_Iteration=2000;%每帧畸变相位屏迭代Number_Iteration次
gama=0.8;

% load ('A_effect316.mat');
% load('Aberration_WaveByDM316.mat');%畸变波前
% Number_Iteration=4000;%每帧畸变相位屏迭代Number_Iteration次
% gama=0.5;

% load('A_effect489.mat');
% load('Aberration_WaveByDM489.mat');
% Number_Iteration=8000;%每帧畸变相位屏迭代Number_Iteration次
% gama=0.2;


test_num = 100;


Number_DM=size(A_effect,3);%驱动器个数32  zhou.变形镜
Size_Wave=size(A_effect,2);%波前大小64
Size_Image=1.5*Size_Wave;%96
Size_FFT=3*Size_Wave;%192


Frame=size(Aberration_Wave,3)*test_num;%共计Frame帧畸变相位屏 1 
rms=zeros(1,Frame);%0
p_m=ComputeFarfield(zeros(Size_Wave,Size_Wave),Size_FFT,Size_Image);%用于归一化的平面波及远场
p_m_max=1;%max(p_m(:));

Wave_Correct=zeros(Size_Wave,Size_Wave);%变形镜补偿波前
Wave_Error=zeros(Size_Wave,Size_Wave);%残差波前
Image_Farfield=zeros(Size_Image,Size_Image);%校正后远场

J=zeros(Frame,Number_Iteration);%性能指标
RMS=zeros(Frame,Number_Iteration);
STR=zeros(Frame,Number_Iteration);%评价标准

Aberration_Wave_ideal=zeros(Size_Wave,Size_Wave,1);
Image_Farfield=ComputeFarfield(Aberration_Wave_ideal,Size_FFT,Size_Image);
Image_Farfield=Image_Farfield/p_m_max;     
STR_ideal=max(Image_Farfield(:)); %理想情况下，峰值

Abs_disturb=0.2;%扰动电压幅度

V_Correct=zeros(1,Number_DM);% zhou correct 正确的32
dv_disturb=zeros(1,Number_DM);% disturb 干扰
Vp_temp=zeros(1,Number_DM);%temp水温 zhou.正向扰动
Vn_temp=zeros(1,Number_DM);%负向扰动。

%%%%%%%%%%%%%%%%%%%%%%闭环校正%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'程序在运行，请等待...'); %打开进度条
figure(1);
subplot(2,2,1);
for nnn=1:Frame;%共计Frame帧畸变相位屏
    rms(nnn)=std2(Aberration_Wave(:,:,1));%初始rms  std()方差
    clf reset;  %清除图像窗口，重置所有对象设置

    V_Correct=V_Correct*0;
    Image_Aberration=ComputeFarfield(Aberration_Wave(:,:,1),Size_FFT,Size_Image);%用于归一化的平面波及远场
    Image_Aberration=Image_Aberration/p_m_max;
    
    subplot(2,2,1);
    axis([1 Size_Image 1 Size_Image]);
    axis image;
    imagesc(Image_Aberration);
    title('Aberration Image');
    drawnow; %刷新屏幕
    
    for iteration=1:Number_Iteration;%每帧畸变相位屏迭代Number_Iteration次
        waitbar(nnn*iteration/(Number_Iteration*Frame));  %进度条
        %%%%%%%%%%%%%%%%%%%%%%%%保存数据%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+V_Correct(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        RMS(nnn,iteration)=std2(Wave_Error);
        
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
%         
        STR(nnn,iteration)=max(Image_Farfield(:));
        
        J(nnn,iteration)=ComputeMR(Image_Farfield);
        
        %%%%%%%%%%%%%%%%%%%%%%%%生成扰动向量%%%%%%%%%%%%%%%%%%
        for i=1:Number_DM;
            if rand(1)>0.5
                dv_disturb(i)=Abs_disturb;
            else
                dv_disturb(i)=-Abs_disturb;
            end
        end

        p=sum(dv_disturb)/Number_DM;
%         dv_disturb=dv_disturb-p;
        Vp_temp=V_Correct+dv_disturb;
        Vn_temp=V_Correct-dv_disturb;
        %%%%%%%%%%%%%%%%%%%%%%计算正扰动的性能指标%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+Vp_temp(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_p=ComputeMR(Image_Farfield);

        %%%%%%%%%%%%%%%%%%%%%%计算负扰动的性能指标%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+Vn_temp(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_n=ComputeMR(Image_Farfield);

        %%%%%%%%%%%%%%%%%%%%%%计算校正电压%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       V_Correct=V_Correct-gama*(J_p-J_n)*dv_disturb*2;
%        V_Correct=V_Correct-mean(V_Correct); % mean 计算平均数
    end
end
subplot(2,2,2);
axis([1 Size_Image 1 Size_Image]);
axis image;
imagesc(Image_Farfield);
title('Image After Correction');
drawnow;

RMS = mean(RMS,1);

STR = mean(STR,1);
STR = STR/STR_ideal;

subplot(2,2,3);
plot(RMS,'r*');
title('RMS');
drawnow;
hold on

J = mean(J,1);

subplot(2,2,4);
plot(J,'b*');
title('J');
drawnow;
hold on

save('J_SPGD_1028.mat',"J","RMS");

close(h);
Mean_STR=mean(STR);
Mean_J=mean(J);
% save('.\\data2\xzMR(2)\xzMR4.mat','Mean_STR','Mean_J','J','STR','Image_Aberration','Image_Farfield','Aberration_Wave','gama');