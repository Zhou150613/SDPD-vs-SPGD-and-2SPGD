%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：基于梯度下降算法算法的Number_DM单元点源自适应光学仿真程序
% 功能：利用各种性能指标作为优化目标函数进行闭环校正，利用峰值Strehl ratio作为评价函数。
% 利用Atm_Turbulence函数生成符合Kolmogorov大气湍流功率谱的畸变波前；
% 服从泊松分布的微小量作为扰动电压；
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

select = 32;

switch(select)
    case 32
        load('A_effect32.mat');%变形镜
%         load('Aberration_WaveByDM32.mat');%畸变波前
        load('Aberration_Wave1028.mat')
        gama = 100;
        Number_Iteration = 800;
    case 156
        load ('A_effect156.mat');
%         load('Aberration_WaveByDM156.mat');%畸变波前
        load('Aberration_Wave1028.mat')
        gama=160;
        Number_Iteration=2000;%每帧畸变相位屏迭代Number_Iteration次
    case 316    
        load ('A_effect316.mat');
        load('Aberration_WaveByDM316.mat');%畸变波前
        gama=150;
        Number_Iteration=4000;%每帧畸变相位屏迭代Number_Iteration次
    case 489
        load ('A_effect489.mat');
        load('Aberration_WaveByDM489.mat');%畸变波前
        gama=80;
        Number_Iteration=8000;%每帧畸变相位屏迭代Number_Iteration次
end



Number_DM=size(A_effect,3);%驱动器个数32  zhou.变形镜
Size_Wave=size(A_effect,2);%波前大小64
Size_Image=1.5*Size_Wave;%96
Size_FFT=3*Size_Wave;%192
Frame=size(Aberration_Wave,3);%共计Frame帧畸变相位屏 1 
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
for nnn=1:1;%共计Frame帧畸变相位屏
    rms(nnn)=std2(Aberration_Wave(:,:,nnn));%初始rms  std()方差
    clf reset;  %清除图像窗口，重置所有对象设置
    dv_i = 1; %变形镜编号
    V_Correct=V_Correct*0;
    Image_Aberration=ComputeFarfield(Aberration_Wave(:,:,nnn),Size_FFT,Size_Image);%用于归一化的平面波及远场
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
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,nnn);
        RMS(nnn,iteration)=std2(Wave_Error);
        
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
%         
        STR(nnn,iteration)=max(Image_Farfield(:));
        
        J(nnn,iteration)=ComputeMR(Image_Farfield);
        
        %%%%%%%%%%%%%%%%%%%%%%%%生成扰动向量%%%%%%%%%%%%%%%%%%
        dv_disturb=zeros(1,Number_DM);% disturb 干扰
        if(dv_i == (Number_DM+1))  %所有变形镜优化完毕后，进行下一轮优化 
            dv_i = 1;
            gama = 0.8*gama;
        end
        dv_disturb(dv_i)=Abs_disturb;
        
        p=sum(dv_disturb)/Number_DM;
%         dv_disturb=dv_disturb-p;
        Vp_temp=V_Correct+dv_disturb;
        Vn_temp=V_Correct-dv_disturb;
        %%%%%%%%%%%%%%%%%%%%%%计算正扰动的性能指标%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+Vp_temp(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,nnn);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_p=ComputeMR(Image_Farfield);

        %%%%%%%%%%%%%%%%%%%%%%计算负扰动的性能指标%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+Vn_temp(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,nnn);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_n=ComputeMR(Image_Farfield);
        
        J_delta = (J_p-J_n);
        
        V_Correct=V_Correct-gama*(J_p-J_n)*dv_disturb*2;
        
        if abs(J_delta)< 0.05
            dv_i = dv_i +1; %优化下一个变形镜
        else    
%             dv_i = dv_i +1; %优化下一个变形镜
        end    
        %%%%%%%%%%%%%%%%%%%%%%计算校正电压%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
%        V_Correct=V_Correct-mean(V_Correct); % mean 计算平均数
    end
end

subplot(2,2,2);
axis([1 Size_Image 1 Size_Image]);
axis image;
imagesc(Image_Farfield);
title('Image After Correction');
drawnow;

subplot(2,2,3);
plot(RMS,'r-');
title('RMS');
legend('定增益梯度下降');
drawnow;
hold on

subplot(2,2,4);
plot(J,'r-');
title('J');
legend('定增益梯度下降');
drawnow;
hold on

save('J_gd.mat',"J","RMS");

STR_GD32 = STR/STR_ideal;

J_GD = J;
RMS_GD = RMS;

load('J_GD_A.mat');load('STR_GD32_A.mat');
J_GD_A = J;STR_GD32_A = STR_GD;
switch(select)
    case 32
        x = 0:1:800-1;
        load('J_SPGD32_average.mat');load('RMS_SPGD32_average.mat');load('STR_SPGD32_average.mat');
        J_SPGD32 = J;RMS_SPGD32 = RMS;STR_SPGD32 = STR;
%         load('J_2SPSA32_average.mat');load('RMS_2SPSA32_average.mat');
%         J_2SPSA32 = J;RMS_2SPSA32 = RMS;        

    case 156

        x = 0:1:2000-1;
        load('J_SPGD156_average.mat');
        load('RMS_SPGD156_average.mat');
        load('STR_SPGD156_average.mat');
        J_SPGD156 = J;RMS_SPGD156 = RMS; STR_SPGD156 = STR;
    case 316
        x = 0:1:4000-1;
        load('J_SPGD316_average.mat');
        load('RMS_SPGD316_average.mat');
        load('STR_SPGD316_average.mat');
        J_SPGD316 = J;RMS_SPGD316 = RMS;STR_SPGD316 = STR;
    case 489
    
end        
        




figure(2);
subplot(2,2,1);
% plot(x,RMS_GD,x,RMS_SPGD32,x,RMS_2SPSA32);

switch(select)
    case 32
        plot(x,RMS_GD,x,RMS_SPGD32);
        title('32维变形镜下 RMS');
        legend('基于PID的串行','SPGD 32维变形镜下1000次平均','2SPSA');
    case 156
        plot(x,RMS_GD,x,RMS_SPGD156);
        title('156维变形镜下 RMS');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        
    case 316
        plot(x,RMS_GD,x,RMS_SPGD316);
        title('316维变形镜下 RMS');
        legend('基于PID的串行','SPGD 316维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 316维变形镜下1000次平均');
    case 489
        title('489维变形镜下 RMS');
        legend('基于PID的串行','SPGD 489维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 489维变形镜下1000次平均');    
end     

subplot(2,2,2);
% plot(x,J_GD,x,J_SPGD32,x,J_2SPSA32);

switch(select)
    case 32
        plot(x,J_GD,x,J_SPGD32,x,J_GD_A);
        title('32维变形镜下 J');
        legend('基于PID的串行','SPGD 32维变形镜下1000次平均',"定增益的串行单维扰动");
    case 156
        plot(x,J_GD,x,J_SPGD156);
        title('156维变形镜下 J');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
    case 316
        plot(x,J_GD,'r--',x,J_SPGD316);
        title('316维变形镜下 J');
        legend('基于PID的串行','SPGD 316维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 316维变形镜下1000次平均');
    case 489
        title('489维变形镜下 J');
        legend('基于PID的串行','SPGD 489维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 489维变形镜下1000次平均');    
end

subplot(2,2,3);
switch(select)
    case 32
        plot(x,STR_GD32,'r--',x,STR_SPGD32,x,STR_GD32_A,'g*');
        title('32维变形镜下 斯特列尔比');
        legend('SDPD','SPGD',"****");
    case 156
        plot(x,STR_GD32,'b--',x,STR_SPGD156);
        title('156维变形镜下 斯特列尔比');
        legend('基于PID的串行单维扰动梯度下降','SPGD 156维变形镜下1000次平均');
        legend('基于PID的串行单维扰动梯度下降','SPGD 156维变形镜下1000次平均');
%         line([1,2000],[0.8,0.8])
    case 316
        
        plot(x,STR_GD32,'b--',x,STR_SPGD316);
        title('316维变形镜下 斯特列尔比');

        legend('基于PID的串行单维扰动梯度下降','SPGD 316维变形镜下1000次平均');
        legend('基于PID的串行单维扰动梯度下降','SPGD 316维变形镜下1000次平均');
    case 489
        title('489维变形镜下 J');
        legend('基于PID的串行单维扰动梯度下降','SPGD 489维变形镜下1000次平均');
        legend('基于PID的串行单维扰动梯度下降','SPGD 489维变形镜下1000次平均');    
end
ylabel('斯特列尔比值');
xlabel('迭代次数');

% close(h);
Mean_STR=mean(STR);
Mean_J=mean(J);


% figure(20);
% x = 0:1:128-1;
% 
% load('D:\我的课题\spsa_and_2spsa\image_SGPD.mat');
% load('D:\我的课题\spsa_and_2spsa\Image_camera_yuanshi.mat');
% load('D:\我的课题\spsa_and_2spsa\image_camera.mat');
% Image_camera_yuanshi = circshift(Image_camera_yuanshi,11)
% h1 = plot(x,image_SGPD,'r');
% hold on
% h2 = plot(x,Image_camera_yuanshi,'g-o');
% h3 = plot(x,image_camera,'b-.');title('波峰');
% ylabel('像素值');
% xlabel('像素位置');
% % plot(x,image_SGPD,'r',x,Image_camera_yuanshi,'g--',x,image_camera,'b-o');
% legend([h1(1),h2(1),h3(3)],'SPGD算法校正','校正前','串行单维扰动梯度下降算法');  % 设置图例


% 
% figure(20);
% x = 0:1:128-1;
% 
% load('C:\Users\zhou\Desktop\spsa_and_2spsa\Image_camera_yuanshi.mat');
% load('C:\Users\zhou\Desktop\spsa_and_2spsa\image_SGPD.mat');
% 
% load('C:\Users\zhou\Desktop\spsa_and_2spsa\image_camera.mat');
% % Image_camera_yuanshi = circshift(Image_camera_yuanshi,11)
% plot(x,sum(Image_camera_yuanshi/4.5),'g-o',x,sum(image_SGPD/4.5),'r',x,sum(image_camera/4.5),'b-.');
% % title('波峰');
% ylabel('');
% xlabel('');
% legend('Before correction','SPGD','SDPD');
% 
% drawnow;

% save('Image_camera_yuanshi.mat')
% save('image_SGPD.mat')
% save('image_camera.mat')

% save('.\\data2\xzMR(2)\xzMR4.mat','Mean_STR','Mean_J','J','STR','Image_Aberration','Image_Farfield','Aberration_Wave','gama');