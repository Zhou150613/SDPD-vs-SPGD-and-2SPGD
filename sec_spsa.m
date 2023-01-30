%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 鍚嶇О锛氬熀浜庣粡鍏窼PGD绠楁硶鐨凬umber_DM鍗曞厓鐐规簮鑷?傚簲鍏夊浠跨湡绋嬪簭
% 鍔熻兘锛氬埄鐢ㄥ悇绉嶆?ц兘鎸囨爣浣滀负浼樺寲鐩爣鍑芥暟杩涜闂幆鏍℃锛屽埄鐢ㄥ嘲鍊糞trehl ratio浣滀负璇勪环鍑芥暟銆?
% 鍒╃敤Atm_Turbulence鍑芥暟鐢熸垚绗﹀悎Kolmogorov澶ф皵婀嶆祦鍔熺巼璋辩殑鐣稿彉娉㈠墠锛?
% 鏈嶄粠娉婃澗鍒嗗竷鐨勫井灏忛噺浣滀负鎵板姩鐢靛帇锛?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
                                                                            
select = 80;
test_num = 10;
switch(select)
    case 32
        load ('A_effect32.mat');%鍙樺舰闀? zhou
        load('Aberration_WaveByDM32.mat')
        A = 160;
        a =80;
        Number_Iteration=800;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?    
    case 52
        load ('A_effect52.mat');%鍙樺舰闀? zhou
        load('Aberration_WaveByDM52.mat')
        A = 160;
        a =35;
        Number_Iteration=800;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?    
    case 69
        load ('A_effect69.mat');%鍙樺舰闀? zhou
        load('Aberration_WaveByDM69.mat')
        A = 160;
        a =80;
        Number_Iteration=800;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?            
    case 80
        load ('A_effect80.mat');%鍙樺舰闀? zhou
        load('Aberration_WaveByDM80.mat')
        A = 160;
        a =45;
        Number_Iteration=800;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?            
    case 156
        load('A_effect156.mat');
        load('Aberration_WaveByDM156.mat')
        A = 200;
        a =200; 
        Number_Iteration=2000;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?
    case 316
        load ('A_effect316.mat');
        load('Aberration_WaveByDM316.mat');%鐣稿彉娉㈠墠
        A = 300;
        a =60;
        Number_Iteration=4000;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?
    case 489
        load('A_effect489.mat');
        load('Aberration_WaveByDM489.mat')
        gama=0.1;        
end        
        
% load ('A_effect32.mat');%鍙樺舰闀? zhou
% load ('A_effect52.mat');
% load ('A_effect80.mat');
% load ('A_effect112.mat');
% load ('A_effect156.mat');

Number_DM=size(A_effect,3);%椹卞姩鍣ㄤ釜鏁?32  zhou.鍙樺舰闀?
Size_Wave=size(A_effect,2);%娉㈠墠澶у皬64
Size_Image=1.5*Size_Wave;%96
Size_FFT=3*Size_Wave;%192

Frame=size(Aberration_Wave,3)*test_num;%鍏辫Frame甯х暩鍙樼浉浣嶅睆 1 
rms=zeros(1,Frame);%0
p_m=ComputeFarfield(zeros(Size_Wave,Size_Wave),Size_FFT,Size_Image);%鐢ㄤ簬褰掍竴鍖栫殑骞抽潰娉㈠強杩滃満
p_m_max=1;%max(p_m(:));

Wave_Correct=zeros(Size_Wave,Size_Wave);%鍙樺舰闀滆ˉ鍋挎尝鍓?
Wave_Error=zeros(Size_Wave,Size_Wave);%娈嬪樊娉㈠墠
Image_Farfield=zeros(Size_Image,Size_Image);%鏍℃鍚庤繙鍦?


J=zeros(Frame,Number_Iteration);%鎬ц兘鎸囨爣
RMS=zeros(Frame,Number_Iteration);
STR=zeros(Frame,Number_Iteration);%璇勪环鏍囧噯

Aberration_Wave_ideal=zeros(Size_Wave,Size_Wave,1);
Image_Farfield=ComputeFarfield(Aberration_Wave_ideal,Size_FFT,Size_Image);
Image_Farfield=Image_Farfield/p_m_max;     
STR_ideal=max(Image_Farfield(:)); %理想情况下，峰值


V_Correct_blocking_threshold = 500;


% parameters for gain sequences
% alpha = 0.602; gamma = 0.101;
alpha = 1; gamma = 1;   
c = 45; c_tilde = c;
% a = 230; A = 1000; c = 0.01; c_tilde = c;
w = 0.002; d = 1;
% generate gain, perturbation, weight sequences
% for k = 0
a_0 = a / (1+A)^alpha;
c_0 = c/(1+A); c_tilde_0 = c_tilde/(1+A); w_0 = 1;
% for k = 1,...,n
a_ks = a ./ (((1:Number_Iteration) + 1 + A) .^ alpha); % a_k = a/(k+1+A)^alpha

% c_ks = c ./ (((1:Number_Iteration) + 1 + A) .^ gamma); % c_k = c/(k+1)^gamma
% c_tilde_ks = c_tilde ./ (((1:Number_Iteration) + 1 + A) .^ gamma); % c_tilde_k = c_tilde/(k+1)^gamma
c_ks = 0.1* ones(1,Number_Iteration);
c_tilde_ks= 0.1* ones(1,Number_Iteration);


% parameters space for 
bound = 5000;
V_Correct_max = +bound * ones(Number_DM,1); % lower bounds on theta
V_Correct_min = -bound * ones(Number_DM,1); % upper bounds on theta

H_bar_0 = eye(Number_DM);

Abs_disturb=0.2;%鎵板姩鐢靛帇骞呭害

V_Correct_0=zeros(Number_DM,1);% zhou correct 姝ｇ‘鐨?32
dv_disturb=zeros(Number_DM,1);% disturb 骞叉壈
Vp_temp=zeros(Number_DM,1);%temp姘存俯 zhou.姝ｅ悜鎵板姩
Vn_temp=zeros(Number_DM,1);%璐熷悜鎵板姩銆?


V_Correct=zeros(Number_DM,1);% zhou correct 姝ｇ‘鐨?32

% Q * H * Q' = M * B * M'
% Q: permutation matrix
% H: estimated Hession
% M: lower triangle matrix
% B: block diagonal matrix of size 2
Q_0_vec = 1:Number_DM; % store Q in vector form
change_0 = ones(1,Number_DM); % block structure indicator for B
H_bar_0_factor = eye(Number_DM); % store M and B

%%%%%%%%%%%%%%%%%%%%%%闂幆鏍℃%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'绋嬪簭鍦ㄨ繍琛岋紝璇风瓑寰?...'); %鎵撳紑杩涘害鏉?
figure(1);
subplot(2,2,1);

algo_idx =1;

for nnn=1:Frame;%鍏辫Frame甯х暩鍙樼浉浣嶅睆
    rms(nnn)=std2(Aberration_Wave(:,:,1));%鍒濆rms  std()鏂瑰樊
    clf reset;  %娓呴櫎鍥惧儚绐楀彛锛岄噸缃墍鏈夊璞¤缃?

    V_Correct=V_Correct*0;
    V_Correct_hat_0=V_Correct*0;
    V_Correct_hat_k = V_Correct;
    Image_Aberration=ComputeFarfield(Aberration_Wave(:,:,1),Size_FFT,Size_Image);%鐢ㄤ簬褰掍竴鍖栫殑骞抽潰娉㈠強杩滃満
    Image_Aberration=Image_Aberration/p_m_max;
    
    dv_disturb_0  = 2 * round(rand(Number_DM,1)) - 1;
    p=sum(dv_disturb_0)/Number_DM;
%         dv_disturb=dv_disturb-p;
    vp_hat_0_plus=V_Correct_0+c_0 * dv_disturb_0;
    vn_hat_0_minus=V_Correct_0-c_0 * dv_disturb_0;
    %%%%%%%%%%%%%%%%%%%%%%璁＄畻姝ｆ壈鍔ㄧ殑鎬ц兘鎸囨爣%%%%%%%%%%%%%%%%%%%%%
    Wave_Correct=0*Wave_Correct;
    for i=1:Number_DM
        Wave_Correct=Wave_Correct+vp_hat_0_plus(i)*A_effect(:,:,i);
    end
    Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
    Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
    Image_Farfield=Image_Farfield/p_m_max;

    J_hat_p=ComputeMR(Image_Farfield);

    %%%%%%%%%%%%%%%%%%%%%%璁＄畻璐熸壈鍔ㄧ殑鎬ц兘鎸囨爣%%%%%%%%%%%%%%%%%%%%%
    Wave_Correct=0*Wave_Correct;
    for i=1:Number_DM
        Wave_Correct=Wave_Correct+vn_hat_0_minus(i)*A_effect(:,:,i);
    end
    Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
    Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
    Image_Farfield=Image_Farfield/p_m_max;

    J_hat_n=ComputeMR(Image_Farfield);

    G_hat_0 = (J_hat_p - J_hat_n) / (2 * c_0) * dv_disturb_0;
            % Hessian estimate
    H_bar_k = H_bar_0;
    if (algo_idx == 2) || (algo_idx == 4) % for FA_2SP and FA_E2SP
        Q_k_vec = Q_0_vec;
        H_bar_k_factor = H_bar_0_factor;
        change_k = change_0;
    end
    s_0 = H_bar_0 \ (-G_hat_0);    
            
%% update theta_hat_k
% blocking
    V_Correct_new = V_Correct_hat_0 + a_0 * s_0;
    if norm(V_Correct_new - V_Correct_hat_0) < V_Correct_blocking_threshold
        V_Correct_hat_k = V_Correct_new;
        V_Correct_hat_k = min(V_Correct_hat_k, V_Correct_max);
        V_Correct_hat_k = max(V_Correct_hat_k, V_Correct_min);
    else 
        V_Correct_hat_k = V_Correct_hat_0;
    end
    Wave_Correct=0*Wave_Correct;
    for i=1:Number_DM
        Wave_Correct=Wave_Correct+V_Correct_hat_k(i)*A_effect(:,:,i);
    end
    Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);

    Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
    Image_Farfield=Image_Farfield/p_m_max;

    subplot(2,2,1);
    axis([1 Size_Image 1 Size_Image]);
    axis image;
    imagesc(Image_Aberration);
    title('Aberration Image');
    drawnow; %鍒锋柊灞忓箷
    
    for iteration=1:Number_Iteration;%姣忓抚鐣稿彉鐩镐綅灞忚凯浠umber_Iteration娆?
        waitbar(nnn*iteration/(Number_Iteration*Frame));  %杩涘害鏉?

                    % perturbations
        delta_k = 2 * round(rand(Number_DM, 1)) - 1;
        delta_tilde_k = 2 * round(rand(Number_DM, 1)) - 1;
        % gain sequences
        a_k = a_ks(iteration);
        c_k = c_ks(iteration);
        c_tilde_k = c_tilde_ks(iteration); 
        delta_tilde_k = 2 * round(rand(Number_DM, 1)) - 1;       

                    % adaptive weight
        if (algo_idx == 1) || (algo_idx == 2)
            w_k = w / (iteration + 1).^d;
        elseif (algo_idx == 3) || (algo_idx == 4)
            % optimal weight
%             w_k = c_k^2 * c_tilde_k^2 / ...
%                 (c_0^2*c_tilde_0^2 + sum(c_ks(1:iter_idx).^2 .* c_tilde_ks(1:iter_idx).^2));
            w_k = w / (iteration + 1).^d;
        end

        %% gradient estimates (two-sided)
        % generate pertubations
        vp_hat_k_plus = V_Correct_hat_k + c_k * delta_k;
        vn_hat_k_minus = V_Correct_hat_k - c_k * delta_k;

                    %%%%%%%%%%%%%%%%%%%%%%璁＄畻姝ｆ壈鍔ㄧ殑鎬ц兘鎸囨爣%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+vp_hat_k_plus(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_p=ComputeMR(Image_Farfield);

        %%%%%%%%%%%%%%%%%%%%%%璁＄畻璐熸壈鍔ㄧ殑鎬ц兘鎸囨爣%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+vn_hat_k_minus(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_n=ComputeMR(Image_Farfield);

            % estimate gradient if using 2SP
        G_hat_k = (J_p - J_n) / (2 * c_k) * delta_k;

            %% Hessian estimates (one-sided)
            vp_tilde_k_plus = vp_hat_k_plus + c_tilde_k * delta_tilde_k;
            vn_tilde_k_minus = vn_hat_k_minus + c_tilde_k * delta_tilde_k;

        %%%%%%%%%%%%%%%%%%%%%%璁＄畻姝ｆ壈鍔ㄧ殑鎬ц兘鎸囨爣%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+vp_tilde_k_plus(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_tidle_p=ComputeMR(Image_Farfield);

        %%%%%%%%%%%%%%%%%%%%%%璁＄畻璐熸壈鍔ㄧ殑鎬ц兘鎸囨爣%%%%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+vn_tilde_k_minus(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
        
        J_tilde_n=ComputeMR(Image_Farfield);

        if (algo_idx == 1) || (algo_idx == 2) % SL_2SP or FA_2SP
            d_k = 1 - w_k;
            delta_y_k = (J_tidle_p - J_p) - (J_tilde_n - J_n);
            b_k = w_k * delta_y_k / (4 * c_k * c_tilde_k);
            u_k = delta_tilde_k; v_k = delta_k;
            u_k_norm = norm(u_k); v_k_norm = norm(v_k);
            u_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k + u_k_norm/v_k_norm*v_k);
            v_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k - u_k_norm/v_k_norm*v_k);   
         elseif (algo_idx == 3) || (algo_idx == 4) % SL_E2SP or FA_E2SP
            d_k = 1;
            delta_y_k = (J_tidle_p - J_p) - (J_tilde_n - J_n);
            b_k = w_k * (delta_y_k / (2 * c_k * c_tilde_k) - delta_k'* H_bar_k * delta_tilde_k) / 2;
            u_k = delta_tilde_k; v_k = delta_k;
            u_k_norm = norm(u_k); v_k_norm = norm(v_k);
            u_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k + u_k_norm/v_k_norm*v_k);
            v_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k - u_k_norm/v_k_norm*v_k);         
        end
        
        %% quasi-Newton step
        if (algo_idx == 1) || (algo_idx == 3) % SL_2SP or SL_E2SP
            H_bar_k = d_k * H_bar_k + b_k * (u_tilde_k * u_tilde_k' - v_tilde_k * v_tilde_k');
        elseif (algo_idx == 2) || (algo_idx == 4) % FA_2SP or FA_E2SP
            % Hbar = d_k * H_bar;
            H_bar_k_factor(1:(Number_DM+1):(Number_DM^2)) = d_k * H_bar_k_factor(1:(Number_DM+1):(Number_DM^2));
            block_k_idx = find(change_k == 2);
            H_bar_k_factor(1+(block_k_idx-1)*(Number_DM+1)+1) = d_k * H_bar_k_factor(1+(block_k_idx-1)*(Number_DM+1)+1);
            % Hbar = Hbar + b_k * (uktilde * uktilde')
            [H_bar_k_factor, Q_k_vec, change_k] = SYMUPD(H_bar_k_factor, Q_k_vec, change_k, b_k, u_tilde_k, Number_DM);
            % Hbar = Hbar - b_k * (vktilde * vktilde')
            [H_bar_k_factor, Q_k_vec, change_k] = SYMUPD(H_bar_k_factor, Q_k_vec, change_k, -b_k, v_tilde_k, Number_DM);
            if algo_idx == 4 % FA_E2SP
                H_bar_k = d_k * H_bar_k + b_k * (u_tilde_k * u_tilde_k' - v_tilde_k * v_tilde_k');
            end
        end            
 

        %% modified-Newton step
        if (algo_idx == 1) || (algo_idx == 3) % SL_2SP or SL_E2SP
            % make it positive definite
            H_barbar_k = sqrtm(H_bar_k * H_bar_k' + 10^(-8) * exp(-iteration) * eye(Number_DM));

        elseif (algo_idx == 2) || (algo_idx == 4) % FA_2SP or FA_E2SP
            % get the factorization matrix from H_bar_k_factor  鍒嗚В鐭╅樀   
            % M_k = P_k,B_k = B_k,Vk = Qk
            [M_k, B_k] = RESTORE(H_bar_k_factor, change_k, Number_DM);
            [V_k, Lambda_k] = eig(B_k); % B = V * Lambda * V'   B = Qk *Lambda*Qk'
            eig_vec = diag(Lambda_k);
%             eig_bar_vec = max(eig_vec, 10^(-4));
            eig_bar_vec = max(max(abs(eig_vec), 10^(-8)),  10^(-8) * Number_DM * max(abs(eig_vec)));
        end
        %% descent direction
        if (algo_idx == 1) || (algo_idx == 3) % SL_2SP or SL_E2SP
            s_k = H_barbar_k \ (-G_hat_k);
        elseif (algo_idx == 2) || (algo_idx == 4) % FA_2SP or FA_E2SP
            % descent direction: (Q' * M * V) * D * (V' * M' * Q) * s = -g
            % M * V * D * V' * M' * (Q*s) = -(Q*g)
            Qs_k = M_k \ -G_hat_k(Q_k_vec); % V * Lambda * V' * M' * (Q*s) = inv(M) * (-(Q*g))
            Qs_k = V_k' * Qs_k; % D * V' * M' * (Q*s) = V' * [inv(M) * (-(Q*g))]
            Qs_k = Qs_k ./ eig_bar_vec; % V' * M' * (Q*s) = inv(Lambda) * [V' * inv(M) * (-(Q*g))]
            Qs_k = V_k * Qs_k; % M' * (Q*s) = V * [inv(Lambda) * V' * inv(M) * (-(Q*g))]
            Qs_k = M_k' \ Qs_k; % (Q*s) = inv(M') * [V * inv(Lambda) * V' * inv(M) * (-(Q*g))]
            trans_Q_k_vec = zeros(1,Number_DM);
            trans_Q_k_vec(Q_k_vec) = 1:Number_DM; % vector form of Q'
            s_k = Qs_k(trans_Q_k_vec); % s = Q' * [inv(M') * V * inv(Lambda) * V' * inv(M) * (-(Q*g))]
        end     

        %% update theta_hat_k
        % blocking
        V_Correct_hat_k_new = V_Correct_hat_k + a_k * s_k;
        if norm(V_Correct_hat_k_new - V_Correct_hat_k) < V_Correct_blocking_threshold
            V_Correct_hat_k = V_Correct_hat_k_new;
            V_Correct_hat_k = min(V_Correct_hat_k, V_Correct_max);
            V_Correct_hat_k = max(V_Correct_hat_k, V_Correct_min);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%淇濆瓨鏁版嵁%%%%%%%%%%%%%%%%%%
        Wave_Correct=0*Wave_Correct;
        for i=1:Number_DM
            Wave_Correct=Wave_Correct+V_Correct_hat_k(i)*A_effect(:,:,i);
        end
        Wave_Error=Wave_Correct+Aberration_Wave(:,:,1);
        RMS(nnn,iteration)=std2(Wave_Error);
        
        Image_Farfield=ComputeFarfield(Wave_Error,Size_FFT,Size_Image);
        Image_Farfield=Image_Farfield/p_m_max;
      
%         
        STR(nnn,iteration)=max(Image_Farfield(:));
        
        J(nnn,iteration)=ComputeMR(Image_Farfield);
       

        p=sum(dv_disturb)/Number_DM;


    end
end

STR_2SPSA = STR/STR_ideal;
STR_2SPSA = mean(STR_2SPSA,1);
RMS = mean(RMS,1);
J = mean(J,1);

subplot(2,2,2);
axis([1 Size_Image 1 Size_Image]);
axis image;
imagesc(Image_Farfield);
title('Image After Correction');
drawnow;

subplot(2,2,3);
plot(RMS,'r*');
title('RMS');
drawnow;
hold on

subplot(2,2,4);
plot(J,'g*');
title('J');
drawnow;
hold on

J_2SPSA = J;
RMS_2SPSA = RMS;


switch(select)
    case 32
        x = 0:1:800-1;
        load('J_SPGD32_average.mat');load('RMS_SPGD32_average.mat');load('STR_SPGD32_average.mat');
        J_SPGD32 = J;RMS_SPGD32 = RMS;STR_SPGD32 = STR;
%         load('J_2SPSA32_average.mat');load('RMS_2SPSA32_average.mat');
%         J_2SPSA32 = J;RMS_2SPSA32 = RMS;        
    case 52
        x = 0:1:800-1;
        load('J_SPGD_52_average.mat');load('RMS_SPGD52_average.mat');load('STR_SPGD52_average.mat');
        J_SPGD52 = J;RMS_SPGD52 = RMS;STR_SPGD52 = STR;
    case 69
        x = 0:1:800-1;
        load('J_SPGD_69_average.mat');load('RMS_SPGD69_average.mat');load('STR_SPGD69_average.mat');
        J_SPGD69 = J;RMS_SPGD69 = RMS;STR_SPGD69 = STR;        
    case 80
        x = 0:1:800-1;
        load('J_SPGD_80_average.mat');load('RMS_SPGD80_average.mat');load('STR_SPGD80_average.mat');
        J_SPGD80 = J;RMS_SPGD80 = RMS;STR_SPGD80 = STR;        
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
        plot(x,RMS_2SPSA,x,RMS_SPGD32);
        title('32维变形镜下 RMS');
        legend('基于PID的串行','SPGD 32维变形镜下1000次平均','2SPSA');
    case 52
        plot(x,RMS_2SPSA,x,RMS_SPGD52);
        title('32维变形镜下 RMS');
        legend('基于PID的串行','SPGD 52维变形镜下1000次平均','2SPSA');        
    case 69
        plot(x,RMS_2SPSA,x,RMS_SPGD69);
        title('32维变形镜下 RMS');
        legend('基于PID的串行','SPGD 69维变形镜下1000次平均','2SPSA');             
    case 80
        plot(x,RMS_2SPSA,x,RMS_SPGD80);
        title('32维变形镜下 RMS');
        legend('基于PID的串行','SPGD 52维变形镜下1000次平均','2SPSA');                
    case 156
        plot(x,RMS_2SPSA,x,RMS_SPGD156);
        title('156维变形镜下 RMS');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        
    case 316
        plot(x,RMS_2SPSA,x,RMS_SPGD316);
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
        plot(x,J_2SPSA,x,J_SPGD32);
        title('32维变形镜下 J');
        legend('基于PID的串行','SPGD 32维变形镜下1000次平均');
    case 52
        plot(x,J_2SPSA,x,J_SPGD52);
        title('52维变形镜下 J');
        legend('基于PID的串行','SPGD 52维变形镜下1000次平均');
    case 69
        plot(x,J_2SPSA,x,J_SPGD69);
        title('69维变形镜下 J');
        legend('基于PID的串行','SPGD 69维变形镜下1000次平均');        
    case 80
        plot(x,J_2SPSA,x,J_SPGD80);
        title('80维变形镜下 J');
        legend('基于PID的串行','SPGD 52维变形镜下1000次平均');             
    case 156
        
        plot(x,J_2SPSA,x,J_SPGD156);
        title('156维变形镜下 J');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
    case 316
        plot(x,J_2SPSA,x,J_SPGD316);
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
        plot(x,STR_2SPSA,"--",x,STR_SPGD32);
%         title('32维变形镜下 SR');
        xlabel('Number of Iterations');
        ylabel('SR');
        legend('2SPSA','SPGD',"location","southeast");
    case 52    
        plot(x,STR_2SPSA,"--",x,STR_SPGD52);
%         title('32维变形镜下 SR');
        xlabel('Number of Iterations');
        ylabel('SR');
        legend('2SPSA','SPGD',"location","southeast");       
    case 69    
        plot(x,STR_2SPSA,"--",x,STR_SPGD69);
%         title('32维变形镜下 SR');
        xlabel('Number of Iterations');
        ylabel('SR');
        legend('2SPSA','SPGD',"location","southeast");               
    case 80    
        plot(x,STR_2SPSA,"--",x,STR_SPGD80);
%         title('32维变形镜下 SR');
        xlabel('Number of Iterations');
        ylabel('SR');
        legend('2SPSA','SPGD',"location","southeast");               
    case 156
        plot(x,STR_2SPSA,x,STR_SPGD156);
        hold on;
%         title('156维变形镜下 SR');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 156维变形镜下1000次平均');
%         line([1,2000],[0.8,0.8])
    case 316
        plot(x,STR_2SPSA,x,STR_SPGD316);
%         title('316维变形镜下 SR');
        legend('基于PID的串行','SPGD 316维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 316维变形镜下1000次平均');
    case 489
%         title('489维变形镜下 SR');
        legend('基于PID的串行','SPGD 489维变形镜下1000次平均');
        legend('基于PID的串行','SPGD 489维变形镜下1000次平均');    
end   
save('J_sec_spsa32_average.mat',"J","RMS");

close(h);
Mean_STR=mean(STR);
Mean_J=mean(J);
% save('.\\data2\xzMR(2)\xzMR4.mat','Mean_STR','Mean_J','J','STR','Image_Aberration','Image_Farfield','Aberration_Wave','gama');