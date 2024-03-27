function [sp,  RFD, mseval,RfacF,r3,r4]=phase_rt_opt(f_obj,S,g2,supp,iter,hio_iter,opt,er_pre,showim)
%opt 是否进行复杂度引导 Whether complexity boot is performed
%er_pre 是否含有初始模块，若无则采用进行原ER迭代 Whether the initial module is included. If no, the original ER iteration is used

%% General assignment
filtercount=10;
filtnum=0;
store=0;
toperrs=single(ones(1,10)*1000);
[Rsize,Csize] = size(S);
R2D=zeros(Rsize,Csize,filtercount,'single');
RfacF = zeros(ceil(iter/2),1,'single');  counter1=0; errorF=1;
ssimG = zeros(1,iter,'single');
psnrG = zeros(1,iter,'single');

%% 定义支撑域及其掩膜mask
Rcenter = ceil(Rsize/2);
Ccenter = ceil(Csize/2);
Rsupport = supp(1);
Csupport = supp(2);
half_Rsupport = ceil(Rsupport/2);
half_Csupport = ceil(Csupport/2);
support = zeros(Rsize,Csize,'single');
support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport-1,Ccenter-half_Csupport+1:Ccenter+half_Csupport-1) = 1;
sp=support;

%% 振幅求梯度合
grad_f = zfgradsum(S);
%% 初始ssim值
max_ssim=0;
ssim_g=g2;

%% HIO-ER迭代过程
for iteration = 1:iter
    
    %% Use best result from last filter
    if opt==1 && mod(iteration,ceil(iter/filtercount))==0
        %if toperrs(filtnum)>1
        %    R2D(:,:,filtnum)=mat2gray(R2D(:,:,filtnum))*1e-4;
        %end
         g2=R2D(:,:,filtnum);  
    end
    
    %% iteration with Support & Positivity constraint
    g1=g2;
    phi=angle(fft2(g1));
    %figure(20),imagesc(g1),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
    G=S.*exp(phi*1i);
    g2=real(ifft2(G));
    if er_pre && iteration<=ceil(iter/filtercount)
        g2(g2<0|sp==0) = 0;
        %初始ssim辅助加速
        [ssim_val,~]=ssim(double(abs(fft2(g2))),S);
        if ssim_val>max_ssim 
            max_ssim=ssim_val;
            ssim_g=g2;
        end
        if iteration==ceil(iter/filtercount)
            g2=ssim_g;
        end
    elseif iteration<=hio_iter % && ((~er_pre) || iteration>ceil(iter/filtercount))
        beta=1-0.1*floor(iteration/(hio_iter/10));
        g2(g2<0|sp==0) = g1(g2<0|sp==0) - beta* g2(g2<0|sp==0);
    else
        g2(g2<0|sp==0) = 0;
    end
    
    %% Calculate errors    
    if rem(iteration,1)==0 % && opt==0

        %% Calculate error
        %Ktemp = abs(fft2(g2));
        %errorF = sum(sum(abs(S-Ktemp))) / sum(sum(abs(S)));%mse(Ktemp-S);
        midgrad=gradsum(g2);
        errorF = abs( log(midgrad) -log(grad_f) );
        counter1=counter1+1; RfacF(counter1) = errorF;

        %% Determine iterations with best error
        filtnum=ceil(iteration*filtercount/iter);
        if errorF<= toperrs(filtnum) && iteration>store+1
            toperrs(filtnum)=errorF;
            R2D(:,:,filtnum)=g2;%best result from this filter
            store=iteration;
        end
        
        %% 不同迭代次数恢复图记录
        if iteration==100
            r1=mat2gray(g2);
        elseif iteration==200
            r2=mat2gray(g2);
        elseif iteration==300
            r3=mat2gray(g2);
        elseif iteration==400
            r4=mat2gray(g2);
        elseif iteration==500
            r5=mat2gray(g2);
        end
%         if mod(iteration,1)==0
%             ssimG(iteration/1)=ssim(g2,f_obj);
%             psnrG(iteration/1)=psnr(g2,f_obj);
%         end
             
        %% Figure shows progress
        if showim==1       
            figure (1),
            subplot(1,2,1), imagesc(mat2gray(g2)),colormap hot, axis image, title(int2str(iteration));colorbar;set(gcf,'color','w');axis off;
            subplot(1,2,2), plot(RfacF), axis([0 iteration 0 1]), title(int2str(errorF*100)); 
            drawnow
        end
    end    
    
   %% 用mse值控制迭代次数
   %Ktemp = abs(fft2(g2));
   %mseval_new=mse(Ktemp-S);
%     if iteration>hio_iter% && usemse==1
%         Ktemp = abs(fft2(g2));
%         mseval_new=mse(Ktemp-S);
%         %disp(iteration);disp(mseval-mseval_new);
%         if (iteration>hio_iter+1 && abs(mseval_new)<1e-2)
%             break;
%             iteration
%         end
%         mseval=mseval_new;     
%     end

end

Ktemp = abs(fft2(g2));
mseval=mse(Ktemp-S);
RFD=mat2gray(g2);
%重建的图像
figure(2), imagesc(r1),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
figure(3), imagesc(r2),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
figure(4), imagesc(r3),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
figure(5), imagesc(r4),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
figure(6), imagesc(r5),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
figure(7),imagesc(RFD),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
%绘制数据图表
% figure(8), plot(RfacF), axis([0 ceil(iteration) 0 1]);
% figure(13), plot(ssimG), axis([0 ceil(iteration) 0 1]);
% figure(14), plot(psnrG), axis([0 ceil(iteration) 0 40]);
