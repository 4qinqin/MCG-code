%% ģ������ɢ��
f_org=imread('images/smile_gray.jpg');
f_org=double(f_org);
psf=rand(2350,2350).*(rand(2350,2350)>0.8);
f_diffused =conv2(f_org,psf,'same');
f_diffused=mat2gray(f_diffused);
%v= var(f_diffused(:)) / 10^(5/10);%�� var(f_diffused(:)) / 10^(snr/10)
%f_diffused=imnoise(f_diffused,'gaussian',0,v);
figure(9),imagesc(f_diffused),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;


%% ����ء����㹦����
f_corr=zixiangguan(f_diffused);
f_corr=f_corr(1051:1250,1051:1250);
%f_corr=f_corr(924:1124,1436:1636);%��Ӧʵ�ĳߴ�
f_corr=imadjust(f_corr,[max(min(f_corr)) 1],[0 1]);
figure(10),imagesc(f_corr),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;
S=sqrt(abs(fft2(f_corr)));
B=fftshift(S);%fftshift������Ҷ�任����Ƶ�ʳɷ��Ƶ�Ƶ�����ģ���Ϊfft2�任�У��źŵ���Ƶ�ʳɷ����ź����Ͻǡ�
C=log(1+abs(B));%�Ӷ����Ա�����ʾͼ��
figure(11),imagesc(C),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;

%% ԭͼ����λ�ã����ڽ���Աȣ�
f_obj=single(mat2gray(f_org(1052:1251,1050:1249)));
figure(12),imagesc(f_obj),colormap hot, axis image;colorbar;set(gcf,'color','w');axis off;

%% ����ظ���λ�ָ����ظ������������޸ģ�
iter=1;
mseval=zeros(1,iter);
msevalnew=zeros(1,iter);
i=1;
for i=1:iter
    mseval_temp=zeros(1,10);
    msevalnew_temp=zeros(1,10);
    for j=1:1
        %% �ռ����ʼ�²�
        [Rsize,Csize] = size(S);
        rng('shuffle','twister');
        g2=rand(Rsize,Csize,'single');
        
        %% ��λ�ָ���
        supp=[40 40];%֧��������
         
        [mask,  RFD1, fmse, RfacF1,r3,r4] = phase_rt_opt(f_obj,S,g2,supp,500,300,1,1,0);
        %����Ϊ[sp,  RFD]=phase_rt_opt(f_obj,S,g2,supp,iter,hio_iter,opt,er_pre,showim)
        %RFD=imadjust(RFD,[0.2,0.6],[0,1]);
        mseval_temp(j)=mse(f_obj.*mask,RFD1.*mask);
        ssim_temp = ssim(f_obj.*mask,RFD1.*mask);
        
        
        [mask,  RFD2, fmse, RfacF2,r3_2,r4_2] = phase_rt_opt(f_obj,S,g2,supp,500,300,0,0,0);
        %����Ϊ[sp,  RFD]=phase_rt_opt(f_obj,S,g2,supp,iter,hio_iter,opt,er_pre,showim)
        %RFD=imadjust(RFD,[0.2,0.6],[0,1]);
        msevalnew_temp(j)=mse(f_obj.*mask,RFD2.*mask);
        ssimnew_temp2= ssim(f_obj.*mask,RFD2.*mask);
        
    end
    mseval(i)=min(mseval_temp);
    msevalnew(i)=min(msevalnew_temp);
    display(i);
end

x=60:140;
figure, plot(x,f_obj(60:140,100),'-k',x,r3_2(60:140,100),'g',x,RFD2(60:140,100),'-r',x,r3(60:140,100),'m',x,RFD1(60:140,100),'b');
legend('Ground truth','HIO-ER 300','HIO-ER 500','MCG 300','MCG 500');
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', 'times')
% figure, plot(x,ssimG2,'-r',x,ssimG1,'-b'),axis([0 500 0 1]);legend('HIO-ER','MCG');%set(legend,'Location','NorthEastOutside');
% figure, plot(x,psnrG2,'-r',x,psnrG1,'-b');legend('HIO-ER','MCG');
% figure, plot(x,RfacF2,'-r',x,RfacF1,'-b');legend('HIO-ER','MCG');