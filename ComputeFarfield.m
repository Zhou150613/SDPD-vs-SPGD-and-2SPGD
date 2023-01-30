%%%%%%%%%%%%%%%%%%%%FFT计算远场%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%wave为波面，fftsize为傅立叶变换大小，outputsize为截取多大的像面输出%%%%%%%
function Img=ComputeFarfield(wave,fftsize,outputsize)
SIZE=size(wave,2);
zhenfu=ones(SIZE,SIZE);
% zhenfu(1:SIZE/2,1:SIZE/2)=0.01;
% zhenfu(SIZE/2:SIZE,SIZE/2:SIZE)=0.5;
wf=zeros(SIZE,SIZE);
r=SIZE/2+0.5;
for p=1:SIZE;
    x=2*(p-1)/(SIZE-1)-1;
    for q=1:SIZE;
        y=2*(q-1)/(SIZE-1)-1;
        r=sqrt(x*x+y*y);
        if r<=1
            wf(p,q)=zhenfu(p,q).*exp(1i*wave(p,q));%*(-1)^(i+j);
        else
            wf(p,q)=0;
        end
    end
end;
X=fftshift(fft2(wf,fftsize,fftsize));

% temp=abs(X).*abs(X);
% CON=zeros(fftsize,fftsize);
% for p=1:fftsize
%     for q=1:fftsize
%         CON(p,q)=p;
%     end
% end
% x=0; y=0; P=0;
% P=sum(sum(temp));
% x=sum(sum(CON.*temp));
% y=sum(sum((CON').*temp));
% m_Xc=round(x/P);
% m_Yc=round(y/P);
% if(m_Xc-outputsize/2>0&m_Xc+outputsize/2-1<=fftsize&m_Yc-outputsize/2>0&m_Yc+outputsize/2-1<=fftsize)
%     Img=temp(m_Xc-outputsize/2:m_Xc+outputsize/2-1,m_Yc-outputsize/2:m_Yc+outputsize/2-1);
% else
    x=outputsize/2;
    y=fftsize/2;
    for p=1:outputsize;
        for q=1:outputsize;
            Img(p,q)=X(y-(x-p),y-(x-q));
        end
    end
    Img=abs(Img).*abs(Img);
% end