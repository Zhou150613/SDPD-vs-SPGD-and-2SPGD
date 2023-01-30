function MR=ComputeMR(far)
SIZE=size(far,2);

CON=zeros(SIZE,SIZE);
for p=1:SIZE
    for q=1:SIZE
        CON(p,q)=p;
    end
end
x=0; y=0; P=0;
P=sum(sum(far));
x=sum(sum(CON.*far));
y=sum(sum((CON').*far));
m_Xc=round(x/P);%SIZE/2;
m_Yc=round(y/P);%SIZE/2;

MR=0;
for p=1:SIZE;
     for q=1:SIZE;
         MR=MR+(sqrt((p-m_Xc)^2+(q-m_Yc)^2))*far(p,q);
     end
end
MR=MR/P;