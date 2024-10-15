% Call: eme(A,M,L)
% The 1st measure, EME, of enhancement calculation��ǿ����ĵ�һ������EME
% of the image X of size mxn by using blocks of size LxLͨ��ʹ�ô�СΪLxL�Ŀ�Դ�СΪMxM��ͼ��X��������
%
% EME=eme(double(I_OUT),256,5);
% EME1=eme(double(I_f),256,5);
function EME =eme(X,l);

%	L=5; 
[m,n]=size(X);
k1=round(n/l);
k2=round(m/l);
x=round(m/k1-0.5);
y=round(n/k2-0.5);
sum=0;
i_max(:,:)=zeros(k1,k2);
i_min(:,:)=zeros(k1,k2);
t=1.00e-05;
for i=1:1:k1
    for j=1:1:k2
        i_max(i,j)=double(max(max(X(x*(i-1)+1:x*i,y*(j-1)+1:y*j))));
        i_min(i,j)=double(min(min(X(x*(i-1)+1:x*i,y*(j-1)+1:y*j))));
        sum=20*log(i_max(i,j)/(i_min(i,j)+t))+sum;
    end
end
EME=sum/(k1*k2);
end

