clear fullist;
clear N;
maxbin=300;


N=zeros(1,300);

for k=0:0
    for i=0:0
         for j=0:100
            l=(i^2+j^2+k^2)^(1/2);
            bin=floor((l+.25)*2);
            if bin<maxbin+1 && bin>0
                N(bin)=N(bin)+1;
                fullist(bin,N(bin),1)=i;
                fullist(bin,N(bin),2)=j;
                fullist(bin,N(bin),3)=k;
            end
        end
    end
    for i=1:100
        for j=-100:100
            l=(i^2+j^2+k^2)^(1/2);
            bin=floor((l+.25)*2);
            if bin<maxbin+1 && bin>0
                N(bin)=N(bin)+1;
                fullist(bin,N(bin),1)=i;
                fullist(bin,N(bin),2)=j;
                fullist(bin,N(bin),3)=k;
            end
        end
    end
end

for k=1:100
    for j=-100:100
        for i=-100:100
            l=(i^2+j^2+k^2)^(1/2);
            bin=floor((l+.25)*2);
            if bin<maxbin+1 && bin>0
                N(bin)=N(bin)+1;
                fullist(bin,N(bin),1)=i;
                fullist(bin,N(bin),2)=j;
                fullist(bin,N(bin),3)=k;
            end
        end
    end
end

maxvectors=300;

for i=1:maxbin
   clear X;
   order=randperm(N(i));
   if N(i)>maxvectors
       X(1,:)=fullist(i,order(1:maxvectors),1);
       X(2,:)=fullist(i,order(1:maxvectors),2);
       X(3,:)=fullist(i,order(1:maxvectors),3);
   else
   X(1,:)=fullist(i,order,1);
   X(2,:)=fullist(i,order,2);
   X(3,:)=fullist(i,order,3);
   end
   Y=X';
   filename=strcat('d:\qvectors\qvector_',num2str(i));
   save(filename,'Y','-ascii')
end