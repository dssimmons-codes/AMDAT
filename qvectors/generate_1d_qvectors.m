for i=1:300
    X(1)=i/2;
    X(2)=0;
    X(3)=0;
    filename=strcat('d:\qvectors1d\qvector_',num2str(i));
    save(filename,'X','-ascii')
end