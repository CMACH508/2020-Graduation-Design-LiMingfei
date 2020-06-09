function  s= sampling(low,up,m,n)
s=[];
t=0;
while 1
 temp=randi([low,up],1);
 if(isempty(find(s==temp)))
 s=[s temp];
 t=t+1;
 end
 if(t>=m*n)
 break;
 end
end
s=reshape(s,m,n);