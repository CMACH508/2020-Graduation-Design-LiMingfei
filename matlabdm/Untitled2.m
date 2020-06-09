clear 
clc 
eps=0.01;%boundary
beta=200;beta0=5;miu=0.95;rho=100;
coordinates = [
1 565.0 575.0;
2 25.0 185.0;
3 345.0 750.0;
4 945.0 685.0;
5 845.0 655.0;
6 880.0 660.0;
7 25.0 230.0;
8 525.0 1000.0;
9 580.0 1175.0;
10 650.0 1130.0;
11 1605.0 620.0;
12 1220.0 580.0;
13 1465.0 200.0;
14 1530.0 5.0;
15 845.0 680.0;
16 725.0 370.0;
17 145.0 665.0;
18 415.0 635.0;
19 510.0 875.0;
20 560.0 365.0;
21 300.0 465.0;
22 520.0 585.0;
23 480.0 415.0;
24 835.0 625.0;
25 975.0 580.0;
26 1215.0 245.0;
27 1320.0 315.0;
28 1250.0 400.0;
29 660.0 180.0;
30 410.0 250.0;
31 420.0 555.0;
32 575.0 665.0;
33 1150.0 1160.0;
34 700.0 580.0;
35 685.0 595.0;
36 685.0 610.0;
37 770.0 610.0;
38 795.0 645.0;
39 720.0 635.0;
40 760.0 650.0;
41 475.0 960.0;
42 95.0 260.0;
43 875.0 920.0;
44 700.0 500.0;
45 555.0 815.0;
46 830.0 485.0;
47 1170.0 65.0;
48 830.0 610.0;
49 605.0 625.0;
50 595.0 360.0;
51 1340.0 725.0;
52 1740.0 245.0;
];
coordinates(:,1) = []; 
amount = size(coordinates,1); 

%compute_dist_matric
coor_x_tmp1 = coordinates(:,1) * ones(1,amount);
coor_x_tmp2 = coor_x_tmp1';
coor_y_tmp1 = coordinates(:,2) * ones(1,amount);
coor_y_tmp2 = coor_y_tmp1';
dist_matrix = sqrt((coor_x_tmp1 - coor_x_tmp2).^2 + (coor_y_tmp1 - coor_y_tmp2).^2);
vec = ones(1,amount^2);veck= ones(1,amount^2);
row=ones(1,amount);
col=ones(1,amount);

while beta>beta0
	hv=getv(dist_matrix,vec,row,col,beta,amount,rho);
	x=ones(1,amount);
	y=ones(1,amount);
	tempx=zeros(1,amount);tempy=zeros(1,amount);
	while func(dist_matrix,row,col,vec,beta,amount,rho)>=0.001
		for i=1:amount  
			for j=1:amount
				tempx(i)=tempx(i)+hv(amount*(j-1)+i);
				tempy(i)=tempy(i)+hv(amount*(i-1)+j);
			end
			x(i)=row(i)*(tempx(i)-1);
			y(i)=col(i)*(tempy(i)-1);
		end
		row=row+miu*x;
		col=col+miu*y;
	end
	crit=0;
	for i=1:amount^2
		crit=crit+sqrt((hv(i)-vec(i))^2);
	end
	if crit<eps   
		beta = 0.9*beta;  
	else vec=vec+gettheta(dist_matrix,vec,row,col,beta,hv,n,rho)*(hv-vec);
	end
end

while true
	for i=1:amount^2
		if vec(i)>=0.9 
			vec(i)=1;
		else vec(i)=0;
		end
	end 
	flag=1;
	for i=1:amount
        tmp=0;
        for j=1:amount
            tmp=tmp+v(i+(j-1)*n);
		end
        if tmp~=1
            flag=0;
		end
	end
    for j=1:amount
        tmp=0;
        for i=1:amount
            tmp=tmp+v(i+(j-1)*n);
		end
        if tmp~=1
            flag=0;
		end
	end
	if flag==1
		result=0;
		for i=1:amount
			for j=1:amount
				result=result+dist_matrix(i,j)*vec(i+(j-1)*n);
			end
		end
		result%output the lowest cost
	else
		hv=getv(dist_matrix,vec,row,col,beta,amount,rho);
		while func(dist_matrix,row,col,vec,beta,amount,rho)>=0.001
			for i=1:amount  
				for j=1:amount
					tempx(i)=tempx(i)+hv(amount*(j-1)+i);
					tempy(i)=tempy(i)+hv(amount*(i-1)+j);
				end
				x(i)=row(i)*(tempx(i)-1);
				y(i)=col(i)*(tempy(i)-1);
			end
			row=row+miu*x;
			col=col+miu*y;
		end
		crit=0;
		for i=1:amount^2
			crit=crit+sqrt((hv(i)-vec(i))^2);
		end
		if crit<eps			
			rho=rho+2;
		else vec=vec+gettheta(dist_matrix,vec,row,col,beta,hv,n,rho)*(hv-vec);
		end
	end
end
