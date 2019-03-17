% drilling lost circulation data

% fileid = fopen('DS_lost_data.txt');
% C= textscan(fileid, '%f %f');
% fclose(fileid);


% The basic parameter of the well

[x, y] = textread('DS_lost_data.txt');

% be aware the x is time
% be aware that y is the cummulative lost rate...

%%%%%%%%%%%%      Drilling data from the field
%%%%%%%%%%%%           Inversion problem
%%%%%%%%%%%%%

rw = 311.2/2 / 1000;  % the radius of the well
t_y = 5.62;           % yield stress
m = 0.8;              % flow index
delta_p = 4.2 * 10^6; % over pressure
k = 1.85;             % flow index

% The operation type curve
X = y / (pi*rw^2);
X = log(X);

Y = x*(m/(2*m+1))*(1/rw)^((1+m)/m)*(delta_p/k)^(1/m);
Y = log(Y);
figure(2)
plot(X,Y,'Ob')

%%%%%%% Linear regression for Field data
%%%%%%%
%%%%%%% Linear regression

[m,n] = size(X);
Z = [ones(m,1),X];
theta = pinv(Z'*Z)*Z'*Y;
hold on

% plot the pic of drilling data
p1 = polyfit(X,Y,1);
plot(X,polyval(p1,X))
slope_drilling = p1(1,1);

% theoretical analysis from the physical model
% theoretical analysis from the physical model
% plot theoretical data from the 
figure (3);
m = 0.65;
a =[0.001,0.008,0.001 , 0.004, 0.01 , 0.02, 0.1 , 0.08];
[b,v] = size(a);
slope_theoretical = zeros(v,1);

for j = 1: v
    y = 1.1 : (1+1/a(j));
    n = floor ((1+1/a(j))) - 1;
    % count how many data points we have
    I = zeros(n,1);
    J = zeros(n,1);
    
    % this solution for large order might not be accurate!
    % we use simpson integrition
    
        for i = 1:n
        f=@(x) 2.^((m+1)./m).* x .* ((x.^(1-m)-1)./(1-m)).^(1./m) ./ (1 - a(j).*(x -1)).^(1./m);
        I(i,1) = simpsons( f,1,y(i), 10000);    
        %%% we can not use Simpson here, the equally step will have HUGE error
        J(i,1) = y(i)^2 - 1;
        end
        
    plot (log(I),log(J));
    c = polyfit(log(I),log(J),1);
    slope_theoretical(j) = c(1,1);
    hold on
end

figure(3)

%%% 

% curve fit for the final version

p1 = polyfit(log(I),log(J),1);
plot(log(I),polyval(p1,log(I)),'O')
