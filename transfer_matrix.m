function [X,A,R,T,Phi]= transfer_matrix(theta0, n0, max_lambda, min_lambda,number_layers, n, k, d, pol)
%% n, k are structures containing the refractive index of each layer
%% By Carlos Ríos, 2008, rewritten in August 2015
%% Update March 2017: This function now also calculates the phase change upon reflection. 

%% Initialize calculation
matrix_lambda=max_lambda-min_lambda;

R=zeros(matrix_lambda,1);
T=zeros(matrix_lambda,1);
A=zeros(matrix_lambda,1);
X=zeros(matrix_lambda,1);
Phi=zeros(matrix_lambda,1);

admitanceinvacuum=2.6544e-3;

for r=1:1:matrix_lambda+1

    theta1=asin(n0*sin(theta0)/(n(1).n_val(r)-1i*k(1).k_val(r)));
    y0=n0*admitanceinvacuum;
    y1=(n(1).n_val(r)-1i*k(1).k_val(r))*admitanceinvacuum;

    
    lambda=0.000000001*(r+min_lambda-1);
    y=zeros();
    eta=zeros();
    delta=zeros();
    theta=zeros();    
    M=[1 0; 0 1];
        
    if pol=='TE'
        eta0=y0*cos(theta0);
        eta1=y1*cos(theta1);
    elseif pol=='TM'
        eta0=y0/cos(theta0);
        eta1=y1/cos(theta1);
    else
        fprintf('There is an error, introduce either TE or TM! \n')
        return
    end
          
    for m=2:1:number_layers    
        y(1)=y1;
        eta(1)=eta1;
        theta(1)=theta1;
        theta(m)=asin(n(m-1).n_val(r)*sin(theta(m-1))/(n(m).n_val(r)-1i*k(m).k_val(r)));
        y(m)=(n(m).n_val(r)-1i*k(m).k_val(r))*admitanceinvacuum;  


        if pol=='TE'
            eta(m)=y(m)*cos(theta(m));
        elseif pol=='TM'
            eta(m)=y(m)/cos(theta(m));
        else
        fprintf('There is an error, introduce either TE or TM! \n')
        return
        end

        delta(m-1)=2*pi*d(m-1)*sqrt(n(m-1).n_val(r)^2-k(m-1).k_val(r)^2-n0*(sin(theta0))^2-2*1i*n(m-1).n_val(r)*k(m-1).k_val(r))/(lambda);
        S=[cos(delta(m-1)) (1i*sin(delta(m-1))/eta(m-1)); 1i*eta(m-1)*sin(delta(m-1)) cos(delta(m-1))]; % eta2=eta0

        M=M*S;
        clear S
     end

        M=M*[1; eta(m)];
        C=M(2);
        B=M(1);

        Y=C/B;

        ro=(eta0-Y)/(eta0+Y);
        R(r)=ro*conj(ro);
        T(r)=4*eta0*real(eta(m))/((eta0*B+C)*conj(eta0*B+C));
        T(r)=4*eta0*real(eta(m))/((eta0*B+C)*conj(eta0*B+C));
        A(r)=4*eta0*real(B*conj(C)-eta(m))/((eta0*B+C)*conj(eta0*B+C));
        Phi(r)=atan2(imag(eta(m)*(B*conj(C)-C*conj(B))),real((eta(m)^2)*B*conj(B)-C*conj(C)));
        %Phi(r)=atan(imag(eta(m)*(B*conj(C)-C*conj(B)))/real((eta(m)^2)*B*conj(B)-C*conj(C)));

        A+T+R;

        X(r)=1e9*lambda;


end