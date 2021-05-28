%% --------- Ergodic Capacity for Nt>Nr
% Upper Bound of ECC
% Lower Bound of ECC
% Asymptotic ECC 
%% --------- Simulation Parameters 
Power = [0:6:30]; % transmit power
Monte_Carlo = 1e4; % Monte-Carlo Simulation Times
Nr = 6; % Number of receive antennas
Nt = 9; % Number of transmit antennas
M = 16; % Number of reflecting elements
Mh = 4; % Number of reflecting elements per row
Mv = 4; % Number of reflecting elements per column
Pt = 1; % Number of LoS paths in the transmitter-RIS channel
Pr = 1; % Number of LoS paths in the RIS-receiver channel
Ph = 2; % Number of LoS paths in the transmitter-receiver channel
T_LoS = 0;
for i = [1:1:Pt]
    theta1 = 2*pi*rand(1);
    phi1 = 2*pi*rand(1);
    theta2 = 2*pi*rand(1);
    T_LoS = T_LoS + 1/sqrt(Pt)*kron(exp(1j*[0:1:(Mv-1)]'*pi*sin(theta1)),exp(1j*[0:1:(Mh-1)]'*pi*cos(theta1)*sin(phi1)))*transpose(exp(1j*[0:1:(Nt-1)]'*pi*sin(theta2)));
end
R_LoS = 0;
for i = [1:1:Pr]
    theta1 = 2*pi*rand(1);
    phi1 = 2*pi*rand(1);
    theta2 = 2*pi*rand(1);
    R_LoS = R_LoS + 1/sqrt(Pr)*(exp(1j*[0:1:(Nr-1)]'*pi*sin(theta2)))*transpose(kron(exp(1j*[0:1:(Mv-1)]'*pi*sin(theta1)),exp(1j*[0:1:(Mh-1)]'*pi*cos(theta1)*sin(phi1))));
end
H_LoS = 0;
for i = [1:1:Ph]
    theta1 = 2*pi*rand(1);
    theta2 = 2*pi*rand(1);
    H_LoS = H_LoS + 1/sqrt(Ph)*(exp(1j*[0:1:(Nr-1)]'*pi*sin(theta1)))*transpose(exp(1j*[0:1:(Nt-1)]'*pi*sin(theta2)));
end
Phi = diag(exp(1j*2*pi*rand(M,1))); % Phase shifts used at the RIS
Kr = 1/2; % Rician factor of the RIS-receiver channel
Kh = 1/2; % Rician factor of the transmitter-receiver channel
Exact_value = ones(1,length(Power));
Exact_value1 = ones(1,length(Power));
Capacity_UB = ones(1,length(Power));
Capacity_UB1 = ones(1,length(Power));
Capacity_LB = ones(1,length(Power));
Capacity_LB1 = ones(1,length(Power));
Capacity_As = ones(1,length(Power));
Capacity_As1 = ones(1,length(Power));
for index = [1:1:length(Power)]
    P = 10^(Power(index)/10);
    Exact_value_Monte = ones(1,Monte_Carlo);
    Exact_value_Monte1 = ones(1,Monte_Carlo);
    Capacity_UB_Monte = ones(1,Monte_Carlo);
    Capacity_LB_Monte = ones(1,Monte_Carlo);
    T = T_LoS;
    Psi = Kr/(1+Kr)*T'*T + Kh/(1+Kh)*eye(Nt);
    V1 = svd(Psi);
    Lambda = transpose(V1((Nt-Nr+1):1:Nt));
    for item = [1:1:Monte_Carlo]
        disp(['Power=',num2str(index),',Monte=',num2str(item)]);
        R = sqrt(Kr/(1+Kr))*(1/sqrt(2)*randn(Nr,M)+1j*1/sqrt(2)*randn(Nr,M))+sqrt(1/(1+Kr))*R_LoS;
        H = sqrt(Kh/(1+Kh))*(1/sqrt(2)*randn(Nr,Nt)+1j*1/sqrt(2)*randn(Nr,Nt))+sqrt(1/(1+Kh))*H_LoS;
        T = T_LoS;
        G = R*Phi*T+H;
        Psi = Kr/(1+Kr)*T'*T + Kh/(1+Kh)*eye(Nt);
%         V1 = svd(Psi);
%         Lambda = transpose(V1((Nt-Nr+1):1:Nt));
        G_bar = sqrt(1/(1+Kr))*R_LoS*Phi*T + sqrt(1/(1+Kh))*H_LoS;
        S = 1/sqrt(2)*randn(Nr,Nt)+1j*1/sqrt(2)*randn(Nr,Nt);
        G1 = G_bar*Psi^(-1/2) + S;
        X = 1/sqrt(2)*randn(Nr,Nt)+1j*1/sqrt(2)*randn(Nr,Nt);
        Psi = Kr/(1+Kr)*(T'*T)+Kh/(1+Kh)*eye(Nt);
        Exact_value_Monte(item) = log2(det(eye(Nr)+P/Nt*G*G'));
        Capacity_UB_Monte(item) = det(eye(Nr)+P/Nt*G*G');
%         Capacity_LB_Monte(item) = log(det(G1*Psi*G1'));
        Capacity_LB_Monte(item) = log(det(G1*G1'))+log(prod(Lambda));
        G = G_bar+X*Psi^(1/2);
        Exact_value_Monte1(item) = log2(det(eye(Nr)+P/Nt*G*G'));
    end
    %% Calculate the Upper Bound (Theorem2, eq. (14))
    G_bar = (sqrt(1/(1+Kr))*R_LoS*Phi*T + sqrt(1/(1+Kh))*H_LoS)';
    Psi = Kr/(1+Kr)*T'*T + Kh/(1+Kh)*eye(Nt);
    Tmp = 1;
    q = Nr; 
    for p_index = [1:1:Nr]
        v = combntns([1:1:Nt],p_index);
        p = p_index;
        tmp_det = 0;
        for p_item = [1:1:size(v,1)]
            lambda = real(eig(inv(Psi(v(p_item,:),v(p_item,:)))*(G_bar(v(p_item,:),:)*G_bar(v(p_item,:),:)')));
            theta = lambda(find(abs(lambda)>0.001));
            L = length(theta);
            Delta = eye(L);
            for i = 1:1:L
                for j = 1:1:L
                    Delta(i,j) = (q-L+j+theta(i))*(theta(i))^(j-1);
                end
            end
            V = 1;
            for i = 1:1:L
                for j = (i+1):1:L
                    V = V*(theta(j)-theta(i));
                end
            end
            tmp_det = tmp_det + (gamma(q-L+1)/gamma(q-p+1)*det(Delta)/V*real(det(Psi(v(p_item,:),v(p_item,:)))));
        end
        Tmp = Tmp + tmp_det*(P/Nt)^p_index;
    end   
    %% Calculate the Lower Bound (Theorem4, eq. (28))
    G_bar = sqrt(1/(1+Kr))*R_LoS*Phi*T + sqrt(1/(1+Kh))*H_LoS;
    Psi = Kr/(1+Kr)*T'*T + Kh/(1+Kh)*eye(Nt);
    q = Nt;
    p = Nr;
    lambda = real(eig((G_bar*Psi^(-1/2))*(G_bar*Psi^(-1/2))'));
    theta = lambda(find(abs(lambda)>0.001));
    L = length(theta);
    Delta = 0;
    precision = 150;
    for i = 1:1:L
        Tmp1 = [];
        for j = 1:1:L
            if j~=i
                a = (theta(j)).^([0:1:(L-1)]);
            else
                a = ones(1,L);
                for index1 = [1:1:L]
                    x = theta(j);
                    tmpx = 0;
                    for y = [1:1:precision]
                        tmpx = tmpx + (1-exp(-1*x)*sum((x.^[0:1:(y-1)])./factorial([0:1:(y-1)])))/(q-L+index1+y-1);
                    end

                    a(index1) = x^(index1-1)*tmpx;
                end
            end
            Tmp1 = [Tmp1;a];
        end
        Delta = Delta + det(Tmp1);
    end
    V = 1;
    for i = 1:1:L
        for j = (i+1):1:L
            V = V*(theta(j)-theta(i));
        end
    end
    Q = 0;
    for k = [0:1:(p-1)]
        Q = Q + psi(q-k);
    end
    r2 = Q+Delta/V;
    Psi = Kr/(1+Kr)*T'*T + Kh/(1+Kh)*eye(Nt);
    V1 = svd(Psi);
    lambda = transpose(V1((Nt-Nr+1):1:Nt));
    %% All Simulation Data
    Capacity_UB1(index) = log2(Tmp);
    Capacity_LB1(index) = Nr*log2(1+P/Nt*exp(1/Nr*(log(prod(Lambda))+r2)));
    Exact_value(index) = mean(Exact_value_Monte);
    Exact_value1(index) = mean(Exact_value_Monte1);
    Capacity_UB(index) = log2(mean(Capacity_UB_Monte));
    Capacity_LB(index) = Nr*log2(1+P/Nt*exp(1/Nr*mean(Capacity_LB_Monte)));
    Capacity_As(index) = Nr*(log2(P)-(log2(Nt)-1/log(2)/Nr*mean(Capacity_LB_Monte)));
    Capacity_As1(index) = Nr*(log2(P)-(log2(Nt)-1/log(2)/Nr*(log(prod(Lambda))+r2)));
end
plot(Power,Exact_value,'-ok');hold on;
plot(Power,Exact_value1,'-xk');hold on;
plot(Power,Capacity_UB,'-vk');hold on;
plot(Power,Capacity_UB1,'-sr');hold on;
plot(Power,Capacity_LB,'-sk');hold on;
plot(Power,Capacity_LB1,'-*r');hold on;
plot(Power,Capacity_As,'-k');hold on;
plot(Power,Capacity_As1,'-.r');

