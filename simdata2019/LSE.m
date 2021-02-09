clear all;
%% Prepare dataset
Files = ["LSE_de3211.mat" ; "LSE_da3211.mat" ; "LSE_dadoublet.mat" ; "LSE_drdoublet.mat"; "LSE_dr3211.mat"];
[a, ~] = size(Files);
Z = [];
U = [];
X = [];
De = [];
Da = [];
Dr = [];
Tc = [];
alpha2 = [];
datasets_evaluated_num = -1;

for i=1:5
    load(Files(i));
    Z = [Z, Z_k1k1(:,2:end)];
    U = [U, U_k(:,2:end)];
    X = [X, XX_k1k1];
    De = [De; de(2:end)];
    Da = [Da; da(2:end)];
    Dr = [Dr; dr(2:end)];
    Tc = [Tc; Tc1(2:end)];
    datasets_evaluated_num = datasets_evaluated_num +1;
    clear Z_k1k1 U_k XX_k1k de da dr Tc1;
end
%Make it work for the old code
Z_k1k1 = Z;
U_k = U;
XX_k1k1 = X;
de = De;
da = Da;
dr = Dr;
Tc1 = Tc;
%% Given Aircraft parameters
mass = 4500; %Kg
Ixx = 11187.8; %kg m2
Iyy = 22854.8; %Kg m2
Izz = 31974.8; %Kg m2
Ixz = 19301.1; %Kg m2
b = 13.3250; %Wingspan in meters
S = 24.9900; %Wing area in m2
chord = 1.9910; %Chord
N = size(XX_k1k1,2);
ignore_data =1;

rho =1.14; %TODO: Change this appropriately
Cx = zeros(N,1);
Cy = zeros(N,1);
Cz = zeros(N,1);
Cl = zeros(N,1);
Cm = zeros(N,1);
Cn = zeros(N,1);
p_old = 0;
q_old = 0;
r_old = 0;
p_new = 0;
q_new = 0;
r_new = 0;
p_dot = 0;
q_dot = 0;
r_dot = 0;

XX_10 = mean(XX_k1k1(10,:));
XX_11 = mean(XX_k1k1(11,:));
XX_12 = mean(XX_k1k1(12,:));
XX_13 = mean(XX_k1k1(13,:));
XX_14 = mean(XX_k1k1(14,:));
XX_15 = mean(XX_k1k1(15,:));

for i=1:N-1
if i>ignore_data
        denom_tmp = (0.5*rho*Z_k1k1(10,i)^2*S);
        Axx(i,1) = U_k(1,i) - XX_10;
        Ayx(i,1) = U_k(2,i) - XX_11;
        Azx(i,1) = U_k(3,i) - XX_12;
        p_new(i,1) = U_k(4,i+1) - XX_13;
        q_new(i,1) = U_k(5,i+1) - XX_14;
        r_new(i,1) = U_k(6,i+1) - XX_15;
        p_old(i,1) = U_k(4,i-1) - XX_13;
        q_old(i,1) = U_k(5,i-1) - XX_14;
        r_old(i,1) = U_k(6,i-1) - XX_15;
        p_dot(i,1) = (p_new(i,1) - p_old(i,1))/(2*dt);
        q_dot(i,1) = (q_new(i,1) - q_old(i,1))/(2*dt);
        r_dot(i,1) = (r_new(i,1) - r_old(i,1))/(2*dt);
        Cx(i,1) = mass*Axx(i,1)/denom_tmp;
        Cy(i,1) = mass*Ayx(i,1)/denom_tmp;
        Cz(i,1) = mass*Azx(i,1)/denom_tmp;
        Cl(i,1) = (p_dot(i,1)*Ixx + (q_new(i,1)*r_new(i,1))*(Izz-Iyy) - (p_new(i,1)*r_new(i,1) + r_dot(i,1))*Ixz)/(denom_tmp * b);
        Cm(i,1) = (q_dot(i,1)*Iyy + (r_new(i,1)*p_new(i,1))*(Ixx-Izz) + (p_new(i,1)*p_new(i,1) - r_new(i,1)*r_new(i,1))*Ixz)/(denom_tmp * chord);
        Cn(i,1) = (r_dot(i,1)*Izz + (p_new(i,1)*q_new(i,1))*(Iyy-Ixx) + (q_new(i,1)*r_new(i,1) - p_dot(i,1))*Ixz)/(denom_tmp * b);
end
end
%% Choose Evluation data
Percent_excited_data = 0.7;
data_evaluation = [];
N_perdata =12000;
Eval_start = 900;
Eval_end = 1800;
for i =0:(datasets_evaluated_num)
    Excited_data_start = Eval_start + i*(N_perdata);
    Excited_data_end = Eval_end + i*(N_perdata);
    data_evaluation = [data_evaluation Excited_data_start+randperm(Excited_data_end-Excited_data_start,Percent_excited_data*(Excited_data_end-Excited_data_start))];  
    %Percent_non_excited_data = 0.1;
    %data_evaluation = [data_evaluation Excited_data_end+randperm(N_perdata+i*(N_perdata)-Excited_data_end,Percent_non_excited_data*(Excited_data_end-Excited_data_start))];
end
size_data_eval = size(data_evaluation,2);
%% Choose validation data
data_validation = [];
N_perdata =12000;
for i =0:(datasets_evaluated_num)
    Excited_data_start = Eval_start + i*(N_perdata);
    Excited_data_end = Eval_end + i*(N_perdata);
    Excited_data = Excited_data_start:Excited_data_end;
    data_validation = [data_validation setdiff(Excited_data,data_evaluation)];  
    %Percent_non_excited_data = 0.1;
    %data_validation = [data_validation Excited_data_end+randperm(N_perdata+i*(N_perdata)-Excited_data_end,Percent_non_excited_data*(Excited_data_end-Excited_data_start))];
end
size_data_val = size(data_validation,2);
%% Cx estimate coefficients
ACx = ones(6,size_data_eval);
ACx(2,:) = (Z_k1k1(11,data_evaluation)); %alpha
ACx(3,:) = (Z_k1k1(11,data_evaluation)).^2; %alpha^2
ACx(4,:) = (U_k(5,data_evaluation) - XX_k1k1(14,data_evaluation)).*chord./Z_k1k1(10,data_evaluation); %q
ACx(5,:) = de(data_evaluation); %de
ACx(6,:) = Tc1(data_evaluation); %throttle
max_ACx = max(abs(ACx),[],2);
in = find(max_ACx == 0);
max_ACx(in) = 1;
ACx = ACx./max_ACx;
max_Cx = max(abs(Cx));
ParamCx = (ACx*ACx')\ACx*Cx(data_evaluation)/max_Cx;
ParamCx = ParamCx./max_ACx*max_Cx;
Cx0 = ParamCx(1,1);
Cxa = ParamCx(2,1);
Cxa2 = ParamCx(3,1);
Cxq = ParamCx(4,1);
Cxde = ParamCx(5,1);
CxTc = ParamCx(6,1);
%Data Validation
ACx = ones(6,size_data_val);
ACx(2,:) = (Z_k1k1(11,data_validation)); %alpha
ACx(3,:) = (Z_k1k1(11,data_validation)).^2; %alpha^2
ACx(4,:) = (U_k(5,data_validation) - XX_k1k1(14,data_validation)).*chord./Z_k1k1(10,data_validation); %q
ACx(5,:) = de(data_validation); %de
ACx(6,:) = Tc1(data_validation); %throttle
%Residual Analysis
[mean_Cx, var_Cx, Rsq_Cx] = ResidualAnalysis(ParamCx, ACx, Cx(data_validation), ignore_data);
P_Cx = cov(ACx*ACx')*var_Cx;
figure;
for i = 1:6
     subplot(6,1,i)
     bar(P_Cx(i,:));
     title('Variance of \theta_{Cx}')
 end
%% Cz estimate coefficients
ACz = ones(5,size_data_eval);
ACz(2,:) = (Z_k1k1(11,data_evaluation)); %alpha
ACz(3,:) = (U_k(5,data_evaluation) - XX_k1k1(14,data_evaluation)).*chord./Z_k1k1(10,data_evaluation); %q
ACz(4,:) = de(data_evaluation); %de
ACz(5,:) = Tc1(data_evaluation); %throttle
max_ACz = max(abs(ACz),[],2);
in = find(max_ACz == 0);
max_ACz(in) = 1;
ACz = ACz./max_ACz;
max_Cz = max(abs(Cz));
ParamCz = (ACz*ACz')\ACz*Cz(data_evaluation)/max_Cz;
ParamCz = ParamCz./max_ACz*max_Cz;
Cz0 = ParamCz(1,1);
Cza = ParamCz(2,1);
Czq = ParamCz(3,1);
Czde = ParamCz(4,1);
CzTc = ParamCz(5,1);
%Data validation
ACz = ones(5,size_data_val);
ACz(2,:) = (Z_k1k1(11,data_validation)); %alpha
ACz(3,:) = (U_k(5,data_validation) - XX_k1k1(14,data_validation)).*chord./Z_k1k1(10,data_validation); %q
ACz(4,:) = de(data_validation); %de
ACz(5,:) = Tc1(data_validation); %throttle
%Residual Analysis
[mean_Cz, var_Cz, Rsq_Cz] = ResidualAnalysis(ParamCz, ACz, Cz(data_validation), ignore_data);
P_Cz = cov(ACz*ACz')*var_Cz;
%% Cm estimate coefficients
ACm = ones(5,size_data_eval);
ACm(2,:) = (Z_k1k1(11,data_evaluation)); %alpha
ACm(3,:) = (U_k(5,data_evaluation) - XX_k1k1(14,data_evaluation)).*chord./Z_k1k1(10,data_evaluation); %q
ACm(4,:) = de(data_evaluation); %de
ACm(5,:) = Tc1(data_evaluation); %throttle
max_ACm = max(abs(ACm),[],2);
in = find(max_ACm == 0);
max_ACm(in) = 1;
ACm = ACm./max_ACm;
max_Cm = max(abs(Cm));
ParamCm = (ACm*ACm')\ACm*Cm(data_evaluation)/max_Cm;
ParamCm = ParamCm./max_ACm*max_Cm;
Cm0 = ParamCm(1,1);
Cma = ParamCm(2,1);
Cmq = ParamCm(3,1);
Cmde = ParamCm(4,1);
CmTc = ParamCm(5,1);
%Data validation 
ACm = ones(5,size_data_val);
ACm(2,:) = (Z_k1k1(11,data_validation)); %alpha
ACm(3,:) = (U_k(5,data_validation) - XX_k1k1(14,data_validation)).*chord./Z_k1k1(10,data_validation); %q
ACm(4,:) = de(data_validation); %de
ACm(5,:) = Tc1(data_validation); %throttle
%Residual Analysis
[mean_Cm, var_Cm, Rsq_Cm] = ResidualAnalysis(ParamCm, ACm, Cm(data_validation), ignore_data);
P_Cm = cov(ACm*ACm')*var_Cm;
%% Cy estimate coefficients
ACy = ones(5,size_data_eval);
ACy(2,:) = (Z_k1k1(12,data_evaluation)); %beta
ACy(3,:) = (U_k(4,data_evaluation) - XX_k1k1(13,data_evaluation)).*b./Z_k1k1(10,data_evaluation); %p
ACy(4,:) = (U_k(6,data_evaluation) - XX_k1k1(15,data_evaluation)).*b./Z_k1k1(10,data_evaluation); %r
ACy(5,:) = da(data_evaluation); %da
ACy(6,:) = dr(data_evaluation); %dr
max_ACy = max(abs(ACy),[],2);
in = find(max_ACy == 0);
max_ACy(in) = 1;
ACy = ACy./max_ACy;
max_Cy = max(abs(Cy));
ParamCy = (ACy*ACy')\ACy*Cy(data_evaluation)/max_Cy;
ParamCy = ParamCy./max_ACy*max_Cy;
Cy0 = ParamCy(1,1);
Cyb = ParamCy(2,1);
Cyp = ParamCy(3,1)/2;
Cyr = ParamCy(4,1)/2;
Cyda = ParamCy(5,1);
Cydr = ParamCy(6,1);
%Data validation 
ACy = ones(5,size_data_val);
ACy(2,:) = (Z_k1k1(12,data_validation)); %beta
ACy(3,:) = (U_k(4,data_validation) - XX_k1k1(13,data_validation)).*b./Z_k1k1(10,data_validation); %p
ACy(4,:) = (U_k(6,data_validation) - XX_k1k1(15,data_validation)).*b./Z_k1k1(10,data_validation); %r
ACy(5,:) = da(data_validation); %da
ACy(6,:) = dr(data_validation); %dr
%Residual Analysis
[mean_Cy, var_Cy, Rsq_Cy] = ResidualAnalysis(ParamCy, ACy, Cy(data_validation), ignore_data);
P_Cy = cov(ACy*ACy')*var_Cy;
%% Cl estimate coefficients
ACl = ones(5,size_data_eval);
ACl(2,:) = (Z_k1k1(12,data_evaluation)); %beta
ACl(3,:) = (U_k(4,data_evaluation) - XX_k1k1(13,data_evaluation)).*b./Z_k1k1(10,data_evaluation); %p
ACl(4,:) = (U_k(6,data_evaluation) - XX_k1k1(15,data_evaluation)).*b./Z_k1k1(10,data_evaluation); %r
ACl(5,:) = da(data_evaluation); %da
ACl(6,:) = dr(data_evaluation); %dr
max_ACl = max(abs(ACl),[],2);
in = find(max_ACl == 0);
max_ACl(in) = 1;
ACl = ACl./max_ACl;
max_Cl = max(abs(Cl));
ParamCl = (ACl*ACl')\ACl*Cl(data_evaluation)/max_Cl;
ParamCl = ParamCl./max_ACl*max_Cl;
Cl0 = ParamCl(1,1);
Clb = ParamCl(2,1);
Clp = ParamCl(3,1)/2;
Clr = ParamCl(4,1)/2;
Clda = ParamCl(5,1);
Cldr = ParamCl(6,1);
%Data validation 
ACl = ones(5,size_data_val);
ACl(2,:) = (Z_k1k1(12,data_validation)); %beta
ACl(3,:) = (U_k(4,data_validation) - XX_k1k1(13,data_validation)).*b./Z_k1k1(10,data_validation); %p
ACl(4,:) = (U_k(6,data_validation) - XX_k1k1(15,data_validation)).*b./Z_k1k1(10,data_validation); %r
ACl(5,:) = da(data_validation); %da
ACl(6,:) = dr(data_validation); %dr
%Residual Analysis
[mean_Cl, var_Cl, Rsq_Cl] = ResidualAnalysis(ParamCl, ACl, Cl(data_validation), ignore_data);
P_Cl = cov(ACl*ACl')*var_Cl;
%% Cn estimate coefficients
ACn = ones(5,size_data_eval);
ACn(2,:) = (Z_k1k1(12,data_evaluation)); %beta
ACn(3,:) = (U_k(4,data_evaluation) - XX_k1k1(13,data_evaluation)).*b./Z_k1k1(10,data_evaluation); %p
ACn(4,:) = (U_k(6,data_evaluation) - XX_k1k1(15,data_evaluation)).*b./Z_k1k1(10,data_evaluation); %r
ACn(5,:) = da(data_evaluation); %da
ACn(6,:) = dr(data_evaluation); %dr
max_ACn = max(abs(ACn),[],2);
in = find(max_ACn == 0);
max_ACn(in) = 1;
ACn = ACn./max_ACn;
max_Cn = max(abs(Cn));
ParamCn = (ACn*ACn')\ACn*Cn(data_evaluation)/max_Cn;
ParamCn = ParamCn./max_ACn*max_Cn;
Cn0 = ParamCn(1,1);
Cnb = ParamCn(2,1);
Cnp = ParamCn(3,1)/2;
Cnr = ParamCn(4,1)/2;
Cnda = ParamCn(5,1);
Cndr = ParamCn(6,1);
%Data validation 
ACn = ones(5,size_data_val);
ACn(2,:) = (Z_k1k1(12,data_validation)); %beta
ACn(3,:) = (U_k(4,data_validation) - XX_k1k1(13,data_validation)).*b./Z_k1k1(10,data_validation); %p
ACn(4,:) = (U_k(6,data_validation) - XX_k1k1(15,data_validation)).*b./Z_k1k1(10,data_validation); %r
ACn(5,:) = da(data_validation); %da
ACn(6,:) = dr(data_validation); %dr
%Residual Analysis
[mean_Cn, var_Cn, Rsq_Cn] = ResidualAnalysis(ParamCn, ACn, Cn(data_validation), ignore_data);
P_Cn = cov(ACn*ACn')*var_Cn;
%% Plot
figure;
subplot(3,1,1)
plot(data_evaluation,da(data_evaluation),'.')
hold on
plot(data_validation,da(data_validation),'.')
subplot(3,1,2)
plot(data_evaluation,de(data_evaluation),'.')
hold on
plot(data_validation,de(data_validation),'.')
subplot(3,1,3)
plot(data_evaluation,dr(data_evaluation),'.')
hold on
plot(data_validation,dr(data_validation),'.')