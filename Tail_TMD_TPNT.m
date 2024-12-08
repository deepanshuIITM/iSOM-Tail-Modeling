% TAIL MODELLING FOR TUNED MASS DAMPER
clc
clear all
close all

%% Generate e6 samples
for k=1:1
n = 10^6;CI=0.5;
[yd2,Data] = mass_damp(n);

%% fit pdf
disp = sort(yd2); yd3 = disp(CI*length(disp)+1:n); yt=-yd3(1);
for i = 1:length(yd3)
    e_cdf1(i,1) =  i/(length(yd3)+1);
end

figure(k)
subplot(1,2,1)
plot(yd3,log(5-norminv(1-CI-(1-CI)*e_cdf1))); hold on;
subplot(1,2,2)
semilogy(yd3,1-CI-(1-CI)*e_cdf1); hold on;

%% 500 samples
rows = 500; 
[yd, Data_500] = mass_damp(rows);

%% SORTING and Plotting Emperical CDF
disp1 = sort(yd);
yd4 = disp1(CI*length(disp1)+1:length(yd));
for i = 1:length(yd4)
    e_cdf(i,1) =  i/(length(yd4)+1);
end

%% fit pdf
figure(k)
subplot(1,2,1)
ln_TPNT = log(5-norminv(1-CI-(1-CI)*e_cdf)); 
p_tp(:,k) = lsqfit_constr(yd4,ln_TPNT);
t = (linspace(-yt,1.1*yd3(end),20))';
pol = [ones(length(t),1) t t.^2 t.^3]*p_tp(:,k);
plot(yd4',ln_TPNT,'ks','markerSize',2.5); hold on;
plot(t', pol-pol(1)+ln_TPNT(1),'-ms','markerSize',2.5);

%% FIT CDF
b_cdf = normcdf(5-exp(pol-pol(1)+ln_TPNT(1)));
figure(k)
subplot(1,2,2)
semilogy(yd4',(1-CI-(1-CI)*e_cdf),'ks','markerSize',2.5); hold on;
semilogy(t',b_cdf,'-ms','markerSize',2.5);

%% SOM FOR ADDING MORE POINTS
init_sam = 300;
[yd_init, Data] = mass_damp(init_sam);
% Data = Data_500(1:init_sam,1:end); yd_init = Data_500(1:init_sam,1:end);
cnames = {'r1','r2','yd'};
sD  = som_data_struct(Data,'comp_names',cnames,'name','Data');
sTrain = som_normalize(sD,'var');
% sTrain  = som_data_struct(Data,'comp_names',cnames,'name','Data');
% sM_tail = som_make(Data,'algorithm','batch','msize',[37 35]);
sMap = modifiedsom_lininit(sTrain,'msize',[31 29]);
sM_tail = modifiedsom_batchtrain(sMap,sTrain);
figure(12)
som_show(som_denormalize(sM_tail))

Tail_data= Data(find(Data(:,end)>-yt),:);

%% Run SOM on Tail data
Xtrain = Tail_data;
sD1  = som_data_struct(Xtrain,'name','Tail','comp_names',cnames);
sTrain = som_normalize(sD1,'var');
% sM_tail = som_make(sTrain ,'algorithm','batch','msize',[33 31]);
sMap = modifiedsom_lininit(sTrain,'msize',[19 17]);
sM_tail = modifiedsom_batchtrain(sMap,sTrain,'trainlen',200);

figure(21)
som_show(som_denormalize(sM_tail))

%% CALCULATE TAIL
a=1;%0.38;
iData= sM_tail.codebook(find(sM_tail.codebook(:,end)>-a*yt),:);
Data_new = [iData; Data];
Xtrain_new = Data_new(:,1:end-1);
Data_new(find(Data_new(:,end)> -a*yt),end)=1;
Data_new(find(Data_new(:,end)<= -a*yt),end)=0;
label_new = int2str(Data_new(:,end));
sTrain_new  = som_data_struct(Xtrain_new,'name','Train','labels',label_new);

%% CALCULATE TAIL
j=0;
while(j< 500-init_sam)  
    mu = 1; sigma = 0.025; r1 = normrnd(mu,sigma,1,1);r2 = normrnd(mu,sigma,1,1);
    tData = [r1 r2];
    sTest  = som_data_struct(tData,'name','Test');
%     sM_tail1 = modified_som_autolabel_NL(sTest,sTrain_new,'vote');
    sM_tail1 = som_autolabel(sTest,sTrain_new,'vote');
if cell2mat(sM_tail1.labels)=='1'
    j=j+1;
    iData1(j,:) = tData;
end
end

%%
yd_fin = mass_damp(size(iData1,1),iData1);
y_new = unique([yd_fin(find(yd_fin> -yt));Tail_data(:,end)]);
accu(k) = (length(y_new)-size(Tail_data,1))/length(yd_fin);

%% ecdf
yi = sort(y_new);
for i = 1:length(yi)
    b_cdf(i) =  i/(length(yi)+1);
end
% Fitting 3rd degree monotonic polynomial
ln_TPNT_SOM = log(5-norminv(1-CI-(1-CI)*b_cdf)); 
p3(:,k) = lsqfit_constr(yi,ln_TPNT_SOM);

%% plot
pol_som = [ones(length(t),1) t t.^2 t.^3]*p3(:,k);
figure(k)
subplot(1,2,1)
plot(yi',ln_TPNT_SOM,'gs','markerSize',2.5);
plot(t', pol_som-pol_som(1)+ln_TPNT_SOM(1),'-rs','markerSize',2.5);
grid on; title('log-TPNT'); xlabel('y'); ylabel('log(5+ \beta)'); 

%% FIT CDF
s_cdf = normcdf(5-exp(pol_som-pol_som(1)+ln_TPNT_SOM(1)));
figure(k)
subplot(1,2,2)
semilogy(yi',(1-CI-(1-CI)*b_cdf),'gs','markerSize',2.5);hold on;
semilogy(t',s_cdf,'-rs','markerSize',2.5); 
grid on; title('P_f'); xlabel('y'); ylabel('P_f'); 

%% Computing error for 500 data
a = yd3; beta_act = log(5-norminv(1-CI-(1-CI)*e_cdf1)); 
beta_init =[ones(length(a),1) a a.^2 a.^3]*p_tp(:,k);
err_init_beta = beta_act./(beta_init-beta_init(1)+ln_TPNT(1));
figure(k)
subplot(1,2,1)
yyaxis right
plot(a,err_init_beta,'k-.','linewidth',0.1); hold on
cdf_init = normcdf(5-exp(beta_init-beta_init(1)+ln_TPNT(1)));
err_init_cdf = log(1-CI-(1-CI)*e_cdf1)./log(cdf_init);
figure(k)
subplot(1,2,2)
yyaxis right
semilogy(a,err_init_cdf,'k-.','linewidth',0.1);hold on

%% Computing error for SOM data 
beta_som =[ones(length(a),1) a a.^2 a.^3]*p3(:,k);
err_som_beta = beta_act./(beta_som-beta_som(1)+ln_TPNT_SOM(1));
cdf_som = normcdf(5-exp(beta_som-beta_som(1)+ln_TPNT_SOM(1)));
err_som_cdf = log(1-CI-(1-CI)*e_cdf1)./log(cdf_som);

figure(k)
subplot(1,2,1)
yyaxis right
plot(a,err_som_beta,'b-.','linewidth',0.1);
ylabel('$$\frac{Actual}{Predicted}$$','interpreter', 'latex','FontSize',10);
xline(1,'r-.',{'failure','Point'});
legend('1e07 Sample','Intial Data','log-TPNT_{500}','SOM Data','SOM log-TPNT','Err_{500}',...
    'Err_{SOM}','location','best')

figure(k)
subplot(1,2,2)
yyaxis right
semilogy(a,err_som_cdf,'b-.','linewidth',0.1);
ylabel('$$\frac{Actual}{Predicted}$$','interpreter', 'latex','FontSize',10);
xline(1,'r-.',{'failure','Point'});
legend('1e07 Sample','Intial Data','log-TPNT_{500}','SOM Data','SOM log-TPNT','Err_{500}',...
    'Err_{SOM}','location','best')

%% Compute Quantiles
pf = [0.1 0.01 0.001 0.0001]; D0 = 27; % D0 = 27;
for r = 1:length(pf)
q_init(k,r) =  a(find(round(cdf_init,4)==pf(r),1))+D0;
q_som(k,r) =  a(find(round(cdf_som,4)==pf(r),1))+D0;
q_orig(k,r) =  a(find(round((1-CI-(1-CI)*e_cdf1),4)==pf(r),1))+D0;
end
err_init(k,:) = 100*abs((q_init(k,:)-q_orig(k,:))./q_orig(k,:));
err_som(k,:) = 100*abs((q_som(k,:)-q_orig(k,:))./q_orig(k,:));
end
