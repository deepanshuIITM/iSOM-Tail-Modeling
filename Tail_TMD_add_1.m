% TAIL MODELLING FOR TMD 2DoF
clc
clear all
close all

%% Generate e6 samples
absorb = 2; damp1 = 1; damp2 =1; spr1=1; spr2=1;
for k=1:1
n = 10^6;CI=0.5; b_cdf1 = []; b_cdf2 = [];  b_cdf3 = []; b_cdf = [];
[yd2,data_e6] = TMD_3DOF(n,absorb,damp1,damp2,spr1,spr2);

%% fit pdf
disp = sort(yd2); yd3 = disp(CI*length(disp)+1:n); yt=-yd3(1);
for i = 1:length(yd3)
    e_cdf1(i,1) =  (length(yd3)+i)/(length(yd2)+1);
end

figure(k)
subplot(1,2,1)
plot(yd3,log(5+norminv(e_cdf1))); hold on;
subplot(1,2,2)
semilogy(yd3,1-e_cdf1); hold on;

%% 500 samples
rows = 500; 
[yd, Data_500] = TMD_3DOF(rows,absorb,damp1,damp2,spr1,spr2);

%% SORTING and Plotting Emperical CDF
disp1 = sort(yd);
yd4 = disp1(CI*length(disp1)+1:length(yd));
for i = 1:length(yd4)
    e_cdf(i,1) =  (length(yd4)+i)/(length(yd)+1);
end

%% fit pdf
ln_TPNT = log(5+norminv(e_cdf)); 
p_tp(:,k) = lsqfit_constr(yd4,ln_TPNT);
t = (linspace(-yt,yd3(end),15))';
pol = [ones(length(t),1) t t.^2 t.^3]*p_tp(:,k);
figure(k)
subplot(1,2,1)
plot(yd4',ln_TPNT,'ks','markerSize',2.5); hold on;
plot(t', pol-pol(1)+ln_TPNT(1),'-ms','markerSize',2.5);

%% FIT CDF
bcdf = normcdf(-5+exp(pol-pol(1)+ln_TPNT(1)));
figure(k)
subplot(1,2,2)
semilogy(yd4',(1-e_cdf),'ks','markerSize',2.5); hold on;
semilogy(t',1-bcdf,'-ms','markerSize',2.5);

%% SOM FOR ADDING MORE POINTS
init_sam = 400;
[yd_init, Data] = TMD_3DOF(init_sam,absorb,damp1,damp2,spr1,spr2);
% Data = Data_500(1:init_sam,1:end); yd_init = Data_500(1:init_sam,1:end);
cnames = {'r1','r2','r3','yd'};
sD  = som_data_struct(Data,'comp_names',cnames,'name','Data');
sTrain = som_normalize(sD,'var');
% sTrain  = som_data_struct(Data,'comp_names',cnames,'name','Data');
% sM_tail = som_make(Data,'algorithm','batch','msize',[37 35]);
sMap = modifiedsom_lininit(sTrain,'msize',[31 29]);
sM_tail = modifiedsom_batchtrain(sMap,sTrain);
figure(12)
som_show(som_denormalize(sM_tail))

Tail_data= Data(find(Data(:,end)>-yt),:);
tail_y = Tail_data(:,end);

%% Add data points in extreme of the tail
n1 = 10^5;
[~,Data1] = TMD_3DOF(n1,absorb,damp1,damp2,spr1,spr2);

Y1 = Data1(:,1:end-1);
Y1 = Y1(find(Y1(:,1)>1.03 & Y1(:,2)>1.065),:);
Y1 = Y1(find(Y1(:,1)<1.05),:);
y_end11 = TMD_3DOF(size(Y1,1),absorb,damp1,damp2,spr1,spr2,Y1);
y_end1 = y_end11(find(y_end11>5));

n2 = 10^4;
[~,Data1] = TMD_3DOF(n2,absorb,damp1,damp2,spr1,spr2);

Y2 = Data1(:,1:end-1);
Y2 = Y2(find(Y2(:,1)>1.03 & Y2(:,2)>1.045),:);
Y2 = Y2(find(Y2(:,1)<1.05),:);
y_end22 = TMD_3DOF(size(Y2,1),absorb,damp1,damp2,spr1,spr2,Y2);
y_end2 = y_end22(find(y_end22>-5));

%% compute cdf
yi1 = sort(y_end1);
for i = 1:length(yi1)
    b_cdf1(i) =  (n1-length(yi1)+i)/(n1+1);
end

yi2 = sort(y_end2);
for i = 1:length(yi2)
    b_cdf2(i) =  (n2-length(yi2)+i)/(n2+1);
end

yi3 = sort(tail_y);
for i = 1:length(yi3)
    b_cdf3(i) =   (init_sam-length(yi3)+i)/(init_sam+1);
end

b_cdf= [b_cdf3 b_cdf2 b_cdf1];
yi = unique(sort([yi3;yi2; yi1]));

%% Fitting 3rd degree monotonic polynomial
ln_TPNT_SOM = log(5+norminv(b_cdf')); 
p3(:,k) = lsqfit_constr(yi,ln_TPNT_SOM);

%% plot
pol_som = [ones(length(t),1) t t.^2 t.^3]*p3(:,k);
figure(k)
subplot(1,2,1)
plot(yi',ln_TPNT_SOM,'gs','markerSize',2.5);
plot(t', pol_som-pol_som(1)+ln_TPNT_SOM(1),'-rs','markerSize',2.5);
grid on; title('log-TPNT'); ylabel('log(5+ \beta)'); 
% xlim([-yt yd3(end)]);
%% FIT CDF
s_cdf = normcdf(-5+exp(pol_som-pol_som(1)+ln_TPNT_SOM(1)));
figure(k)
subplot(1,2,2)
semilogy(yi',1-b_cdf,'gs','markerSize',2.5);hold on;
semilogy(t',1-s_cdf,'-rs','markerSize',2.5); 
grid on; title('P_f'); ylabel('P_f'); 
% ylim([10^-6 1]);
%% Computing error for 500 data
a = yd3; beta_act = log(5+norminv(CI+(1-CI)*e_cdf1)); 
beta_init =[ones(length(a),1) a a.^2 a.^3]*p_tp(:,k);
err_init_beta = beta_act./(beta_init-beta_init(1)+ln_TPNT(1));

cdf_init = normcdf(-5+exp(beta_init-beta_init(1)+ln_TPNT(1)));
err_init_cdf = log(1-e_cdf1)./log(1-cdf_init);

figure(k)
subplot(1,2,1)
yyaxis right
plot(a,err_init_beta,'k-.','linewidth',0.1); hold on

figure(k)
subplot(1,2,2)
yyaxis right
semilogy(a,err_init_cdf,'k-.','linewidth',0.1);hold on

%% Computing error for SOM data 
beta_som =[ones(length(a),1) a a.^2 a.^3]*p3(:,k);
err_som_beta = beta_act./(beta_som-beta_som(1)+ln_TPNT_SOM(1));

cdf_som = normcdf(-5+exp(beta_som-beta_som(1)+ln_TPNT_SOM(1)));
err_som_cdf = log(1-e_cdf1)./log(1-cdf_som);

figure(k)
subplot(1,2,1)
yyaxis right
plot(a,err_som_beta,'b-.','linewidth',0.1);
ylabel('Ratio of log(5+ \beta)');xline(0,'r-.',{'failure','Point'});
legend('1e06 Sample','Intial Data','TPNT_{500}','SOM Data','SOM TPNT','Err_{500}',...
    'Err_{SOM}','location','best')

figure(k)
subplot(1,2,2)
yyaxis right
semilogy(a,err_som_cdf,'b-.','linewidth',0.1);
ylabel('Ratio of log(P_f)');xline(0,'r-.',{'failure','Point'});
legend('1e06 Sample','Intial Data','TPNT_{500}','SOM Data','SOM TPNT','Err_{500}',...
    'Err_{SOM}','location','best')
% end
%% Compute Quantiles
% pf = [0.1 0.01 0.001 0.0001]; D0 = 27;
% 
% for r = 1:length(pf)
% ini_root = roots([p_tp(4) p_tp(3) p_tp(2) (p_tp(1)-(log(5+norminv(1-pf(r)))+pol(1)-ln_TPNT(1)))]);
% q_init(k,r) =  real(ini_root(end)) +D0;
% som_root = roots([p3(4) p3(3) p3(2) (p3(1)-(log(5+norminv(1-pf(r)))+pol_som(1)-ln_TPNT_SOM(1)))]);
% q_som(k,r) = real(som_root(end)) +D0;
% q_orig(k,r) =  a(find(round((1-e_cdf1),5)==pf(r),1))+D0;
% end
% err_init(k,:) = 100*abs((q_init(k,:)-q_orig(k,:))./q_orig(k,:));
% err_som(k,:) = 100*abs((q_som(k,:)-q_orig(k,:))./q_orig(k,:));

end
%% BOXPLOT
% figure(15)
% for p=1:size(err_som,2)-2
% subplot(1,2,p)
% boxplot([err_init(:,p+2) err_som(:,p+2)],{'log-TPNT','log-TPNT-SOM'},'colors','bk');
% title(sprintf('Error in y for P_f =1e%d',-p-2))
% end
% 
% % Compactness
% figure(16)
% for p=1:size(p3,1)
% subplot(2,2,p)
% boxplot([p_tp(p,:)' p3(p,:)'],{'log-TPNT','log-TPNT-SOM'},'colors','bk');
% title(sprintf('a_%d',p-1))
% end
% sgtitle('Variation in log-TPNT fit','FontSize',11)
