% TAIL MODELLING FOR TMD-1
clc
clear all
close all

%% Generate e6 samples
for k=1:1
n = 10^6;CI=0.5; b_cdf1 = []; b_cdf2 = [];  b_cdf3 = []; b_cdf = [];
[yd2,data2] = mass_damp(n);

%% fit pdf
disp = sort(yd2); yd3 = disp(CI*length(disp)+1:n); yt=-yd3(1);
for i = 1:length(yd3)
    e_cdf1(i,1) =  (length(yd3)+i)/(length(yd2)+1);
end

% figure(k)
% subplot(1,2,1)
% plot(yd3,log(5+norminv(e_cdf1)),'linewidth',1.5); hold on;
% subplot(1,2,2)
% semilogy(yd3,1-e_cdf1,'linewidth',1.5); hold on;

%% 500 samples
xx = 0.9 + 0.1*lhsdesign(500,2);
rows = 500; 
[yd, Data_500] = mass_damp(rows,xx);
% [yd, Data_500] = mass_damp(rows);

%% SORTING and Plotting Emperical CDF
disp1 = sort(yd);
yd4 = disp1(CI*length(disp1)+1:length(yd));
for i = 1:length(yd4)
    e_cdf(i,1) =  (length(yd4)+i)/(length(yd)+1);
end

%% fit pdf
ln_TPNT = log(5+norminv(e_cdf)); 
p_tp(:,k) = lsqfit_constr(yd4,ln_TPNT);
t = (linspace(-yt,1.15*yd3(end),15))';
pol = [ones(length(t),1) t t.^2 t.^3]*p_tp(:,k);
% figure(k)
% subplot(1,2,1)
% plot(yd4',ln_TPNT,'ks','markerSize',2.5); hold on;
% plot(t', pol,'-k','linewidth',1);

%% FIT CDF
bcdf = normcdf(-5+exp(pol));
% figure(k)
% subplot(1,2,2)
% semilogy(yd4',(1-e_cdf),'ks','markerSize',2.5); hold on;
% semilogy(t',1-bcdf,'-k','linewidth',1);

%% SOM FOR ADDING MORE POINTS
init_sam = 400;
xx = 0.95 + 0.1*lhsdesign(init_sam,2); %LHS Sampling
[yd_init, Data] = mass_damp(init_sam,xx);
% Data = Data_500(1:init_sam,1:end); yd_init = Data_500(1:init_sam,1:end);
cnames = {'$r_1$','$r_2$','$y_d$'};
sD  = som_data_struct(Data,'comp_names',cnames,'name','Data');
sTrain = som_normalize(sD,'range');
% sTrain  = som_data_struct(Data,'comp_names',cnames,'name','Data');
% sM_tail = som_make(Data,'algorithm','batch','msize',[37 35]);
sMap = modifiedsom_lininit1(sTrain,'msize',[17 17]);%
sM_tail = modifiedsom_batchtrain(sMap,sTrain);
som_plot = som_denormalize(sM_tail, sD);
figure(12)
som_show(som_plot,'umat', [3], 'comp', 'all')

%% RoI 
h_RoI=zeros(som_plot.topol.msize(1)*som_plot.topol.msize(2),1); 
h_RoI(find(som_plot.codebook(:,end)>1.0))=1;

% som_show_add('hit',h_RoI,'Markersize',0.4,'MarkerColor','g','EdgeColor','none','Subplot',2:4);
% som_show_add('hit',h_RoI,'Markersize',1,'MarkerColor','none','EdgeColor','r','Subplot',2:4);
% som_show_add('hit',h_RoI,'Markersize',0.2,'MarkerColor','b','EdgeColor','none','Subplot',2:4);

%%
Tail_data= Data(find(Data(:,end)>-yt),:);
tail_y = Tail_data(:,end);
length(tail_y)

%% Add data points in extreme of the tail
n1 = 5*10^3; % 10^5 earlier
xx = 0.95 + 0.1*lhsdesign(n1,2);
% [~,Data1] = mass_damp(n1);
[~,Data1] = mass_damp(n1,xx);

Y1 = Data1(:,1:end-1);
Y1 = Y1(find(Y1(:,1)>1.04 & Y1(:,2)>1.04),:);  %% Y1 = Y1(find(Y1(:,1)>1.025 & Y1(:,2)>1.025),:)
y_end11 = mass_damp(size(Y1,1),Y1); 
Y1 = Y1(find(y_end11>1.3),:);
y_end1 = y_end11(find(y_end11>1.3)); %20,1.9 earlier

n2 = 5*10^3;
xx = 0.95 + 0.1*lhsdesign(n2,2);
% [~,Data1] = mass_damp(n2);
[~,Data1] = mass_damp(n2,xx);

Y2 = Data1(:,1:end-1);
Y2 = Y2(find(Y2(:,1)<0.96 & Y2(:,2)<0.96),:);
y_end22 = mass_damp(size(Y2,1),Y2); 
Y2 = Y2(find(y_end22>1),:);
y_end2 = y_end22(find(y_end22>1)); %0,1.2

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
yi = [yi3; yi2; yi1];

%% Fitting 3rd degree monotonic polynomial
ln_TPNT_SOM = log(5+norminv(b_cdf')); 
p3(:,k) = lsqfit_constr(yi,ln_TPNT_SOM);

%% plot
pol_som = [ones(length(t),1) t t.^2 t.^3]*p3(:,k);
% figure(k)
% subplot(1,2,1)
% plot(yi(1:size(tail_y))',ln_TPNT_SOM(1:size(tail_y)),'gs','markerSize',2.5);
% plot(yi(1+size(tail_y,1):size(yi,1))',ln_TPNT_SOM(1+size(tail_y,1):size(yi,1)),...
%     'bs','markerSize',2.5);
% plot(t', pol_som,'-r','linewidth',1);
% grid on; 
% % title('log-TPNT'); ylabel('log(5+ \beta)'); 
% ylim([1.6 2.3]);

%% FIT CDF
s_cdf = normcdf(-5+exp(pol_som)); som_cdf = 1-b_cdf;
% figure(k)
% subplot(1,2,2)
% semilogy(yi(1:size(tail_y))',som_cdf(1:size(tail_y)),'gs','markerSize',2.5);hold on;
% semilogy(yi(1+size(tail_y,1):size(yi,1))',som_cdf(1+size(tail_y,1):size(yi,1)),...
%     'bs','markerSize',2.5);
% semilogy(t',1-s_cdf,'-r','linewidth',1); 
% grid on;
% % title('P_f'); ylabel('P_f'); 
% ylim([10^-5 1]);

%% Computing error for 500 data
a = yd3; beta_act = log(5+norminv(e_cdf1)); 
beta_init =[ones(length(a),1) a a.^2 a.^3]*p_tp(:,k);
err_init_beta = beta_act./(beta_init);

cdf_init = normcdf(-5+exp(beta_init));
err_init_cdf = log(1-e_cdf1)./log(1-cdf_init);

% figure(k)
% subplot(1,2,1)
% %yyaxis right
% %plot(a,err_init_beta,'k-.','linewidth',0.1); hold on
% 
% figure(k)
% subplot(1,2,2)
% %yyaxis right
% %semilogy(a,err_init_cdf,'k-.','linewidth',0.1);hold on

%% Computing error for SOM data 
beta_som =[ones(length(a),1) a a.^2 a.^3]*p3(:,k);
err_som_beta = beta_act./(beta_som);

cdf_som = normcdf(-5+exp(beta_som));
err_som_cdf = log(1-e_cdf1)./log(1-cdf_som);

% figure(k)
% subplot(1,2,1)
% %yyaxis right
% %plot(a,err_som_beta,'b-.','linewidth',0.1);
% grid on; 
% % title('log-TPNT', 'interpreter','latex'); 
% xlabel('$y_d$', 'interpreter','latex','FontName','times','FontSize',14); 
% ylabel('log(5+ $\beta$)', 'interpreter','latex','FontName','times','FontSize',14); 
% xline(1,'r-.',{'Failure','Point'}, 'interpreter','latex','FontName','times','FontSize',12);
% legend('1e05 Sample','Intial Sample','log-TPNT$_{500}$','iSOM$_{400}$','Adaptive Sample','log-TPNT-iSOM','location','best', 'interpreter','latex',...
%     'FontName','times','FontSize',11)
% 
% figure(k)
% subplot(1,2,2)
% %yyaxis right
% %semilogy(a,err_som_cdf,'b-.','linewidth',0.1);
% grid on; 
% % title('$P_f$', 'interpreter','latex'); 
% xlabel('$y_d$', 'interpreter','latex','FontName','times','FontSize',14); 
% ylabel('$P_f$', 'interpreter','latex','FontName','times','FontSize',14); 
% xline(1,'r-.',{'Failure','Point'}, 'interpreter','latex','FontName','times','FontSize',12);
% legend('1e05 Sample','Intial Sample','log-TPNT$_{500}$','iSOM$_{400}$','Adaptive Sample','log-TPNT-iSOM','location','best', 'interpreter','latex',...
%     'FontName','times','FontSize',11)

%%
R=0.01; J=0.01;y0=27;
mu = 1; sigma = 0.025; 
[r1, r2] = meshgrid(linspace(0.9,1.1,100),linspace(0.9,1.1,100));

for i = 1:100
    for j = 1:100
        yd_plot(i,j) = abs(1-(1/r2(i,j))^2)/(sqrt((1-R*(1/r1(i,j))^2 -(1/r1(i,j))^2-(1/r2(i,j))^2 +...
            (1/r1(i,j))^2*(1/r2(i,j))^2)^2 +4*(J^2)*((1/r1(i,j)) -(1/r1(i,j))*(1/r2(i,j)^2))^2)*y0); 
    end
end

[~, Data]=mass_damp(500);
x = Data(:,1); y = Data(:,2); z = Data(:,3);
% figure(101)
% surf(r1,r2,yd_plot,'FaceAlpha',0.5, 'EdgeColor', 'none')
% hold on;
% scatter3(x,y,z,10,'ro','filled');
% scatter3([Y1(:,1); Y2(:,1)],[Y1(:,2); Y2(:,2)],[y_end1; y_end2],10,'bo','filled');
% legend('Actual LSF','Initial Sample','Adaptive samples','fontname','times new roman');
% xlabel('$r_1$', 'interpreter','latex'); ylabel('$r_2$', 'interpreter','latex'); zlabel('$y_d$', 'interpreter','latex');

%% Compute Quantiles
pf = [0.01 0.001 0.00013 0.00001]; D0 = 27;
for r = 1:length(pf)
ini_root = roots([p_tp(4,k) p_tp(3,k) p_tp(2,k) (p_tp(1,k)-(log(5+norminv(1-pf(r)))))]);
q_init(k,r) =  real(ini_root(1))*D0;
som_root = roots([p3(4,k) p3(3,k) p3(2,k) (p3(1,k)-(log(5+norminv(1-pf(r)))))]);
q_som(k,r) = real(som_root(1))*D0;
q_orig(k,r) =  a(find(round((1-e_cdf1),5)==pf(r),1))*D0;
end
err_init(k,:) = 100*abs((q_init(k,:)-q_orig(k,:))./q_orig(k,:));
err_som(k,:) = 100*abs((q_som(k,:)-q_orig(k,:))./q_orig(k,:));

end

%% BOXPLOT
figure(101)
for p=1:size(err_som,2)-2
subplot(1,2,p)
boxplot([err_init(:,p+2) err_som(:,p+2)],{'log-TPNT','log-TPNT-SOM'},'colors','bk');
title(sprintf('Error in y for P_f =1e%d',-p-3))
end

% Compactness
figure(102)
for p=1:size(p3,1)
subplot(2,2,p)
boxplot([p_tp(p,:)' p3(p,:)'],{'log-TPNT','log-TPNT-SOM'},'colors','bk');
title(sprintf('a_%d',p-1))
end
% sgtitle('Variation in log-TPNT fit','FontSize',11)

%% Calculation of posterior probability Using Bayes Rule
% mu = mean(Data2(:,1:end-1)); cor_mat = cov(Data2(:,1:end-1));
% x = Y1(:,1:end-1);
% for m=1:size(x,1)
% lik_p(m) = (1/(2*pi*det(cor_mat))^(length(mu)/2))*exp(-(x(m,:)-mu)*inv(cor_mat)*(x(m,:)-mu)'/2);
% end
% figure(22)
% plot(yd3,log(5+norminv(e_cdf1)),'linewidth',1.5); hold on;
% plot(t', pol_som,'-r','linewidth',2); hold on;
% t1 = (linspace(ln_TPNT(1),1.1*ln_TPNT(end),20))';
% ln_TPNT_SOM1 = log(5+norminv(b_cdf')); 
% p3a(:,k) = lsqfit_constr(ln_TPNT_SOM1,yi);
% pol_som1 = [ones(length(t1),1) t1 t1.^2 t1.^3]*p3a(:,k);
% plot(pol_som1,t1','-k','linewidth',2);
% curvature(ln_TPNT_SOM,yi);
% legend('Actual','\beta = f(yd)','y(d) = f(\beta)','SOM data')
% ylim([1.3 2.4]);