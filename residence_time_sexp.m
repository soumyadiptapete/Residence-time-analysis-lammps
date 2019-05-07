%extract residence time

clc;close all;clear all;

chains=20;side_chains=10;


%getting DP and lambda from path address automatically
path=cd;

DP_str='ab';lambda_str='cd';
k=strfind(path,'Nafion_');
x=length('Nafion_');
DP_str(1)=path(k(2)+x);
if path(k(2)+x+1)~='/'
    DP_str(2)=path(k(2)+x+1);
else
    DP_str(2)=[];
end
DP=str2num(DP_str)

k=strfind(path,'lambda_');
x=length('lambda_');
lambda_str(1)=path(k+x);
if path(k+x+1)~='/'
    lambda_str(2)=path(k+x+1);
else
    lambda_str(2)=[];
end
lambda=str2num(lambda_str)


atoms_per_chain=692-DP;

sulfur_atoms=(10-DP)*chains;


hydronium_mols=DP*chains;
water_mols=lambda*chains*10-hydronium_mols;

water_start=atoms_per_chain*chains+1;
water_end=water_start+water_mols*3-1;


hydrogen_atoms=hydronium_mols*3;
oxygen_atoms=hydronium_mols;


water=water_start:water_end; %water atoms index
oxy_w=water_start+1:3:water_end-1; %oxygen of water atoms index



path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda),'/Residence_time')
cd(path);

load SulfurCF.mat


sulfur_number=200;

%plotting

for sul_count=1:sulfur_number
        
          y=SulfurCF(sul_count).sumdevCF;
          sumdevCF(:,sul_count)=y;
           if sum(imag(y))>0
            disp(sul_count);
           end

          if sum(isnan(y))==0 && sum(isinf(y))==0
              x=(1:length(SulfurCF(sul_count).sumdevCF))-1;
              y=y/y(1);
              
              hold on;
              semilogx(x,y,'o--')
              hold on;
              clear x y f
          else
              m=isnan(y)+isinf(y);
              sul_count
              x=(1:length(SulfurCF(sul_count).sumdevCF))-1;
              if length(find(m))<length(x)
              x(find(m))=[];
              y(find(m))=[];
               y=y/y(1);
              hold on;
              semilogx(x,y,'o--');
              hold on;
              clear x y f
                
              end
          end
end
 savefig('tau_w_dcf_norm');
 save sumdevCF.mat -v7.3 sumdevCF
clear sumdevCF;
 hold off;

disp('water_fitting');
%fitting
for sul_count=1:sulfur_number
    sul_count
    y=SulfurCF(sul_count).sumdevCF;
    x=(1:length(SulfurCF(sul_count).sumdevCF))-1;
    
    g=fittype('a*exp(-(x/b)^c)');
    
    %initial guess 
    
    b_ig=2;
    a_ig=y(1);c_ig=0.7;
    [f,gof]=fit(x',y,g,'Robust','Bisquare','StartPoint',[a_ig b_ig c_ig],'MaxIter',2000,'MaxFunEvals',2000);
    rsqr_w(sul_count)=gof.adjrsquare;
    tau_w_dcf(sul_count)=f.b;
    beta_w_dcf(sul_count)=f.c; 
    clear x y f
    
end


 save res_time_dcf_sexp.mat -v7.3 tau_w_dcf rsqr_w beta_w_dcf





