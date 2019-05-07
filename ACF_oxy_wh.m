%%%%%%% residence time of water oxygen around sulfur
clc;close all;clear all;

sow_rad=4.1; %1st coordination shell

chains=20;side_chains=10;
sow_rad=4.1; %1st coordination shell


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


%assign sulfur types as per DP
if DP==0
    sulfur_old_type=4;
    sulfur_new_type=4;
elseif DP==3
    sulfur_old_type=4;
    sulfur_new_type=8;
elseif DP==5
    sulfur_old_type=4;
    sulfur_new_type=6;
elseif DP==7
    sulfur_old_type=4;
    sulfur_new_type=8;
elseif DP==10
    sulfur_old_type=4;
    sulfur_new_type=4;
end

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



Max_ts=2000; %number of timesteps
max_lags=1800;
M=2000; %number of time origins

if DP==10
    path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda))
else
    path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda),'/restart')
end
cd(path);

sulfur_number=side_chains*chains;


%%preallocate arrays
for i=1:sulfur_number
    Sulfur(i).array=zeros(Max_ts,length(oxy_w));
 
    SulfurCF(i).devCF=zeros(max_lags+1,length(oxy_w));
    SulfurCF(i).sumdevCF=zeros(max_lags+1,1);
    
    
    SulfurCF(i).ACF=zeros(max_lags+1,length(oxy_w));
    SulfurCF(i).sumACF=zeros(max_lags+1,1);
    
end



try
    dump = fopen('wrap.lammpstrj','r');
catch
    error('Dumpfile not found!');
end

ts_count=0;
while ~feof(dump)
    if ts_count>=Max_ts
        break;
    end
    
    id = fgetl(dump);
    if (strcmpi(id,'ITEM: TIMESTEP'))
        
        x=str2num(fgetl(dump));       
        id=fgetl(dump);
        %timestep(i) = str2num(fgetl(dump));
        if (strcmpi(id,'ITEM: NUMBER OF ATOMS'))
            Natoms = str2num(fgetl(dump));
            id=fgetl(dump);
        end
        if (strcmpi(id,'ITEM: BOX BOUNDS pp pp pp'))
            x_bound = str2num(fgetl(dump));
            y_bound = str2num(fgetl(dump));
            z_bound = str2num(fgetl(dump));
            id=fgetl(dump);
        end
        if (strcmpi(id(1:11),'ITEM: ATOMS'))
            ts_count=ts_count+1
            counter=0;
            atom_data=zeros(Natoms,3); %preallocate atom_data
            for j = 1 : 1: Natoms
                atom_data_temp=zeros(1,6);
                atom_data_temp= str2num(fgetl(dump));
               
                if ts_count==1
                    if (atom_data_temp(2)==sulfur_old_type) || (atom_data_temp(2)==sulfur_new_type)
                        counter=counter+1;
                        sulfur_index(counter)=j; %sulfur array contains index of sulfur atoms populated in the first timestep
                    end
                end
                atom_data_temp([1 2 6])=[];
                atom_data(j,:)=atom_data_temp;
               
            end
            
        end
        

        
       
        for sul_count=1:length(sulfur_index)
            
            %calculate boolean distance for oxygen atoms
            for oxy_count=1:length(oxy_w)
                
                sq_dist= sum((atom_data(sulfur_index(sul_count),1:3)-atom_data(oxy_w(oxy_count),1:3)).^2);
             
                
                if sq_dist> sow_rad^2
                    dist=0;
                else
                    dist=1;
                end
                
                Sulfur(sul_count).array(ts_count,oxy_count)=dist;
            end
          
                        
        end
      
    end
      clear atom_data;
end



path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda),'/Residence_time')
cd(path);
save Sulfur_data -v7.3 Sulfur
% 
%devanthan cf calculation
%water CF
for sul_count=1:sulfur_number
   
    for oxy_count=1:length(oxy_w)
      
        temp= xcorr(Sulfur(sul_count).array(:,oxy_count),max_lags,'none'); %biased and non-normalized auto-correlation
        temp=temp(max_lags+1:2*max_lags+1);
         SulfurCF(sul_count).devCF(:,oxy_count)= temp;    
    end
    SulfurCF(sul_count).sumdevCF=sum(SulfurCF(sul_count).devCF,2);
end
clear temp;


save SulfurCF.mat -v7.3 SulfurCF
% 
% % %normal ACF computation
% 
%water
for sul_count=1:sulfur_number
   
    for oxy_count=1:length(oxy_w)
        temp1=(sum(Sulfur(sul_count).array(:,oxy_count).^2))/length(Sulfur(sul_count).array(:,oxy_count)); %normalizing factor
        temp= xcorr((Sulfur(sul_count).array(:,oxy_count)),max_lags,'unbiased')*(1/temp1); %autocorrelation unbiased and normalized
        temp=temp(max_lags+1:2*max_lags+1);
         SulfurCF(sul_count).ACF(:,oxy_count)=temp;      
    end
    SulfurCF(sul_count).sumACF=sum(SulfurCF(sul_count).ACF,2)/length(oxy_w);
end
% hydronium


clear temp temp1
save SulfurCF.mat -v7.3 SulfurCF







            
            
