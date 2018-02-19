function [alpha_out]=get_Ge_k(lambda,pt)
%Input the desired wavelength you want to look for in microns, will do a
%simple two-point interpolation from the provided data to find the
%effective thermal diffusivity at that point. pt is a plot command for the interpolated value on the data.

if nargin<2
    pt=0;
end

Ge_heat_capacity=321.4; %in J/kg*K
Ge_density=5.323e3; %in kg/m^3

Ge_k=dlmread('ge_300K_1dTTG.txt');

if lambda~=0
    
    lambda_units=lambda*10^(-6);
    
    [~,ind]=min(abs(Ge_k(:,1)-lambda_units));
    
    if abs(Ge_k(ind+1,1)-lambda_units)<abs(Ge_k(ind-1,1)-lambda_units)
        ind2=ind+1;
    else
        ind2=ind-1;
    end
    
    m=(Ge_k(ind2,2)-Ge_k(ind,2))/(Ge_k(ind2,1)-Ge_k(ind,1));
    b=Ge_k(ind,2)-m*Ge_k(ind,1);
    
    k_out=m*lambda_units+b;
    
    alpha_out=k_out/(Ge_density*Ge_heat_capacity);
    alpha_out=alpha_out*10^4; %to get in cm^2/s
    
else
    alpha_out=0;
end

if pt
    figure
    plot(Ge_k(:,1)/(10^(-6)),(Ge_k(:,2)/(Ge_density*Ge_heat_capacity))*10^4,'-','Color',[0 0 0.75],'LineWidth',1.25)
    xlim([0.5 15])
    if lambda~=0
        hold on
        plot(lambda,alpha_out,'ro','MarkerSize',8,'MarkerFaceColor','r');
    end
    xlabel(['Grating spacing [' 956 'm]'])
    ylabel('Effective Thermal Diffusivity [cm^2/s]')
end
end