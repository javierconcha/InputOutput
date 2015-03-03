%% ELM parameters

% ENVI ASCII Plot File [Wed Feb 25 11:13:20 2015]
% Column 1: Wavelength
% Column 2: Darkrad140929.txt:C2~~2##255,0,0
% Column 3: ONTOSRef_140919.txt:C2~~22##0,128,0
% Column 4: PIFrad140929.txt:C2~~4##0,0,255
% Column 5: RrsDTROC4272.txt:C2~~6##0,255,255
% Column 6: RrsDTROC4272.txt:C2~~7##255,0,255
ELMpar = [... 
 0.443000  47.715229  0.009902  77.846487  0.033414  0.033414;
 0.482600  38.305432  0.009372  75.344770  0.038789  0.038789;
 0.561300  21.273056  0.007824  65.115974  0.047262  0.047262;
 0.654600  10.047512  0.003416  57.065069  0.052816  0.052816;
 0.864600   2.932883  0.001856  40.646813  0.061532  0.061532;
 1.609000   0.161773  0.000688  10.359676  0.067213  0.067213;
 2.201000   0.029518  0.000066   2.766715  0.057581  0.057581];


Wavelength      = ELMpar(:,1);
Darkrad140929 	= ELMpar(:,2);
ONTOSRef_140919 = ELMpar(:,3);
PIFrad140929   	= ELMpar(:,4);
RrsDTROC4272 	= ELMpar(:,5);

figure
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ELMpar(:,1),Darkrad140929,'k','LineWidth',lw)
hold on 
plot(ELMpar(:,1),PIFrad140929 ,'LineWidth',lw)
legend('Dark: water ROI','Bright: PIF from L8')
title('Radiance values for ELM-based method','fontsize',fs)
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Radiance [W/m^2/sr/\mum]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5])
grid on

figure
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ELMpar(:,1),ONTOSRef_140919,'k','LineWidth',lw)
hold on 
plot(ELMpar(:,1),RrsDTROC4272,'LineWidth',lw)
legend('Dark: ONTOS field','Bright: PIF L8 refl. product')
title('Reflectance values for ELM-based method','fontsize',fs)
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Rrs [1/\pi]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5])
grid on

%% ROIs 
ROIstat = [...
  0.443000  0.009901  0.010249  0.010188  0.010691  0.010247  0.032758;
  0.482600  0.009343  0.009598  0.009578  0.010056  0.010751  0.040786;
  0.561300  0.007739  0.013285  0.013020  0.019612  0.015439  0.059070;
  0.654600  0.003437  0.008682  0.008326  0.012264  0.008745  0.072488;
  0.864600  0.001848  0.003361  0.003281  0.003740  0.003287  0.092126;
  1.609000  0.000699  0.000749  0.000777  0.000946  0.000744  0.115534;
  2.201000  0.000065  0.000056  0.000089  0.000254  0.000077  0.117337];

wl    = ROIstat(:,1);
ONTOS = ROIstat(:,2);
LONGS = ROIstat(:,3);
LONGN = ROIstat(:,4);
CRANB = ROIstat(:,5);
IBAYN = ROIstat(:,6);
Sand  = ROIstat(:,7);


figure
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ROIstat(:,1),ROIstat(:,2:end),'LineWidth',lw)
legend('ONTOS','LONGS','LONGN','CRANB','IBAYN','Sand')
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Rrs [1/\pi]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5]) 
grid on