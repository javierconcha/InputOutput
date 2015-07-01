function [M_RGB,Rrs] = createtile(input,CHlimit,SMlimit,CDlimit,...
    DPFlimit,c,LUTconcDPF,LUTconc,LUT,tilex,tiley)

rule = strcmp(c{1}(:),input)& ...
    LUTconc(:,1)==CHlimit & LUTconc(:,2)==SMlimit & LUTconc(:,3)==CDlimit & ...
    LUTconcDPF(:)==DPFlimit;

Rrs = LUT(rule,:);

M_RGB = nan(tilex,tiley,3);

M_RGB(:,:,1) = Rrs(4).*ones(tilex,tiley); % R
M_RGB(:,:,2) = Rrs(3).*ones(tilex,tiley); % G
M_RGB(:,:,3) = Rrs(2).*ones(tilex,tiley); % B