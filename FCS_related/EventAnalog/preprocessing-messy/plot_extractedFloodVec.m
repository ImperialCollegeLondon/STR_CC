
bw = FMap > 0;
% imshow(bw)
% title('Image with Circles')


% stats = regionprops('table',bw,'PixelIdxList');
CC = bwconncomp(bw);

hold on;

% figure;
xlim([min(E_rps(:)),max(E_rps(:))]);
ylim([min(N_rps(:)),max(N_rps(:))]);

for i = 1:length(CC.PixelIdxList)
    ind = cell2mat(CC.PixelIdxList(i));
    
    plot(E_rps(ind),N_rps(ind),'o');
    hold on;
    text(E_rps(ind(1)),N_rps(ind(1))-2000,sprintf('%02d',i),'fontsize',8);
    drawnow;
end
%%
colorbar
title('extracted flooding Loc for one subWC')

sPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS\';
filename = [sPath,'Example_HydraulicNodeVector_subWC1'];
savePlot(filename,'XYWH',[150,0,600,600]);%,'needreply','N');



%% plot Info

pcolor(E,N,Topo);shading flat
cptcmap('DEM_print', 'mapping','direct');
caxis([70,300])





