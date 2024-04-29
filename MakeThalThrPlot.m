PlotThalGradientSlices(ThalConnMean,vox_coords_full,parula(256),"Mean seed connectivity",4);
print(['./figures/FigS1A.png'],'-dpng')
load('SeedsThr.mat', 'ThalConnMean')
figure
histogram(ThalConnMean,50)
ylabel('Number of seeds')
xlabel('Mean seed connectivity')
set(gca,'FontSize',14)
print(['./figures/FigS1B.png'],'-dpng')
ThalData2Plot = double(~(ThalConnMean>100));

c = PlotThalGradientSlices(ThalData2Plot,vox_coords_full,[lines(2)],"Coverage",4);

caxis([-.5 1.5])

c.Ticks = 0:1;
c.TickLabels = {'Included','Excluded'};
print(['./figures/FigS1C.png'],'-dpng')