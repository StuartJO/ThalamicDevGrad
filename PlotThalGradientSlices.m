function PlotThalGradientSlices(GRAD_DATA,GRAD_DATA_vox_coord,grad_cmap,TITLE,thr,climits)

if nargin < 6
   climits = [min(GRAD_DATA) max(GRAD_DATA)];
end

load('Week40_brain_thal_data.mat','braindata','thalmask')

inmask = find(thalmask==1);
[X,Y,Z] = ind2sub(size(thalmask),inmask);
ThalMaskDists = pdist2([X Y Z],GRAD_DATA_vox_coord);
[ThalMask_nearest_seedmm,ThalMask_nearest_seed] = min(ThalMaskDists');

ThalVoxClust = changem(ThalMask_nearest_seed,GRAD_DATA,1:length(GRAD_DATA_vox_coord));

thalmaskclust = thalmask;
thalmaskclust(thalmaskclust==0) = NaN;
thalmaskclust(inmask) = ThalVoxClust;

thalmaskclust(inmask(ThalMask_nearest_seedmm>thr)) = NaN;

figure

slices_cmap = [254,178,76;...
253,141,60;...
252,78,42;...
227,26,28;...
189,0,38;...
128,0,38]./255;
axes_slices = axes('Position',[.6 .57 .4 .43]);
brain_back_slice = (squeeze(braindata(103,:,:))');
h_brain_back = imagesc(brain_back_slice);
set(gca,'YDir','normal')
colormap(axes_slices,gray(256))
slices = [116   110   104    98    92    86];
for i = 1:6
    hold on
plot([1 size(brain_back_slice,2)],[slices(i) slices(i)],'Color',slices_cmap(i,:))
end
%axis image
axes_slices.XTick = [];
axes_slices.YTick = [];

yslice_range = 85:145;
xslice_range = 70:150;
for i = 1:length(slices)
slice = slices(i);
brain_slices{i} = flipud(squeeze(braindata(xslice_range,yslice_range,slice))');
thal_slices{i} = flipud(squeeze(thalmaskclust(xslice_range,yslice_range,slice))');
end

brain_back_slice = [brain_slices{1} brain_slices{2} brain_slices{3}; brain_slices{4} brain_slices{5} brain_slices{6}];
thal_slice = [thal_slices{1} thal_slices{2} thal_slices{3}; thal_slices{4} thal_slices{5} thal_slices{6}];

axes_back = axes('Position',[0 0 1 .57]);
h_brain_back = imagesc(brain_back_slice);
colormap(axes_back,gray(256));
%axis image

P = get(axes_back,'pos');
axes_front = axes('Position',P);
h = imagesc(thal_slice);
set(h, 'AlphaData', ~isnan(thal_slice))
set(axes_front,'color','none')
%axis image
axes_front.XTick = [];
axes_front.YTick = [];
axes_back.XTick = [];
axes_back.YTick = [];
if climits(1)~=climits(2)
caxis(climits);
end
colormap(axes_front,grad_cmap);
% axes_back.Position = [0 -.2 1 1];
% axes_front.Position = [0 -.2 1 1];

XtextSpacing = size(brain_back_slice,2)/3;
YtextSpacing = size(brain_back_slice,1)/2;

x = [1 XtextSpacing (XtextSpacing*2) 1 XtextSpacing (XtextSpacing*2)]+10;
y = [1 1 1 YtextSpacing YtextSpacing YtextSpacing]+20;

for i = 1:6
text(x(i),y(i),['\bf',num2str(i)],'Color',slices_cmap(i,:),'FontSize',28)
end

hold on
plot([0 size(brain_back_slice,2)+.5],[YtextSpacing YtextSpacing],'k','LineWidth',2) 

plot([XtextSpacing XtextSpacing],[0 size(brain_back_slice,1)+.5],'k','LineWidth',2) 
plot([XtextSpacing XtextSpacing]*2,[0 size(brain_back_slice,1)+.5],'k','LineWidth',2) 

c = colorbar('Location','northoutside','FontSize',18);

c.Position = [0.025,0.638095238095239,.5-0.05,0.066666666666666];
c.Position = [0.05,0.638095238095239,0.5,0.066666666666666];

annotation(gcf,'textbox',...
    [0.05,0.754761904761905,0.5,0.238095238095242],...
    'VerticalAlignment','middle',...
    'String',TITLE,...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
