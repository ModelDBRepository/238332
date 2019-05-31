

colors=[[255/255,210/255,50/255]; [255/255,0/255,0/255]; [70/255,70/255,180/255]; [0/255,170/255,220/255]; [150/255,80/255,30/255]; [0/255,150/255,60/255]; [255/255,170/255,140/255]];

load P1000

Hin=zeros(7,1);
Hout=zeros(7,1);
types=[1,7,2,3,4,5,6];

for num=1:7
    x=sum(P);
    x=x(find_indexes(types(num)));
    Hin(num)=mea_degree_heterogeneity_hu_wang(x');
    
    x=sum(P'); 
    x=x(find_indexes(types(num)));
    Hout(num)=mea_degree_heterogeneity_hu_wang(x');
    
    subplot(1,2,1)
    hold on
    h=bar(num,Hin(num,1));
    set(h,'FaceColor',colors(types(num),:),'FaceAlpha',1,'BarWidth',1)
    subplot(1,2,2)
    hold on
    h=bar(num,Hout(num,1));
    set(h,'FaceColor',colors(types(num),:),'FaceAlpha',1,'BarWidth',1)
end
set(gcf, 'Position', [100, 100, 1000, 500])
subplot(1,2,1)
ylim([0,0.3])
ylabel('Heterogeneity index')
xlabel('In-degree')

subplot(1,2,2)
ylim([0,0.3])
xlabel('Out-degree')










