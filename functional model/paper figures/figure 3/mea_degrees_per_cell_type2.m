
%% generate real network (tadpole) from probability matrix P
n=1382; % number of cells
types = {'RB','dlc','aIN','cIN','dIN','MN','dla'};
num_types = {1,2,3,4,5,6,7};
type_number = [0,126,104,136,384,236,338,58];
correct = [1,7,2,3,4,5,6];

v=zeros(1,8);
for i=1:8
    v(i)=sum(type_number(1:i));
end
colors=[[1 1 0]; [1 0 0]; [0 0 0.55]; [0 1 1]; [0.6 0.3 0.1]; [0 .5 0]; [1,0.4,0.6]];

names = zeros(1,n);
for i = 1:length(v)-1
    ind=v(i)+1:v(i+1);
    for j=ind
        names(j) = num_types{i};
    end
end

load P1000
[n,m]=size(P);
neurons=1:n;
load pos.txt
pos=pos';

path='/home/andrea/Desktop/tadpole-spinal-cord-models/probabilistic model/anatomical adjacency matrixes/';
m=1000;
indeg_all=zeros(m,n);
outdeg_all=zeros(m,n);
for i=1:m
    i
    load(strcat(path,'A_connectome',num2str(i),'.mat'));
    indeg_all(i,:)=sum(A);
    outdeg_all(i,:)=sum(A');
end
indeg=mean(indeg_all);
outdeg=mean(outdeg_all);
std_indeg=std(indeg_all);
std_outdeg=std(outdeg_all);

%% inward and outward degrees distribution and correlation
count=0;
for num=1:7
    figure(1)
    subplot(2,1,1)
    hold on
    idx=(v(correct(num))+1):v(correct(num))+(v(correct(num)+1)-v(correct(num)))/2;
    c1 = indeg(idx)-std_indeg(idx);
    c2 = indeg(idx)+std_indeg(idx);
    x = pos(idx)+1500*(num-1)-500;
    
    x2 = [x, fliplr(x)];
    inBetween = [c1, fliplr(c2)];
    h=fill(x2, inBetween, colors(correct(num),:),'LineStyle','none');
    set(h,'facealpha',.4)
    plot(x,indeg(idx),'.','markersize',5,'color','k')
    
    ylabel('Inward degree')
    in_lim = max(indeg+std_indeg);
    plot([1500*correct(num),1500*correct(num)],[-5,150],'-','color',[0.5 0.5 0.5])
    xlim([0,10500])
    ylim([-5,150])
    count=count+length(idx);
    set(gca,'xtick',[])
    
    subplot(2,1,2)
    hold on
    idx=(v(correct(num))+1):v(correct(num))+(v(correct(num)+1)-v(correct(num)))/2;
    c1 = outdeg(idx)-std_outdeg(idx);
    c2 = outdeg(idx)+std_outdeg(idx);
    x = pos(idx)+1500*(num-1)-500;
    
    x2 = [x, fliplr(x)];
    inBetween = [c1, fliplr(c2)];
    h=fill(x2, inBetween, colors(correct(num),:),'LineStyle','none');
    plot(x,outdeg(idx),'.','markersize',5,'color','k')
    set(h,'facealpha',.4)
    
    xlabel('RC positions')
    ylabel('Outward degree')
    in_lim = max(outdeg+std_outdeg);
    plot([1500*correct(num),1500*correct(num)],[-5,160],'-','color',[0.5 0.5 0.5])
    xlim([0,10500])
    ylim([-5,160])
    count=count+length(idx);
    set(gca,'xtick',[])  

end

% figure(2)
% plot(indeg((v(4)+1):v(4+1)),outdeg((v(4)+1):v(4+1)),'.','markersize',15,'color',colors(4,:))
% 
% figure(3)
% plot(indeg((v(5)+1):v(5+1)),outdeg((v(5)+1):v(5+1)),'.','markersize',15,'color',colors(5,:))
% 
% indeg=indeg(v(5):v(6));
% outdeg=outdeg(v(5):v(6));
% mdl=polyfit(indeg,outdeg,1);
% hold on
% limit=max(indeg);
% plot([0 limit],[mdl(1)*[0 limit]+mdl(2)],'color','k','markersize',5,'linewidth',1)
% xlabel('in-degree')
% ylabel('out-degree')
% xlim([20,90])



