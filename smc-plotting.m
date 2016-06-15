clear
%%
pth = 'c:/users/huben/Dropbox/ms-smc/outputs/';
%pth = 'D:/Ben/projects/multiscale-smc/outputs/';
%pth = './outputs/';
%pth = 'n:/mathematical Biology/private/staff folders/ben/outputs/ms-smc/outputs/';
dfiles = dir([pth 'dat*txt' ]);
pars = {};
for i=1:size(dfiles,1)
  j = 1+str2num(dfiles(i).name(4:end-3));
  pars{j} = load([pth dfiles(i).name]);
end
efiles = dir([pth 'err*txt' ]);
errs = {};
for i=1:size(efiles,1)
  j = 1+str2num(efiles(i).name(4:end-3));
  errs{j} = load([pth efiles(i).name]);
end
N=size(pars{1},1);
npar = size(pars{1},2);
nrnd = size(pars,2)-1; % last file is still running?
nmet = size(errs{1},2);
wht = load([pth 'wht.txt']);
  
%% Particle distributions
%       Sus    trans      Ker       Delay  DCs  
kmin = [0  ,0 ,0   ,0   , 0  ,1     ,0 ,0, 1.5,0.0, 0.75];
kmax = [250,30,1e-4,1e-5, 2.5,12.5  ,15,5, 5.5,0.025, 1];
figure('units','normalized','outerposition',[0 0 1 1])
cc=parula(nrnd);
cfd = 0.25;
for rn=1:nrnd
  ids =1:size(pars{1},1);
  pall = pars{rn};
  lbls = {'Sus_{cow}','Trans_{cow}','Trans_{sheep}','K_r','K_a','Delay_\mu','Delay_\theta','F_A','F_B','f'};
  lw = 2;
  t_ratio = 0;
  nparam = size(pall,2);

  if (t_ratio)
    subplot(2,3,2)
    hold on
    [k,fd] = mykerest(pall(ids,2)./pall(ids,3),kmin(1),kmax(1));
    plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw)
    plot(pall(ids,2)./pall(ids,3),1e-6*rand([1,min(sum(ids),length(ids))]),'.','MarkerSize',5)
    xlabel('T_{cow}/T_{sheep}');
  end

  for i=1:nparam
    subplot(2,3,i);
    hold all
    [k,fd] = mykerest(pall(ids,i));%,kmin(i+1),kmax(i+1));
    l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
    l.Color(4)=(1-cfd)*rn/nrnd+cfd;
  end
end

for i=1:nparam
  subplot(2,3,i)
  hold all
  l=plot(pars{nrnd}(:,i),abs(normrnd(0,1e-5*range(pars{nrnd}(:,i)),[N,1])),'.','MarkerSize',10,'Color',[0.1,0.1,0.1,0.01]);
  xlabel(lbls(i));
end

subplot(2,3,6)
hold all
lg = {};
for i=1:nrnd
  plot(0:1,i:i+1,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
legend(lg)
%% Error distributions per metric
figure
cfd = 0;
nb = 25;
cc = parula(nrnd);
for m=1:nmet
  subplot(2,3,m)
  hold all
  for r=2:nrnd-1
    histogram(errs{r}(2:end,m),nb,'FaceColor',cc(r,:),'FaceAlpha',(1-cfd)*r/nrnd+cfd,'EdgeColor','none')
    %vline(errs{r}(1,m),'r-.')
    %[h,c] = hist(errs{r}(2:end,m));
    %plot(c,smooth(h/N),'Color',cc(r,:))
  end
end
%%
