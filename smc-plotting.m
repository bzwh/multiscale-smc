clear
%%
pth = 'D:/Ben/Dropbox/ms-smc/mallorn_x/';
pth = 'D:/Ben/Dropbox/ms-smc/custard/';
pth = 'D:/Ben/Dropbox/ms-smc/mallorn-complex/';
pth = 'D:/Ben/Dropbox/ms-smc/mallorn-complex/';
pth = 'D:/Ben/Dropbox/ms-smc/outputs/';
%pth = './mallorn-/';
%pth = './outputs/';
pth = './mallorn-wfit/';
pth = './mallorn-wfm/';
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
m=2;
n=2;
nreg=5;
nparam = 7;%size(pall,2);
nrp = 3;%(size(pars{1},2)-nparam)/nreg;
%       Sus    trans      Ker       Delay  DCs  
kmin = [0  ,0 ,0   ,0   , 0  ,1     ,0 ,0, 1.5,0.0, 0.75];
kmax = [250,30,1e-4,1e-5, 2.5,12.5  ,15,5, 5.5,0.025, 1];
lbls = {'Sus_{cow}','Trans_{cow}','Trans_{sheep}','K_r','K_a','Delay_\mu','Delay_\theta','F_A','F_B','f'};
if nreg>1
  rlbl = {'Cumbria','Devon','Wales','Scotland','RoUK'};
else
  rlbl = {'',' '};
end
cc=parula(nrnd);
cfd = 0.25;
lw = 2;

for reg=1:nreg    % Figure for each region separately
  figure('outerposition',[0,0,1920,1080],'visible','off')
  for rn=1:nrnd   % Plot each round on all subplots
    ids =1:size(pars{1},1);
    pall = pars{rn};
    % Plot posteriors
    for i=1:nrp
      subplot(m,n,i);
      hold all
      [k,fd] = mykerest(pall(ids,i+(reg-1)*3));%,kmin(i+1),kmax(i+1));
      l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
      l.Color(4)=(1-cfd)*rn/nrnd+cfd;
    end
  end
  % Actual particles for final round and labels
  for i=1:nrp
    subplot(m,n,i)
    hold all
    l=plot(pars{nrnd}(:,i+(reg-1)*3),abs(normrnd(0,1e-5*range(pars{nrnd}(:,i)),[N,1])),'.','MarkerSize',10,'Color',[0.1,0.1,0.1,0.01]);
    xlabel(lbls(i));
  end

  subplot(m,n,m*n)
  hold all
  lg = {};
  for i=1:nrnd
    plot(0,i,'color',cc(i,:),'linewidth',2)
    lg{i}=num2str(i);
  end
  axis off
  legend(lg)
  saveas(gcf,[pth 'pars-dens_' rlbl{reg} '.png'])
  close()
end

figure('outerposition',[0,0,1920,1080],'visible','off')
if nparam==2
  mm=2; nn=2;
else
  mm=3; nn=3;
end
for i=1:nparam
  subplot(mm,nn,i)
  hold all
  for rn=1:nrnd
    pall = pars{rn};
    [k,fd] = mykerest(pall(:,i+nreg*3));%,kmin(i+1),kmax(i+1));
    l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
    l.Color(4)=(1-cfd)*rn/nrnd+cfd;
  end
end
for i=1:nparam
  subplot(mm,nn,i)
  hold all
  l=plot(pars{nrnd}(:,i+nreg*3),abs(normrnd(0,1e-5*range(pars{nrnd}(:,i)),[N,1])),'.','MarkerSize',10,'Color',[0.1,0.1,0.1,0.01]);
  xlabel(lbls(i+3));
end
subplot(mm,nn,mm*nn)
hold all
lg = {};
for i=1:nrnd
  plot(0,i,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
axis off
legend(lg)
saveas(gcf,[pth 'pars-dens_others.png'])
close()
%% WFM
m=2;n=2;
II=[1,2,3,4];
prc = load('../../data/uk/wfm-cows.txt');
figure('outerposition',[0,0,1920,1080],'visible','off')
lbls = {'kE','mE','kI','mI'};
for i=1:4
  subplot(m,n,II(i))
  hold all
  for rn=1:nrnd
    pall = pars{rn};
    [k,fd] = mykerest(pall(:,end-13+i*2));%,kmin(i+1),kmax(i+1));
    l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
    l.Color(4)=(1-cfd)*rn/nrnd+cfd;
    plot(k,gampdf(k,prc(i,1),prc(i,2)),'c:')
  end
  xlabel(lbls{i})
end
subplot(m,n,m*n)
hold all
lg = {};
for i=1:nrnd
  %plot(0,i,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
%axis off
legend(lg)
saveas(gcf,[pth 'pars-dens_wfmc.png'])
close()
% Sheep--------------------------------------------------------------------
prs = load('../../data/uk/wfm-lamb.txt');
figure('outerposition',[0,0,1920,1080],'visible','off')
lbls = {'kE','mE','kI','mI'};
for i=1:4
  subplot(m,n,II(i))
  hold all
  for rn=1:nrnd
    pall = pars{rn};
    [k,fd] = mykerest(pall(:,end-12+i*2));%,kmin(i+1),kmax(i+1));
    l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
    l.Color(4)=(1-cfd)*rn/nrnd+cfd;
    plot(k,gampdf(k,prs(i,1),prs(i,2)),'c:')
  end
  xlabel(lbls{i})
end
subplot(m,n,m*n)
hold all
lg = {};
for i=1:nrnd
  %plot(0,i,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
%axis off
legend(lg)
saveas(gcf,[pth 'pars-dens_wfms.png'])
close()
% Transmission rate--------------------------------------------------------
prb = [prc(end,:) ; [10,0.006];[10,0.006];prs(end,:)];
figure('outerposition',[0,0,1920,1080],'visible','off')
lbls={'Cow to cow','Sheep to cow','Cow to sheep','Sheep to sheep'};
for i=1:4
  subplot(m,n,II(i))
  hold all
  for rn=1:nrnd
    pall = pars{rn};
    [k,fd] = mykerest(pall(:,end-4+i));%,kmin(i+1),kmax(i+1));
    l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
    l.Color(4)=(1-cfd)*rn/nrnd+cfd;
    plot(k,gampdf(k,prb(i,1),prb(i,2)),'c:')
    xlabel(lbls{i});
  end
end
subplot(m,n,m*n)
hold all
lg = {};
for i=1:nrnd
  %plot(0,i,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
%axis off
legend(lg)
saveas(gcf,[pth 'pars-dens_wfmb.png'])
close()

%%
figure('outerposition',[0,0,1920,1080])
lbls = {'Sus_{cow}','Trans_{cow}','Trans_{sheep}','K_r','K_a','Delay_\mu','Delay_\theta','F_A','F_B','f'};
cc=parula(nrnd);
cfd = 0.25;
nparam = 5;%size(pall,2);
for rn=1:nrnd
  ids =1:size(pars{1},1);
  pall = pars{rn};
  lw = 2;
  t_ratio = 0;

  for i=1:nparam
    subplot(m,n,i);
    hold all
    [k,fd] = mykerest(pall(ids,i+5));%,kmin(i+1),kmax(i+1));
    l=plot(k,pdf(fd,k),'Color',cc(rn,:),'linewidth',lw);
    l.Color(4)=(1-cfd)*rn/nrnd+cfd;
  end
end

for i=1:nparam
  subplot(m,n,i)
  hold all
  l=plot(pars{nrnd}(:,i+5),abs(normrnd(0,1e-5*range(pars{nrnd}(:,i+5)),[N,1])),'.','MarkerSize',10,'Color',[0.1,0.1,0.1,0.01]);
  xlabel(lbls(i+5));
end

subplot(m,n,6)
hold all
lg = {};
for i=1:nrnd
  plot(0,i,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
axis off
legend(lg)
saveas(gcf,[pth 'pars-dens2.png'])
close()
%% Error distributions per metric(region)
lbls = {'Cumbria','Devon','Wales','Scotland','RoUK'};
m = 2;
n = 3;
figure('outerposition',[0,0,1920,1080],'Visible','off')
cfd = 0;
nb = 25;
cc = parula(nrnd);
for i=1:nmet
  subplot(m,n,i)
  hold all
  for r=nrnd:-1:1
    histogram(errs{r}(2:end,i),nb,'FaceColor',cc(r,:),'FaceAlpha',(1-cfd)*r/nrnd+cfd,'EdgeColor','none')
    %vline(errs{r}(1,m),'r-.')
    %[h,c] = hist(errs{r}(2:end,m));
    %plot(c,smooth(h/N),'Color',cc(r,:))
  end
  xlabel(lbls{i});
end
subplot(m,n,6)
hold all
lg = {};
for i=1:nrnd
  plot(0,i,'color',cc(i,:),'linewidth',2)
  lg{i}=num2str(i);
end
axis off
legend(lg)
saveas(gcf,[pth 'errs.png'])
close
%%
np = size(pars{nrnd},2);
N = size(pars{nrnd},1);
nb = 10;

lbls = {'Sus_{cow}','Trans_{cow}','Trans_{sheep}','K_r','K_a'};

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:np
  for j=1:np
    if i==j
      subplot(np,np,(j-1)*np+i)
      [h,c]=hist(pars{nrnd}(:,i),nb);
      plot(c,h/N);
      xlabel(lbls{i})
    elseif i>j
      subplot(np,np,(j-1)*np+i)
      if(0)
        plot(pars{nrnd}(:,j),pars{nrnd}(:,i),'.');
      else
        hh=hist3(pars{nrnd}(:,[j,i]),[nb,nb]);
        h = hh';
        h(size(hh,1) + 1, size(hh,2) + 1) = 0;
        h(h==0) = NaN;
        yb = linspace(min(pars{nrnd}(:,i)),max(pars{nrnd}(:,i)),size(hh,1)+1);
        xb = linspace(min(pars{nrnd}(:,j)),max(pars{nrnd}(:,j)),size(hh,1)+1);
        p=pcolor(xb,yb,h);
        set(p,'Edgecolor','none');
      end
    end
  end
end
saveas(gcf,[pth 'corrs.png'])


