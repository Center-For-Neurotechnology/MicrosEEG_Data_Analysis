clear
addpath(genpath('Y:\StimDataBackup\Code_Stimulation'))
DIRIntan='Y:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\';
load('Y:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\ExtractedMatlabData_pig1_protocol3_trial1_yrle037_6000ua_190610_174636\TrialsStim.mat')

addpath(genpath(['Y:\IntraOp_Micro\Papers_Published_Inprogress\PEDOTPtNRDepthLaminarPaper\Code\']))
RHDDir=dir([DIRIntan,'*.rhd']);
% ChanImpedances=zeros(64,1);
TrialsData=[];TrialsDataPow=[];
TrialsDataPow2=[];
TrialsDataPow3=[];
TrialsDataPow4=[];
TrialsDataPow5=[];
TrialsDataPow6=[];
TrialsDataPow=[];
TrialsDatacsd=[];TrialsDatacsdmua=[];CDSAll=[];TrialsAllAcross=[];Pow6All=[];PowMUAAll=[];
ChanImpedances2=zeros(64,1);Pow4All=[];Pow5All=[];Pow3All=[];
for siS=1:length(RHDDir)
    FiName=RHDDir(siS).name;
    Pathfi=[DIRIntan,'\ExtractedMatlabData_',FiName(1:end-4)];
    DIRIntanInd=[DIRIntan,'\ExtractedMatlabData_',FiName(1:end-4)];

    load([DIRIntanInd,'\SavedTrials'],'CDSraw','Pow6','Pow5',...
        'Pow4','Pow3','Pow2','Pow1','timeoi','CDS', ...
        'Trials','CSDTri','CSDmuaTri', ...
        'freqoi','MUAPowTri')
    load([DIRIntanInd,'\CleanedFilteredData'], 'TrialsAll','ChanSel')

    TrialRejectionEnergy

    if siS==1
        CDSAll(:,:,1:size(CDS(:,:,ValTrial),3))=CDS(:,:,ValTrial);
        TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co=size(CDS(:,:,ValTrial),3)

    elseif siS>1
        CDSAll(:,:,co+(1:size(CDS(:,:,ValTrial),3)))=CDS(:,:,ValTrial);
        TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co=size(CDSAll,3)
        %             pause
    end

    if siS==1
        Pow6All(:,:,1:size(Pow6(:,:,ValTrial),3))=Pow6(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co6=size(Pow6(:,:,ValTrial),3)

    elseif siS>1
        Pow6All(:,:,co6+(1:size(Pow6(:,:,ValTrial),3)))=Pow6(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co6=size(Pow6All,3)
        %             pause
    end

    if siS==1
        Pow5All(:,:,1:size(Pow5(:,:,ValTrial),3))=Pow5(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co5=size(Pow5(:,:,ValTrial),3)

    elseif siS>1
        Pow5All(:,:,co5+(1:size(Pow5(:,:,ValTrial),3)))=Pow5(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co5=size(Pow5All,3)
        %             pause
    end

    if siS==1
        Pow4All(:,:,1:size(Pow4(:,:,ValTrial),3))=Pow4(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co4=size(Pow4(:,:,ValTrial),3)

    elseif siS>1
        Pow4All(:,:,co4+(1:size(Pow4(:,:,ValTrial),3)))=Pow4(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        co4=size(Pow4All,3)
        %             pause
    end
    if siS==1
        PowMUAAll(:,:,1:size(MUAPowTri(:,:,ValTrial),3))=MUAPowTri(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        com=size(MUAPowTri(:,:,ValTrial),3)

    elseif siS>1
        PowMUAAll(:,:,com+(1:size(MUAPowTri(:,:,ValTrial),3)))=MUAPowTri(:,:,ValTrial);
        %             TrialsAllAcross=[TrialsAllAcross;TrialsAll(ValTrial,:)];
        com=size(PowMUAAll,3)
        %             pause
    end

    TrialsData=[TrialsData; repmat([siS TrialsAll(1,2)],size(CDS,1),1) ...
        (1:size(MUAPowTri,1))' squeeze(nanmean(CDS(:,:,ValTrial),3)) squeeze(nanstd(CDS(:,:,ValTrial),[],3))];

    TrialsDatacsd=[TrialsDatacsd; repmat([siS TrialsAll(1,2)],size(CSDTri,1),1) ...
        (1:size(CSDTri,1))' squeeze(nanmean(CSDTri(:,:,ValTrial),3)) squeeze(nanstd(CSDTri(:,:,ValTrial),[],3))];

    TrialsDatacsdmua=[TrialsDatacsdmua; repmat([siS TrialsAll(1,2)],size(CSDmuaTri,1),1) ...
        (1:size(CSDmuaTri,1))' squeeze(nanmean(CSDmuaTri(:,:,ValTrial),3)) squeeze(nanstd(CSDmuaTri(:,:,ValTrial),[],3))];

    TrialsDataPow=[TrialsDataPow; repmat([siS TrialsAll(1,2)],size(MUAPowTri,1),1) ...
        (1:size(MUAPowTri,1))' squeeze(nanmean(MUAPowTri(:,:,ValTrial),3)) squeeze(nanstd(MUAPowTri(:,:,ValTrial),[],3))];

    TrialsDataPow2=[TrialsDataPow2; repmat([siS TrialsAll(1,2)],size(Pow2,1),1) ...
        (1:size(Pow2,1))' squeeze(nanmean(Pow2(:,:,ValTrial),3)) squeeze(nanstd(Pow2(:,:,ValTrial),[],3))];

    TrialsDataPow3=[TrialsDataPow3; repmat([siS TrialsAll(1,2)],size(Pow3,1),1) ...
        (1:size(Pow3,1))' squeeze(nanmean(Pow3(:,:,ValTrial),3)) squeeze(nanstd(Pow3(:,:,ValTrial),[],3))];

    TrialsDataPow4=[TrialsDataPow4; repmat([siS TrialsAll(1,2)],size(Pow4,1),1) ...
        (1:size(Pow4,1))' squeeze(nanmean(Pow4(:,:,ValTrial),3)) squeeze(nanstd(Pow4(:,:,ValTrial),[],3))];

    TrialsDataPow5=[TrialsDataPow5; repmat([siS TrialsAll(1,2)],size(Pow5,1),1) ...
        (1:size(Pow5,1))' squeeze(nanmean(Pow5(:,:,ValTrial),3)) squeeze(nanstd(Pow5(:,:,ValTrial),[],3))];

    TrialsDataPow6=[TrialsDataPow6; repmat([siS TrialsAll(1,2)],size(Pow6,1),1) ...
        (1:size(Pow6,1))' squeeze(nanmean(Pow6(:,:,ValTrial),3)) squeeze(nanstd(Pow6(:,:,ValTrial),[],3))];
    siS
end
%%
baselineRange=1:950;MeasureRange=1010:2000;FS=1000;
[PVLAllcds,h6cds,PVLcds,h3cds,Meanpowcds]=StatsCheckValues(CDSAll,TrialsAllAcross,baselineRange,MeasureRange,FS);

save(['Z:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\StatsInfoCDS'],...
    'PVLAllcds','h6cds','PVLcds','h3cds','CDSAll','TrialsAllAcross','Meanpowcds')

baselineRange=1:45;MeasureRange=51:200;FS=100;
[PVLAllpow6,h6pow6,PVLpow6,h3pow6,Meanpow6]=StatsCheckValues(Pow6All,TrialsAllAcross,baselineRange,MeasureRange,FS);

save(['Z:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\StatsInfoPow6'],...
    'PVLAllpow6','h6pow6','PVLpow6','h3pow6','Pow6All','TrialsAllAcross','Meanpow6')

baselineRange=1:95;MeasureRange=102:400;FS=200;
[PVLAllmua,h6mua,PVLmua,h3mua,Meanpowmua]=StatsCheckValues(PowMUAAll,TrialsAllAcross,baselineRange,MeasureRange,FS);

save(['Y:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\StatsInfoCDS'],...
    'PVLAllmua','h6mua','PVLmua','h3mua','PowMUAAll','TrialsAllAcross','Meanpowmua')


baselineRange=1:45;MeasureRange=51:200;FS=100;
[PVLAllpow4,h6pow4,PVLpow4,h3pow4,Meanpow4]=StatsCheckValues(Pow4All,TrialsAllAcross,baselineRange,MeasureRange,FS);

save(['Y:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\StatsInfoPow4'],...
    'PVLAllpow4','h6pow4','PVLpow4','h3pow4','Pow4All','TrialsAllAcross','Meanpow4')

baselineRange=1:45;MeasureRange=51:200;FS=100;
[PVLAllpow5,h6pow5,PVLpow5,h3pow5,Meanpow5]=StatsCheckValues(Pow5All,TrialsAllAcross,baselineRange,MeasureRange,FS);

save(['Y:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\StatsInfoPow5'],...
    'PVLAllpow5','h6pow5','PVLpow5','h3pow5','Pow5All','TrialsAllAcross','Meanpow5')

pause
%%
%
% [ME,SE,gname]=grpstats(TrialsData(:,4:end),[TrialsData(:,2) TrialsData(:,3)],{'nanmean','sem','gname'});
%
% current=str2num(char(gname(:,1)));
% chan=str2num(char(gname(:,2)));
% PowPlot=TrialsDatacsd; fs=30000; STEP=.5; Signal='CSD';


clf
addpath(genpath('Y:\Projects\Lab_Materials\Analysis_Tools_and_Software\ExampleCodeRunning\shadedErrorBar\'))

addpath(genpath('Y:\Projects\Lab_Materials\Analysis_Tools_and_Software\ExampleCodeRunning\bipolar_colormap'))
colormap(bipolar2(255))
curr=unique(TrialsDataPow3(:,2) );
sBT=.9;
% PowPlot=TrialsData;fs=1000; STEP=2; Signal='LFP';PVplot2=PVLcds;
% PowPlot=TrialsDataPow; fs=200; STEP=.5; Signal='MUA power'; PVplot2=PVLmua;
% PowPlot=TrialsDataPow6; fs=100; STEP=.5; Signal='high gamma power'; PVplot2=PVLpow6;
% PowPlot=TrialsDataPow5; fs=100; STEP=.5; Signal='gamma power'; PVplot2=PVLpow5;
PowPlot=TrialsDataPow4; fs=100; STEP=.5; Signal='beta power'; PVplot2=PVLpow4;
hs=[];
hs=PVplot2;
hs(PVplot2<=0.001)=1;
hs(PVplot2>0.001)=NaN;
for cu=1:8
    FI=find(PowPlot(:,2)==curr(cu));
    PlotV=PowPlot(FI,6:end);
    % plotsle=TrialsDataPow(FI,5:end);
    subplot(2,8,cu)
    hold on
    %     plot(-1+(1:300)/100,PowPlot(FI,5:end)')
    %     plot(-1+(1:600)/200,PowPlot(FI,5:end)')
    %     plot(-1+(1:size(PlotV,2)/2)/fs,...
    %         PlotV(:,(1:size(PlotV,2)/2))'+repmat(STEP*PowPlot(FI,3),1,size(PlotV(:,(1:size(PlotV,2)/2)),2))','color',[0 0 0 .5])
    for ch1=1:22
        hpn=squeeze(hs(ch1,2:end,cu));
        hpn(-1+(1:size(PlotV,2)/2)/fs<0021)=NaN;

        plot(-sBT+(1:size(PlotV,2)/2)'/fs,...
            PlotV(ch1,(1:size(PlotV,2)/2))'.*hpn'+ ...
            repmat(STEP*(ch1),1,size(PlotV(ch1,(1:size(PlotV,2)/2)),2))', ...
            '.','color',[0 .8 0],'markersize',20)
        shadedErrorBar(-sBT+(1:size(PlotV,2)/2)'/fs,...
            PlotV(ch1,(1:size(PlotV,2)/2))'+repmat(STEP*(ch1),1,size(PlotV(ch1,(1:size(PlotV,2)/2)),2))', ...
            PlotV(ch1,(size(PlotV,2)/2+1:end))'/sqrt(10), ...
            {'color',[0 0 0],'linewidth',2},.4)
    end
    %     plot(-1+(1:90000)/30000,PowPlot(FI,5:end)')

    xlim([-.25 .75])
    ylim([1 STEP*(ch1+1)+STEP*2])
    set(gca,'ytick',STEP*(2+(1:length(ChanSel))),'yticklabel',num2str(ChanSel'))
    %      ylim([-5 5])
    %      ylim([-10 10])
    %      ylim([.8 1.2])
    if cu==1
        ylabel('channel')
        xlabel('time (sec)')
    end
    box off
    xline(0,'linestyle',':','color','k','linewidth',2);
    title([num2str(curr(cu)),' microA'])
    subplot(2,8,cu+8)
    imagesc(-sBT+(1:size(PlotV,2))/fs,1:22,PowPlot(FI,5:end))
    %     caxis([-5 5])
    caxis([0.8 1.2])
    %     caxis([-2 2])
    xlim([-.25 .75])
    axis xy
    xline(0,'linestyle',':','color','w','linewidth',2);
    if cu==1
        ylabel('channel')
        xlabel('time (sec)')
    end
    title(Signal)
    box off
    % pause
    % clf
end

%  print(gcf,'-dpng','-r600',['Z:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\ResponsesAlongDepth',Signal])

%%
Channel=1;
clf
curr=unique(TrialsDataPow3(:,2) );
STEP=.5;
STEP=0;
bls=1;
% PowPlot=TrialsDatacsd; fs=30000; STEP=.5; Signal='CSD';bls=1;

COLST=colormap(hsv(8));
PowPlot=TrialsData;fs=1000; STEP=2; Signal='LFP';PVplot=PVLAllcds;bls=1;
% PowPlot=TrialsDataPow; fs=200; STEP=.5; Signal='MUA power'; PVplot=PVLAllmua;bls=1;
% PowPlot=TrialsDataPow6; fs=100; STEP=0.25; Signal='high gamma power'; PVplot=PVLAllpow6;bls=1;
% PowPlot=TrialsDataPow5; fs=100; STEP=.5; Signal='gamma power'; PVplot=PVLAllpow5;bls=1;
% PowPlot=TrialsDataPow4; fs=100; STEP=.5; Signal='beta power'; PVplot=PVLAllpow4;bls=1;
hs=[];
hs=PVplot;
hs(PVplot<=0.01)=1;
hs(PVplot>0.01)=NaN;
for Channel=1:22
    FI=find(PowPlot(:,3)==Channel);
    PlotV=PowPlot(FI,6:end);

    currentSteps=PowPlot(FI,2);
    % plotsle=TrialsDataPow(FI,5:end);
    subplot(3,8,Channel)
    hold on
    %     plot(-1+(1:300)/100,PowPlot(FI,5:end)')
    %     plot(-1+(1:600)/200,PowPlot(FI,5:end)')
    %     plot(-1+(1:size(PlotV,2)/2)/fs,...
    %         PlotV(:,(1:size(PlotV,2)/2))'+repmat(STEP*PowPlot(FI,3),1,size(PlotV(:,(1:size(PlotV,2)/2)),2))','color',[0 0 0 .5])
    hsC=hs(Channel,2:end);
    hsC(-1+(1:size(PlotV,2)/2)/fs<.021)=NaN;
    stem(-bls+(1:size(PlotV,2)/2)'/fs,STEP*8+5+hsC', ...
        '.','linewidth',2,'color',[.8 .8 .8])

    for currst=1:8
        fiC=find(currentSteps==curr(currst));
        shadedErrorBar(-bls+(1:size(PlotV,2)/2)'/fs,...
            PlotV(fiC,(1:size(PlotV,2)/2))'+ ...
            repmat(STEP*currst,1,size(PlotV(fiC,(1:size(PlotV,2)/2)),2))', ...
            PlotV(fiC,(size(PlotV,2)/2+1:end))'/sqrt(10), ...
            {'color',COLST(currst,:),'linewidth',2},.5)
        yline(STEP*currst+bls-1,'linestyle',':','color','k','linewidth',1);
    end

    xlim([-.25 .75])
    ylim([.5 STEP*currst+2*STEP+1])
    set(gca,'linewidth',2)
    %     set(gca,'ytick',STEP*PowPlot(FI,3),'yticklabel',num2str((1:22)'))
    %      ylim([-5 5])
    %      ylim([-10 10])
    %          ylim([.8 1.2])
    if cu==1
        ylabel('channel')
        xlabel('time (sec)')
    end
    box off
    xline(0,'linestyle',':','color','k','linewidth',1);

    title([num2str(ChanSel(Channel)),' chan'])

    %     subplot(2,8,cu+8)
    %     imagesc(-1+(1:size(PlotV,2))/fs,1:22,PowPlot(FI,5:end))
    %     %     caxis([-5 5])
    %     caxis([0.8 1.2])
    %     % caxis([-2 2])
    %     xlim([-.5 1])
    %     axis xy
    %     xline(0,'linestyle',':','color','w','linewidth',2);
    %     if cu==1
    %         ylabel('channel')
    %         xlabel('time (sec)')
    %     end
    %     title(Signal)
    %     box off
    % pause
    % clf
end

%  print(gcf,'-dpng','-r600',['Z:\IntraOp_Micro\MicroPigData\061019_Pig case\Recording files\ResponsesAcrossChannelsStats',Signal])



%%
Channel=1;
clf
curr=unique(TrialsDataPow3(:,2) );
STEP=.5;
STEP=0;
bls=1;

COLST=colormap(hsv(8));
PowPlot=TrialsData;fs=1000; STEP=2; Signal='LFP';PVplot=PVLAllcds;bls=1;
for mwn=1:5
    if mwn==1
        MeanMp=Meanpowcds;
    elseif mwn==2
        MeanMp=Meanpow4;
    elseif mwn==3
        MeanMp=Meanpow5;
    elseif mwn==4
        MeanMp=Meanpow6;
    elseif mwn==5
        MeanMp=Meanpowmua;
    end
    colormap(bipolar2(255))
%     colormap(cool(255))
    [ME,SE,gname]=grpstats(MeanMp(:,8:9),[MeanMp(:,1) MeanMp(:,2)],{'mean','sem','gname'});
    chem=unique(str2num(char(gname(:,1))));
    currw=unique(str2num(char(gname(:,2))));
    subplot(2,3,mwn)
    MTriw=[];MTriwst=[];
    for CW=1:length(currw)
        for ch=1:length(chem)
            FI=find(MeanMp(:,1)==chem(ch) & MeanMp(:,2)==currw(CW));
            MTriw(CW,ch)=nanmean(MeanMp(FI,8),1);
             MTriwst(CW,ch)=nanstd(MeanMp(FI,8),1)/sqrt(10);
        end
    end
%     MTriw=[[MTriw [NaN*ones(size(MTriw,1),1)]];NaN*ones(1,size(MTriw,2)+1)]
    pcolor(1:length(currw)+0,1:length(chem)+0,MTriw')
    shading flat
    colorbar
    if mwn==1
    caxis([-6 3])
    end
    if mwn>1
    caxis([-1 1])
    end

%     caxis([0.1 .9])

end

%%
Channel=1;
clf
curr=unique(TrialsDataPow3(:,2) );
STEP=.5;
STEP=0;
bls=1;

COLST=colormap(hsv(8));
COLChan=colormap(parula(22))*.9;

PowPlot=TrialsData;fs=1000; STEP=2; Signal='LFP';PVplot=PVLAllcds;bls=1;
for mwn=1:5
    if mwn==1
        MeanMp=Meanpowcds;
    elseif mwn==2
        MeanMp=Meanpow4;
    elseif mwn==3
        MeanMp=Meanpow5;
    elseif mwn==4
        MeanMp=Meanpow6;
    elseif mwn==5
        MeanMp=Meanpowmua;
    end
    colormap(bipolar2(255))
%     colormap(cool(255))
    [ME,SE,gname]=grpstats(MeanMp(:,8:9),[MeanMp(:,1) MeanMp(:,2)],{'mean','sem','gname'});
    chem=unique(str2num(char(gname(:,1))));
    currw=unique(str2num(char(gname(:,2))));
    subplot(2,3,mwn)
    MTriw=[];MTriwst=[];
    for CW=1:length(currw)
        for ch=1:length(chem)
            FI=find(MeanMp(:,1)==chem(ch) & MeanMp(:,2)==currw(CW));
            MTriw(CW,ch)=nanmean(MeanMp(FI,9),1);
             MTriwst(CW,ch)=nanstd(MeanMp(FI,9),1)/sqrt(10);
        end
    end

%         for CW=1:length(currw)
%             FI=find(MeanMp(:,1)>1 & MeanMp(:,2)==currw(CW));
% MN=MeanMp(FI,:);
% [p,h,stats]=kruskalwallis(MN(:,8),MN(:,1),'off');
%  multcompare(stats)
% p
% pause 
% clf
%         end

%     MTriw=[[MTriw [NaN*ones(size(MTriw,1),1)]];NaN*ones(1,size(MTriw,2)+1)]
%     pcolor(1:length(currw)+0,1:length(chem)+0,MTriw')
%     shading flat
%     colorbar
%     if mwn==1
%     caxis([-6 3])
%     end
%     if mwn>1
%     caxis([-1 1])
%     end
% for ch=2:22
% errorbar(MTriw(:,ch),MTriwst(:,ch),'CapSize',0,'color',COLChan(ch,:),'linewidth',2)
% hold on
% end
% box off

for cu=1:8
errorbar(MTriw(cu,:),1:22,MTriwst(cu,:), ...
    'horizontal',...
    'CapSize',0,'color',COLST(cu,:),'linewidth',2)
hold on
end
box off
xlim([0 1.0])
%     caxis([0.1 .9])

end

% subplot(1,2,2)
% MTriw=[];
% for CW=1:length(currw)
%     for ch=1:length(chem)
%      FI=find(MeanMp(:,1)==chem(ch) & MeanMp(:,2)==currw(CW));   
% MTriw(CW,ch)=nanmean(MeanMp(FI,9),1);
%     end
% end
% MTriw=[[MTriw [NaN*ones(size(MTriw,1),1)]];NaN*ones(1,size(MTriw,2)+1)]
% pcolor(1:length(currw)+1,1:length(chem)+1,MTriw')
% shading flat
% colorbar

% for Channel=1:22
%     FI=find(MeanMp(:,1)==Channel);
%     PlotV=MeanMp(FI,6:end);
% 
%     currentSteps=PowPlot(FI,2);
%     hsC=hs(Channel,2:end);
%     hsC(-1+(1:size(PlotV,2)/2)/fs<.021)=NaN;
%     %    stem(-bls+(1:size(PlotV,2)/2)'/fs,STEP*currst+5+hsC', ...
%     %        '.','linewidth',2,'color',[.8 .8 .8])
% 
%     for currst=1:8
%         %         fiC=find(currentSteps==curr(currst));
% Meanpowcds
%     end
% 
% end
