

%Identify trials which are bad across channels (preferably
%already eliminated bad channels)
artifact_indicator_matrix=[];
SigAlignedD=[];
thr=3; %Setting a threshold
channels=size(CDS,1);
for ch_itr=1:channels
    r_aligned=detrend(squeeze(CDS(ch_itr,:,:)))';
    %%%%% uncomment if you want to plot the per channel trials
    %             imagesc(r_aligned')
    %              plot(r_aligned')
    %              pause(0.01)
    %              clf
    r_aligned_unbiased=r_aligned-nanmean(nanmean(r_aligned));
    r_norm=diag(r_aligned_unbiased*r_aligned_unbiased');
    r_norm_mean=nanmean(r_norm);
    r_norm_var=nanvar(r_norm);
    r_norm_std=sqrt(r_norm_var);
    artifact_indicator_matrix(ch_itr,:)=r_norm'>r_norm_mean+thr*r_norm_std;
    SigAlignedD(ch_itr,:,:)=r_aligned_unbiased;
    ch_itr
end

TriRem=1; %Sets the threshold to reject bad trials

%Use this to identify the trials to keep, cn
cn=find(sum(artifact_indicator_matrix)<=TriRem);

%Also can be done this way:
ValTrial=1:size(CDS,3); %Or pick the trials that of a particular condition.
Rem=[]; %trials to remove
badTrials=unique([find(sum(artifact_indicator_matrix)>=TriRem) ]);

if contains(FiName(1:end-4),'pig1_protocol3_trial1_yrle037_400ua_190610_174136')==1
badTrials=11:13;
end


if contains(FiName(1:end-4),'pig1_protocol3_trial1_yrle037_200ua_190610_174049')==1
badTrials=[1 7     8    10  11  12  13  14    16];
end

for KL=1:length(badTrials)
    if isempty(find(ValTrial==badTrials(KL)))==0
        Rem=[Rem;badTrials(KL) KL find(ValTrial==badTrials(KL))];
    end
end

%trials to remove
if isempty(Rem)==0
    ValTrial(Rem(:,3))=[];
end


%%
% clf
% for ch=1:size(CDS,1)
% imagesc(squeeze(CDS(ch,:,ValTrial)))
% % plot(squeeze(CDS(ch,:,ValTrial)))
% pause
% clf
% end
