function NEAR_ICA_step2_NF(study_info, age)

pipeline='NEARICA_NF'; 

addpath('ADJUST')
addpath('NEAR_ChannelRejection')
ext='.set';

% enter all the event/condition markers
task_event_markers = {'go','grsp'};

% epoch length in second
task_epoch_length = [-3 3];

% lower and upper threshold (in mV)
volt_threshold = [-150 150];

% enter the list of frontal channels to check
frontal_channels = {'E1', 'E8', 'E14', 'E17','E21','E25','E32'};
% recommended list for EGI 128 channel net

%params.isVisIns = 0; % set to 1 if you want to visualize intermediate cleaning of NEAR Cleaning (bad channels + bad segments)
params.isInterp = 1; % set to 1 if you want to interpolate the removed bad channels (by Spherical Interpolation)
params.isAvg    = 1; % set to 1 if you want to perform average referencing

% Segmentation using fixation intervals - parameters begin %
% N.B: The following parameters can be set to [] if params.isSegt = 0
%params.sname = 'segt_visual_attention.xlsx'; % the visual segmentation coding file
%params.sloc  = []; % location of the xlsx file
%params.look_thr = 4999; % consider only the segments that exceed this threshold+1 in ms to retain; alternatively can be set to [] if no thresholding is preferred
% Segmentation using fixation intervals - parameters end %

% Parameters for NEAR- Bad Segments Correction/Rejection using ASR end %
interp_type = 'spherical'; % other values can be 'v4'. Please refer to pop_interp.m for more details.

%% Initialize output variables
subj_ids={};
subj_ages={};
lof_flat_channels={};
lof_channels={};
lof_periodo_channels={};
% Bad channels identified using LOF
lof_bad_channels={};
% number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
ica_preparation_bad_channels=[];
% length of data (in second) fed into ICA decomposition
length_ica_data=[];
% total independent components (ICs)
total_ICs=[];
% number of artifacted ICs
ICs_removed=[];
% number of epochs before artifact rejection
total_epochs_before_artifact_rejection=[];
% number of epochs after artifact rejection
total_epochs_after_artifact_rejection=[];
% total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
total_channels_interpolated=[];
asr_tot_samples_modified=[];
asr_change_in_RMS=[];

subj_idx=1;

%% Loop over all data files
for a_idx=53:53%size(study_info.participant_info,1)
    subj_age = study_info.participant_info.age{a_idx};
    if strcmp(subj_age,age)
    
        % Get subject ID from study info
        subj_id=study_info.participant_info.participant_id{a_idx};

        subj_ids{subj_idx}=subj_id;
        subj_ages{subj_idx}=subj_age;
        
        % Path containing subject data
        subject_dir=fullfile(study_info.data_dir,'data', 'derivatives', pipeline, age, subj_id);
        subject_data_dir=fullfile(subject_dir, '01_filtered_data');
        subject_output_data_dir=fullfile(study_info.data_dir, 'data','derivatives', pipeline, age, subj_id);
        data_file_name=sprintf('%s_task-%s_eeg.set',subj_id, study_info.task);

        fname=sprintf('%s_task-grasping_eeg_filtered_data.set',subj_id);

        % Load data
        EEG=pop_loadset('filepath', fullfile(subject_data_dir),...
            'filename', fname);
        origEEG=EEG;

        %zapline
        load(fullfile(subject_dir,sprintf('%s_zapline_continuous_old.mat',subj_id)));
        data=data_continuous_old.*1e6;
        EEG.data=data';%add'forold

        
        %     fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
        %     saveas(fig, fullfile(subject_dir,'zapped_psd.png'));

        % Parameters for NEAR - Bad Channels Detection begin %
        isFlat        = 1;
        flatWin       = 5; % tolerance level in s(default: 5)
        isLOF         = 1;
        dist_metric = 'seuclidean'; % Distance metric to compute k-distance
        thresh_lof    = 2.5; % Threshold cut-off for outlier detection on LOF scores
        isAdapt       = 10; % The threshold will be incremented by a factor of 1 if the given threshold detects more than of total channels (eg., 10); if this variable left empty [], no adaptive thresholding is enabled.
        isPeriodogram = 0; % flag variable to enable or disable periodogram method (default: 0)
        frange        = [1 20]; % Frequency Range in Hz
        winsize       = 1; % window length in s
        winov         = 0.66; % 66% overlap factor
        pthresh       = 4.5; % Threshold Factor to predict outliers on the computed energy

        rej_cutoff    = 13;   % A lower value implies severe removal (Recommended value range: 20 to 30)
        rej_mode      = 'off'; % Set to 'off' for ASR Correction and 'on for ASR Removal (default: 'on')
        add_reject    = 'off'; % Set to 'on' for additional rejection of bad segments if any after ASR processing (default: 'off')

        %% NEAR Bad Channel Detection
        [EEG, flat_ch, lof_ch, periodo_ch, LOF_vec] = NEAR_getBadChannels(EEG, isFlat, flatWin, isLOF, thresh_lof, dist_metric, isAdapt, ...
            isPeriodogram, frange, winsize, winov, pthresh, 0);
        save(fullfile(subject_output_data_dir, 'LOF_Values_old.mat'), 'LOF_vec'); % save .mat format
        disp('Bad Channel Detection is performed successfully');
        badChans = sort(unique(union(union(flat_ch, lof_ch),periodo_ch)));

        if(~isempty(badChans))
            if(size(badChans,1) ~= 1)
                badChans = badChans';
            end
        end

        EEG = pop_select(EEG, 'nochannel', badChans);

        lof_flat_channels{subj_idx}='';
        if numel(flat_ch)>0
            lof_flat_channels(subj_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(flat_ch,3), 'UniformOutput', false)',',');
        end
        lof_channels{subj_idx}='';
        if numel(lof_ch)>0
            lof_channels(subj_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(lof_ch,3), 'UniformOutput', false)',',');
        end
        lof_periodo_channels{subj_idx}='';
        if numel(periodo_ch)>0
            lof_periodo_channels(subj_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(periodo_ch,3), 'UniformOutput', false)',',');
        end
        lof_bad_channels{subj_idx}='';
        if numel(badChans)>0
            lof_bad_channels(subj_idx)=join(cellfun(@(x) num2str(x(1)), num2cell(badChans,3), 'UniformOutput', false)',',');
        end

        fig=figure();
        topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
        saveas(fig, fullfile(subject_output_data_dir,'06-lof_removed_old.png'));

        %% Bad epochs correction/removal using ASR
        EEG_copy = EEG;
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off', ...
            'Highpass','off','BurstCriterion',rej_cutoff,'WindowCriterion',add_reject,'BurstRejection',rej_mode,'Distance','Euclidian');

        if(strcmp(rej_mode, 'on'))
            modified_mask = ~EEG.etc.clean_sample_mask;
        else
            modified_mask = sum(abs(EEG_copy.data-EEG.data),1) > 1e-10;
        end

        tot_samples_modified = (length(find(modified_mask)) * 100) / EEG_copy.pnts;
        tot_samples_modified = round(tot_samples_modified * 100) / 100;
        asr_tot_samples_modified(subj_idx)=tot_samples_modified;
        change_in_RMS = -(mean(rms(EEG.data,2)) - mean(rms(EEG_copy.data,2))*100)/mean(rms(EEG_copy.data,2)); % in percentage
        change_in_RMS = round(change_in_RMS * 100) / 100;
        asr_change_in_RMS(subj_idx) =change_in_RMS;
        %if(isVisIns)
        %   try
        %       vis_artifacts(EEG,EEG_copy);
        %   catch
        %       warning('vis_artifacts failed. Skipping visualization.')
        %  end
        % end
        fprintf('\nArtifacted epochs are corrected by ASR algorithm\n');

        %% Save data after running ASR function, if saving interim results was preferred
        EEG = eeg_checkset( EEG );
        EEG = pop_editset(EEG, 'setname', strrep(data_file_name, ext, '_asr_data_old'));
        EEG = pop_saveset( EEG,'filename',strrep(data_file_name, ext, '_asr_data_old.set'),...
            'filepath', [subject_output_data_dir filesep '02_near_data' filesep]); % save .set format

        fig=compute_and_plot_psd(EEG,1:EEG.nbchan);
        saveas(fig, fullfile(subject_output_data_dir,'07-asr_psd_old.png'));

        %% STEP 8: Prepare data for ICA
        EEG_copy=EEG; % make a copy of the dataset
        EEG_copy = eeg_checkset(EEG_copy);

        % Perform 1Hz high pass filter on copied dataset
        transband = 1;
        fl_cutoff = transband/2;
        fl_order = 3.3 / (transband / EEG.srate);

        if mod(floor(fl_order),2) == 0
            fl_order=floor(fl_order);
        elseif mod(floor(fl_order),2) == 1
            fl_order=floor(fl_order)+1;
        end

        EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff,...
            'ftype', 'highpass', 'wtype', 'hamming', 'forder', fl_order,...
            'minphase', 0);
        EEG_copy = eeg_checkset(EEG_copy);

        % Create 1 second epoch
        % insert temporary marker 1 second apart and create epochs
        EEG_copy=eeg_regepochs(EEG_copy,'recurrence', 1, 'limits',[0 1],...
            'rmbase', [NaN], 'eventtype', '999');
        EEG_copy = eeg_checkset(EEG_copy);

        % Find bad epochs and delete them from dataset
        % [lower upper] threshold limit(s) in mV.
        vol_thrs = [-1000 1000];
        % [lower upper] threshold limit(s) in dB.
        %emg_thrs = [-100 30];
        % [lower upper] frequency limit(s) in Hz.
        %emg_freqs_limit = [20 40];

        % Find channel/s with xx% of artifacted 1-second epochs and delete them
        chanCounter = 1; ica_prep_badChans = [];
        numEpochs =EEG_copy.trials; % find the number of epochs
        all_bad_channels=0;

        for ch=1:EEG_copy.nbchan
            % Find artifaceted epochs by detecting outlier voltage
            EEG_copy = pop_eegthresh(EEG_copy,1, ch, vol_thrs(1), vol_thrs(2),...
                EEG_copy.xmin, EEG_copy.xmax, 0, 0);
            EEG_copy = eeg_checkset( EEG_copy );

            % 1         : data type (1: electrode, 0: component)
            % 0         : display with previously marked rejections? (0: no, 1: yes)
            % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)

            % Find artifaceted epochs by using thresholding of frequencies in the data.
            % this method mainly rejects muscle movement (EMG) artifacts
            %         EEG_copy = pop_rejspec( EEG_copy, 1,'elecrange',ch ,'method','fft',...
            %             'threshold', emg_thrs, 'freqlimits', emg_freqs_limit,...
            %             'eegplotplotallrej', 0, 'eegplotreject', 0);

            % method                : method to compute spectrum (fft)
            % threshold             : [lower upper] threshold limit(s) in dB.
            % freqlimits            : [lower upper] frequency limit(s) in Hz.
            % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
            % eegplotreject         : 0 = Do not reject marked trials (but store the  marks).

            % Find number of artifacted epochs
            EEG_copy = eeg_checkset( EEG_copy );
            EEG_copy = eeg_rejsuperpose( EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
            artifacted_epochs=EEG_copy.reject.rejglobal;

            % Find bad channel / channel with more than 20% artifacted epochs
            if sum(artifacted_epochs) > (numEpochs*20/100)
                ica_prep_badChans(chanCounter) = ch;
                chanCounter=chanCounter+1;
            end
        end

        % If all channels are bad, save the dataset at this stage and ignore the remaining of the preprocessing.
        if numel(ica_prep_badChans)==EEG.nbchan || numel(ica_prep_badChans)+1==EEG.nbchan
            all_bad_channels=1;
            warning(['No usable data for datafile', data_file_name]);
        else
            % Reject bad channel - channel with more than xx% artifacted epochs
            EEG_copy = pop_select( EEG_copy,'nochannel', ica_prep_badChans);
            EEG_copy = eeg_checkset(EEG_copy);
        end

        if numel(ica_prep_badChans)==0
            ica_preparation_bad_channels{subj_idx}='0';
        else
            ica_preparation_bad_channels{subj_idx}=num2str(ica_prep_badChans);
        end

        if all_bad_channels == 1
            length_ica_data(subj_idx)=0;
            total_ICs(subj_idx)=0;
            ICs_removed{subj_idx}='0';
            total_epochs_before_artifact_rejection(subj_idx)=0;
            total_epochs_after_artifact_rejection(subj_idx)=0;
            total_channels_interpolated(subj_idx)=0;
            subj_idx=subj_idx+1;
            continue % ignore rest of the processing and go to next datafile
        end

        % Find the artifacted epochs across all channels and reject them before doing ICA.
        EEG_copy = pop_eegthresh(EEG_copy,1, 1:EEG_copy.nbchan, vol_thrs(1),...
            vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax,0,0);
        EEG_copy = eeg_checkset(EEG_copy);

        % 1         : data type (1: electrode, 0: component)
        % 0         : display with previously marked rejections? (0: no, 1: yes)
        % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)

        % Find artifaceted epochs by using power threshold in 20-40Hz frequency band.
        % This method mainly rejects muscle movement (EMG) artifacts.
        %     EEG_copy = pop_rejspec(EEG_copy, 1,'elecrange', 1:EEG_copy.nbchan,...
        %         'method', 'fft', 'threshold', emg_thrs ,'freqlimits', emg_freqs_limit,...
        %         'eegplotplotallrej', 0, 'eegplotreject', 0);

        % method                : method to compute spectrum (fft)
        % threshold             : [lower upper] threshold limit(s) in dB.
        % freqlimits            : [lower upper] frequency limit(s) in Hz.
        % eegplotplotallrej     : 0 = Do not superpose rejection marks on previous marks stored in the dataset.
        % eegplotreject         : 0 = Do not reject marked trials (but store the  marks).

        % Find the number of artifacted epochs and reject them
        EEG_copy = eeg_checkset(EEG_copy);
        EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        reject_artifacted_epochs=EEG_copy.reject.rejglobal;
        EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);

        %fig=compute_and_plot_psd(EEG_copy, 1:EEG_copy.nbchan);
        %saveas(fig, fullfile(subject_output_data_dir,'08-ica_copy_epochs_psd.png'));

        %% STEP 9: Run ICA
        length_ica_data(subj_idx)=EEG_copy.trials; % length of data (in second) fed into ICA
        EEG_copy = eeg_checkset(EEG_copy);
        EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1,...
            'stop', 1E-7, 'interupt','off');

        EEG_copy = eeg_checkset(EEG_copy);
        EEG_copy = pop_editset(EEG_copy, 'setname',  strrep(data_file_name, ext, '_ica_old'));
        EEG_copy = pop_saveset(EEG_copy, 'filename', strrep(data_file_name, ext, '_ica_old.set'),...
            'filepath', [subject_output_data_dir filesep '03_ica_data' filesep ]); % save .set format

        % Find the ICA weights that would be transferred to the original dataset
        ICA_WINV=EEG_copy.icawinv;
        ICA_SPHERE=EEG_copy.icasphere;
        ICA_WEIGHTS=EEG_copy.icaweights;
        ICA_CHANSIND=EEG_copy.icachansind;

        % If channels were removed from copied dataset during preparation of ica, then remove
        % those channels from original dataset as well before transferring ica weights.
        EEG = eeg_checkset(EEG);
        EEG = pop_select(EEG,'nochannel', ica_prep_badChans);

        % Transfer the ICA weights of the copied dataset to the original dataset
        EEG.icawinv=ICA_WINV;
        EEG.icasphere=ICA_SPHERE;
        EEG.icaweights=ICA_WEIGHTS;
        EEG.icachansind=ICA_CHANSIND;
        EEG = eeg_checkset(EEG);

        %% STEP 10: Run adjust to find artifacted ICA components
        badICs=[];

        if size(EEG_copy.icaweights,1) == size(EEG_copy.icaweights,2)
            figure()
            badICs = adjusted_ADJUST(EEG_copy, [[subject_output_data_dir filesep '03_ica_data' filesep] strrep(data_file_name, ext, '_adjust_report_old')]);
            close all;
        else % if rank is less than the number of electrodes, throw a warning message
            warning('The rank is less than the number of electrodes. ADJUST will be skipped. Artefacted ICs will have to be manually rejected for this participant');
        end

        % Mark the bad ICs found by ADJUST
        for ic=1:length(badICs)
            EEG.reject.gcompreject(1, badICs(ic))=1;
            EEG = eeg_checkset(EEG);
        end
        total_ICs(subj_idx)=size(EEG.icasphere, 1);
        if numel(badICs)==0
            ICs_removed{subj_idx}='0';
        else
            ICs_removed{subj_idx}=num2str(double(badICs));
        end
        % Mark the bad 1/f found by adjust
        % if numel(cwb)==0
        %    bad_1_f{subj_idx}='0';
        %else
        %   bad_1_f{subj_idx}=num2str(double(cwb));
        %end

        %% Save dataset after ICA, if saving interim results was preferred
        EEG = eeg_checkset(EEG);
        EEG = pop_editset(EEG, 'setname',  strrep(data_file_name, ext, '_ica_data_old'));
        EEG = pop_saveset(EEG, 'filename', strrep(data_file_name, ext, '_ica_data_old.set'),...
            'filepath', [subject_output_data_dir filesep '03_ica_data' filesep ]); % save .set format

        EEG = pop_loadset('filepath', [subject_output_data_dir filesep '03_ica_data' filesep ],...
            'filename', strrep(data_file_name, ext, '_ica_data_old.set'));

        %% STEP 11: Remove artifacted ICA components from data
        all_bad_ICs=0;
        ICs2remove=find(EEG.reject.gcompreject); % find ICs to remove

        % If all ICs and bad, save data at this stage and ignore rest of the preprocessing for this subject.
        if numel(ICs2remove)==total_ICs(subj_idx)
            all_bad_ICs=1;
            warning(['No usable data for datafile', data_file_name]);
        else
            EEG = eeg_checkset( EEG );
            EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
        end

        if all_bad_ICs==1
            total_epochs_before_artifact_rejection(subj_idx)=0;
            total_epochs_after_artifact_rejection(subj_idx)=0;
            total_channels_interpolated(subj_idx)=0;
            subj_idx=subj_idx+1;
            continue % ignore rest of the processing and go to next datafile
        end

        fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
        saveas(fig, fullfile(subject_output_data_dir,'09-ica_art_rej_psd_old.png'));

        unepochedEEG=EEG;

        %% STEP 12: Segment data into fixed length epochs
        % loop over task_event_markers
        %     if task_event_markers {'grsp'}==0
        %         warning(['No grasp for subject', subj_id]);
        %         continue
        %     end

        EEG = eeg_checkset(unepochedEEG);
        EEG = pop_epoch( EEG, task_event_markers, task_epoch_length, 'epochinfo', 'yes');

        total_epochs_before_artifact_rejection(subj_idx)=EEG.trials;

        n_trials(subj_idx)=EEG.trials;

        %% Step 14: Artifact rejection
        all_bad_epochs=0;
        chans=[]; chansidx=[];chans_labels2=[];
        chans_labels2=cell(1,EEG.nbchan);
        for i=1:EEG.nbchan
            chans_labels2{i}= EEG.chanlocs(i).labels;
        end
        [chans,chansidx] = ismember(frontal_channels, chans_labels2);
        frontal_channela_idx = chansidx(chansidx ~= 0);
        badChans = zeros(EEG.nbchan, EEG.trials);
        badepoch=zeros(1, EEG.trials);
        if isempty(frontal_channela_idx)==1 % check whether there is any frontal channel in dataset to check
            warning('No frontal channels from the list present in the data. Only epoch interpolation will be performed.');
        end 
%         else
            % find artifaceted epochs by detecting outlier voltage in the specified channels list and remove epoch if artifacted in those channels
%             for ch =1:length(frontal_channela_idx)
%                 EEG = pop_eegthresh(EEG,1, frontal_channela_idx(ch), volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
%                 EEG = eeg_checkset( EEG );
%                 EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
%                 badChans(ch,:) = EEG.reject.rejglobal;
%             end
%             for ii=1:size(badChans, 2)
%                 badepoch(ii)=sum(badChans(:,ii));
%             end
%             badepoch=logical(badepoch);
%         end
% 
%         
%         
%         if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
%             all_bad_epochs=1;
%             warning(['No usable data for datafile', data_file_name]);
%         else
%             EEG = pop_rejepoch( EEG, badepoch, 0);
%             EEG = eeg_checkset(EEG);
%         end

        grsp_trials=[];
        go_trials=[];
        for trial=1:length(EEG.epoch)
            trial_evts=EEG.epoch(trial).eventtype;
            evt_times=EEG.epoch(trial).eventlatency;
            for evt=1:length(trial_evts)
                if strcmp(trial_evts{evt}, 'go') && evt_times{evt}==0
                    go_trials(end+1)=trial;
                elseif strcmp(trial_evts{evt}, 'grsp') && evt_times{evt}==0
                    grsp_trials(end+1)=trial;
                end
            end
        end
        
        if all_bad_epochs==0
            % Interpolate artifaczed data for all reaming channels
            badChans = zeros(EEG.nbchan, EEG.trials);
            % Find artifacted epochs by detecting outlier voltage but don't remove
            for ch=1:EEG.nbchan
                EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), task_epoch_length(1), 0,0,0);
                EEG = eeg_checkset(EEG);
                EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                bad_pre_trials= EEG.reject.rejglobal;
                bad_pre_trials(go_trials)=0;

                EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), 0, task_epoch_length(2),0,0);
                EEG = eeg_checkset(EEG);
                EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                bad_post_trials= EEG.reject.rejglobal;
                bad_post_trials(grsp_trials)=0;
                all_bad_trials = bad_pre_trials | bad_post_trials;
                badChans(ch,:) = all_bad_trials;

            end
            tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
            for et = 1:EEG.trials
                % Select only this epoch (e)
                EEGe = pop_selectevent( EEG, 'epoch', et, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
                badChanNum = find(badChans(:,et)==1); % find which channels are bad for this epoch
                if length(badChanNum) < round((10/100)*EEG.nbchan)% check if more than 10% are bad
                    EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                    tmpData(:,:,et) = EEGe_interp.data; % store interpolated data into matrix
                end
            end
            EEG.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file

            % If more than 10% of channels in an epoch were interpolated, reject that epoch
            badepoch=zeros(1, EEG.trials);
            for ei=1:EEG.trials
                NumbadChan = badChans(:,ei); % find how many channels are bad in an epoch
                if sum(NumbadChan) > round((10/100)*EEG.nbchan)% check if more than 10% are bad
                    badepoch (ei)= sum(NumbadChan);
                end
            end
            badepoch=logical(badepoch);

            % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
            if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                all_bad_epochs=1;
                warning(['No usable data for datafile', data_file_name]);
            else
                EEG = pop_rejepoch(EEG, badepoch, 0);
                EEG = eeg_checkset(EEG);
            end

            % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
            if sum(EEG.reject.rejthresh)==EEG.trials || sum(EEG.reject.rejthresh)+1==EEG.trials
                all_bad_epochs=1;
                warning(['No usable data for datafile', data_file_name]);
            else
                EEG = pop_rejepoch(EEG,(EEG.reject.rejthresh), 0);
                EEG = eeg_checkset(EEG);
            end
        end
        
        % if all epochs are found bad during artifact rejection
        if all_bad_epochs==1
            total_epochs_after_artifact_rejection(subj_idx)=0;
            total_channels_interpolated(subj_idx)=0;
            continue % ignore rest of the processing and go to next datafile
        else
            total_epochs_after_artifact_rejection(subj_idx)=EEG.trials;
        end

        %         trials_removed{subj_idx, e}='';
        %         for n=1:EEG.trials
        %             if numel(n_trials)>0
        %                 trials_removed(subj_idx, e)=join(cellfun(@(x) num2str(x(1)), num2cell(n_trials,3), 'UniformOutput', false)',',');
        %             end
        %         end

        %% Interpolation
        total_channels_interpolated(subj_idx)=0;
        total_channels_interpolated(subj_idx)=length(origEEG.chanlocs)-length(EEG.chanlocs);
        EEG = pop_interp(EEG, origEEG.chanlocs, interp_type);
        fprintf('\nMissed channels are spherically interpolated\n');

        %% Re-referencing
        EEG = pop_reref( EEG, []);

        fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
        saveas(fig, fullfile(subject_output_data_dir,'11-art_rej_reref_psd_old.png'));


        %% Save processed data
        EEG = eeg_checkset(EEG);
        EEG = pop_editset(EEG, 'setname',  strrep(data_file_name, ext, '_rereferenced_data_old'));
        EEG = pop_saveset(EEG, 'filename', strrep(data_file_name, ext, '_rereferenced_data_old.set'),...
            'filepath', [subject_output_data_dir filesep '04_rereferenced_data']); % save .set format
        
        subj_idx=subj_idx+1;
    end
    
end

%% Create the report table for all the data files with relevant preprocessing outputs.
report_table=table(subj_ids',subj_ages',...
    lof_flat_channels', lof_channels', lof_periodo_channels', lof_bad_channels',...
    asr_tot_samples_modified', asr_change_in_RMS', ica_preparation_bad_channels',...
    length_ica_data', total_ICs', ICs_removed', total_epochs_before_artifact_rejection',...
    total_epochs_after_artifact_rejection',...
    total_channels_interpolated');

report_table.Properties.VariableNames={'subject','age',...
    'lof_flat_channels', 'lof_channels','lof_periodo_channels', 'lof_bad_channels'...
    'asr_tot_samples_modified', 'asr_change_in_RMS','ica_preparation_bad_channels'...
    'length_ica_data', 'total_ICs', 'ICs_removed','total_epochs_before_artifact_rejection',...
    'total_epochs_after_artifact_rejection',...
    'total_channels_interpolated'};
writetable(report_table, fullfile(study_info.data_dir, 'data','derivatives', pipeline, age, [sprintf('NEARICA_preprocessing_report_old'), datestr(now,'dd-mm-yyyy'),'.csv']));