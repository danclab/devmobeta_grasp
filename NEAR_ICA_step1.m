function NEAR_ICA_step1(study_info)


%% Clear variable space and run eeglab

ext='.set';

% Enter the path of the channel location file
channel_locations = 'C:\Users\mgautier\Desktop\DEVMOBETA\eeg_data\devmobeta\data\0_2AverageNet128_v1.sfp';

% set sampling rate (in Hz), if you want to down sample
sampling_rate = 500;

% recommended list for EGI 128 channel net: {'E17' 'E38' 'E43' 'E44' 'E48'
% 'E49' 'E113' 'E114' 'E119' 'E120' 'E121' 'E125' 'E126' 'E127' 'E128'
% 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107'}

% Initialize the filters
% High-pass frequency
highpass = 1;
% Low-pass frequency. We recommend low-pass filter at/below line noise
% frequency (see manuscript for detail)
lowpass  = 100;


%% Loop over all data files
for s_idx=1:size(study_info.participant_info,1)
    % Get subject ID from study info
    subject=study_info.participant_info.participant_id{s_idx};
    
    % Where original raw data is located
    subject_raw_data_dir=fullfile(study_info.data_dir, 'data','raw', subject, 'eeg');
    
    % Where to put processed (derived) data
    subject_output_data_dir=fullfile(study_info.data_dir, 'data','derivatives', 'NEARICA', subject);
    
    if exist([subject_output_data_dir filesep '01_filtered_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '01_filtered_data'])
    end

    if exist([subject_output_data_dir filesep '02_near_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '02_near_data'])
    end
        
    if exist([subject_output_data_dir filesep '03_ica_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '03_ica_data'])
    end

    if exist([subject_output_data_dir filesep '04_rereferenced_data'], 'dir') == 0
        mkdir([subject_output_data_dir filesep '04_rereferenced_data'])
    end
    
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', subject);
            
    %% Step 2a: Import data
    raw_fname=dir(fullfile(subject_raw_data_dir,'*.mff'));
    EEG = pop_mffimport({fullfile(subject_raw_data_dir, raw_fname.name)},{'code'},0,0);    
    data_file_name=sprintf('%s_task-%s_eeg.set',subject, study_info.task);    
    
    %% STEP 1.5: Delete discontinuous data from the raw data file
    % (OPTIONAL, but necessary for most EGI files)
    % Note: code below may need to be modified to select the appropriate
    % markers (depends on the import function) remove discontinous data at
    % the start and end of the file
    % boundary markers often indicate discontinuity
    disconMarkers = find(strcmp({EEG.event.type}, 'boundary'));
    if length(disconMarkers)>0
        % remove discontinuous chunk
        EEG = eeg_eegrej( EEG, [1 EEG.event(disconMarkers(1)).latency] );
        EEG = eeg_checkset( EEG );
    end
    % remove data after last task event (OPTIONAL for EGI files... useful when
    % file has noisy data at the end)
    end_flags=find(strcmp({EEG.event.type},'grsp'));
    latency=EEG.event(end_flags(end)).latency;
    % remove everything 3 seconds after the last event
    EEG = eeg_eegrej( EEG, [(latency+(3*EEG.srate)) EEG.pnts] );
    EEG = eeg_checkset( EEG );
    
    start_flags = find(strcmp({EEG.event.type},'strt'));
    latency=EEG.event(start_flags(1)).latency;
    % remove everything until 3 seconds before the first event
    EEG = eeg_eegrej( EEG, [1 (latency-(3*EEG.srate))] );
    EEG = eeg_checkset( EEG );
    
    %% STEP 2: Import channel locations
    EEG=pop_chanedit(EEG, 'load',{channel_locations 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );
    
    % Check whether the channel locations were properly imported. The EEG
    % signals and channel numbers should be same.
    if size(EEG.data, 1) ~= length(EEG.chanlocs)
        error('The size of the data does not match with channel numbers.');
    end
    
    % Plot channel layout
    fig=figure();
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(fig, fullfile(subject_output_data_dir,'01-initial_ch_locations.png'));
    
    %% STEP 3: Adjust anti-aliasing and task related time offset
%     % adjust anti-aliasing filter time offset
%     if filter_timeoffset~=0
%         for aafto=1:length(EEG.event)
%             EEG.event(aafto).latency=EEG.event(aafto).latency+(filter_timeoffset/1000)*EEG.srate;
%         end
%     end        
    
    %% STEP 4: Change sampling rate
    if floor(sampling_rate) > EEG.srate
        error ('Sampling rate cannot be higher than recorded sampling rate');
    elseif floor(sampling_rate) ~= EEG.srate
        EEG = pop_resample( EEG, sampling_rate);
        EEG = eeg_checkset( EEG );
    end
    
    %% STEP 5: Delete outer layer of channels
%     chans_labels=cell(1,EEG.nbchan);
%     for i=1:EEG.nbchan
%         chans_labels{i}= EEG.chanlocs(i).labels;
%     end
%     if delete_outerlayer==1
%         [chans,chansidx] = ismember(study_info.outerlayer_channel, chans_labels);
%         outerlayer_channel_idx = chansidx(chansidx ~= 0);
%         if isempty(outerlayer_channel_idx)==1
%             error(['None of the outer layer channels present in channel locations of data.'...
%                 ' Make sure outer layer channels are present in channel labels of data (EEG.chanlocs.labels).']);
%         else
%             fig=compute_and_plot_psd(EEG,outerlayer_channel_idx);
%             saveas(fig, fullfile(subject_output_data_dir,'02-outer_ch_psd.png'));
%             
%             EEG = pop_select( EEG,'nochannel', outerlayer_channel_idx);
%             EEG = eeg_checkset( EEG );
%         end
%     end
    
    % Plot channel locations
    fig=figure();
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(fig, fullfile(subject_output_data_dir,'03-inner_ch_locations.png'));
    
    % Plot PSD
    fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'04-inner_ch_psd.png'));    
    
    %% STEP 6: Filter data
    % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
    % df = transition band width, dF = normalized transition width, fs = sampling rate
    % dF is specific for the window type. Hamming window dF = 3.3
    
    high_transband = highpass; % high pass transition band
    low_transband = 10; % low pass transition band
    
    hp_fl_order = 3.3 / (high_transband / EEG.srate);
    lp_fl_order = 3.3 / (low_transband / EEG.srate);
    
    % Round filter order to next higher even integer. Filter order is always even integer.
    if mod(floor(hp_fl_order),2) == 0
        hp_fl_order=floor(hp_fl_order);
    elseif mod(floor(hp_fl_order),2) == 1
        hp_fl_order=floor(hp_fl_order)+1;
    end
    
    if mod(floor(lp_fl_order),2) == 0
        lp_fl_order=floor(lp_fl_order)+2;
    elseif mod(floor(lp_fl_order),2) == 1
        lp_fl_order=floor(lp_fl_order)+1;
    end
    
    % Calculate cutoff frequency
    high_cutoff = highpass/2;
    low_cutoff = lowpass + (low_transband/2);
    
    % Performing high pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass',...
        'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % pop_firws() - filter window type hamming ('wtype', 'hamming')
    % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
    
    % Performing low pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass',...
        'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );
    
    % pop_firws() - transition band width: 10 Hz
    % pop_firws() - filter window type hamming ('wtype', 'hamming')
    % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
    
    % Plot PSD
    fig=compute_and_plot_psd(EEG,1:EEG.nbchan);
    saveas(fig, fullfile(subject_output_data_dir,'05-filtered_psd.png'));
    
    %% Save data after running filter and LOF function, if saving interim results was preferred
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', strrep(data_file_name, ext, '_filtered_data'));
    EEG = pop_saveset( EEG,'filename',strrep(data_file_name, ext, '_filtered_data.set'),...
        'filepath', [subject_output_data_dir filesep '01_filtered_data' filesep]); % save .set format
end
end