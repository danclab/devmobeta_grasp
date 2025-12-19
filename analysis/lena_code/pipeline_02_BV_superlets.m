function pipeline_02_superlets()

raw_data_dir='/home/common/bonaiuto/devmobeta/data/';
deriv_data_dir='/home/common/bonaiuto/devmobeta/derivatives/';
addpath('/home/bonaiuto/eeglab2021.1');
eeglab;


participants_file = '/home/common/bonaiuto/devmobeta/derivatives/participants_v2.csv';
participants = readtable(participants_file);  


subs=dir(raw_data_dir);
subs=subs(3:end);

for s_idx=1:length(subs)
    subject=subs(s_idx).name;
    %diary('/home/ldurieux/Desktop/disp_superlet.txt');
    %diary on;

    subject_info = participants(strcmp(participants.subject_id, subject), :);  
    
    sessions=dir(fullfile(raw_data_dir, subject));
    sessions=sessions(3:end);
    
    for ses_idx=1:length(sessions)
      session=sessions(ses_idx).name;
         if strcmp(subject, 'sub-276') && strcmp(session, 'ses-03')
             fprintf('\n Traitement de %s %s\n', subject, session);
             if exist(fullfile(raw_data_dir, subject, session, 'eeg'))==7       
    %for ses_idx=1:length(sessions)
        %session=sessions(ses_idx).name;
        
        %if exist(fullfile(raw_data_dir, subject, session, 'eeg')) == 7
            session_info = participants(strcmp(participants.subject_id, subject) & ...
                                        strcmp(participants.session, session), :);
            if isempty(session_info)
                continue; % Passer à la session suivante si pas d'infos trouvées
            
            end
            system = session_info.eeg_system{1};  
            disp(['Processing session: ' session ' with EEG system: ' system]);

            % Si le système est  BrainVision, on continue à traiter la session
            if any(strcmp(system, {'BrainVision'}))
                pipeline='NEARICA_NF';
                
                
                subject_output_data_dir=fullfile(deriv_data_dir, subject, session, 'eeg', pipeline);
                disp(subject_output_data_dir);
                
                fname=fullfile(subject_output_data_dir, '04_rereferenced_data', sprintf('%s_task-devmobeta_grasp_eeg_rereferenced_data.set',subject));
                if exist(fname)==2
                    disp(fname);
                    foi=linspace(5,40,50);
                    
                    
                    EEG=pop_loadset(fname);
                    
                    disp('Événements détectés :');
                    disp(unique({EEG.event.type}));
                    
                    event_map = struct('go', 'S  2', 'grsp', 'S  3');
                    for epoch_type = fieldnames(event_map)'
                        epoch_name = epoch_type{1}; % 'go' ou 'grsp'
                        event_label = event_map.(epoch_name); % 'S  2' ou 'S  3
                        
                         trial_idx = find(strcmp(event_label, {EEG.event.type}));
                         disp(['Nombre d''essais pour ', epoch_name, ' : ', num2str(length(trial_idx))]);
                    
                         if isempty(trial_idx)
                           warning(['Aucun essai trouvé pour ', epoch_name, ' !']);
                           continue;
                        end
                    %epochs = {'go', 'grsp'};
                   % disp (epochs);
                    %for ep_idx = 1:length(epochs)
                       % epoch_type = epochs{ep_idx};
                       % disp(epoch_type)
                        
                        %trial_idx=find(strcmp({EEG.event.type},epoch_type));
                        
                        trial_tf=zeros(length(trial_idx), length(foi), EEG.nbchan, EEG.pnts);
                        for c=1:EEG.nbchan
                            for t=1:length(trial_idx)
                                tf=abs(faslt(double(squeeze(EEG.data(c,:,trial_idx(t)))), EEG.srate, foi, 4, [5 40], 1));
                                trial_tf(t,:,c,:)=tf;
                                disp('processing')
                            end
                        end
                        save(fullfile(subject_output_data_dir, sprintf('%s_%s_processed_superlet_tf.mat', subject, epoch_name)), 'foi', 'trial_tf', '-v7.3');

                    end
                end
            else
                disp(['Skipping subject ' subject ' because it is not using BrainVision EEG system.']);
                %diary off;
            end
        end      
         end
    end
end
