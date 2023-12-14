function final_zapline(study_info, age)

addpath('NoiseTools');
ext='.set';
pipeline='NEARICA_NF'; 

subj_idx=1;

for a_idx=53:53%size(study_info.participant_info,1)
    subj_age = study_info.participant_info.age{a_idx};
    if strcmp(subj_age,age)

    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{a_idx};

    % Path containing subject data
    subject_dir=fullfile(study_info.data_dir,'data', 'derivatives', pipeline, age, subj_id);
    subject_data_dir=fullfile(subject_dir, '04_rereferenced_data');
    subject_output_data_dir=fullfile(study_info.data_dir, 'data','derivatives', pipeline, age, subj_id);
    data_file_name=sprintf('%s_task-%s_eeg.set',subj_id, study_info.task);

    fname=sprintf('%s_task-grasping_eeg_rereferenced_data_old.set',subj_id);

    % Load data
        if exist(fullfile(subject_data_dir,sprintf('%s_task-grasping_eeg_rereferenced_data_old.set',subj_id)),'file')==2

            EEG=pop_loadset('filepath', subject_data_dir,...
                'filename', fname);   

            load(fullfile(subject_dir,sprintf('%s_zapline_epoched_old_new.mat', subj_id)));

            data=data_epoched_old_new.*1e6;

            EEG.data=permute(data,[2 3 1]);
            %EEG.data=permute(data,[2 1 3]);%for old epoched zapline
            EEG = pop_editset(EEG, 'setname', strrep(data_file_name, ext, '_final_zapped_old_new'));
            EEG = pop_saveset( EEG,'filename',strrep(data_file_name, ext, '_final_zapped_old_new.set'),...
            'filepath', [subject_output_data_dir filesep '05_final_zapped_data' filesep]); % save .set format

            %fig=compute_and_plot_psd(EEG, 1:EEG.nbchan);
            %saveas(fig, fullfile(subject_output_data_dir,'final-zapped.png'));
        else 
            continue 
        end
        
    subj_idx=subj_idx+1;

    end 
end    