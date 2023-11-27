function dir_name = create_unique_directory(base_dir_name)
    % Initialize the directory name with the base name
    dir_name = base_dir_name;
    % Initialize a suffix with an empty string
    suffix = '';
    % Initialize a counter
    i = 1;
    % Letters for suffixes
    alphabet = 'abcdefghijklmnopqrstuvwxyz';

    % Check if the directory exists and create a new name if it does
    while exist([dir_name suffix], 'dir')
        % Move to the next letter in the alphabet for the suffix
        suffix = ['_' alphabet(i)];
        i = i + 1;
        % Wrap around the alphabet if necessary
        if i > length(alphabet)
            error('Ran out of letters to append.');
        end
    end

    % When a unique name is found, create the directory
    mkdir([dir_name suffix]);
    
    % Update dir_name to the unique name
    dir_name = [dir_name suffix];
end
