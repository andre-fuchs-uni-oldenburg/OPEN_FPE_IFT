% Ask the user to select an option
function result = askOptions

answer = questdlg('Which normalization should be used to calculate the spectral density?', ...
'Normalization', ...
'power spectral density (ESD)','energy spectral density (PSD)','energy spectral density (PSD)');

% Handle response
switch answer
    case 'power spectral density (ESD)'
        result = 1;
    case 'energy spectral density (PSD)'
        result = 0;   
end

% % Ask the user to select an option, using the listdlg function
% function result = askOptions(cl_questions, cl_options)
%     % return char, for consistency
%     result = num2str(listdlg('PromptString', cl_questions{1}, ...
%                      'SelectionMode', 'single', ...
%                      'ListString', cl_options, ...
%                      'Name', 'Input for TurbTools', ...
%                      'ListSize', [300 100]));
% end