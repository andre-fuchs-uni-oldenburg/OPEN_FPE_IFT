% Ask input from the user with the inputdlg function. A single cell result is converted to a char
        function result = askYesno(c_question, c_default)
            result = questdlg(c_question, 'Input', 'Yes', 'No', c_default);          
            result = double(strcmp(result,'Yes'));          
        end
        
% Ask the user to select an option, using the listdlg function
        function result = askOptions(cl_questions, cl_options)
            
            % return char, for consistency
            result = num2str(listdlg('PromptString', cl_questions{1}, ...
                             'SelectionMode', 'single', ...
                             'ListString', cl_options, ...
                             'Name', 'Input for TurbTools', ...
                             'ListSize', [300 100]));
            
        end