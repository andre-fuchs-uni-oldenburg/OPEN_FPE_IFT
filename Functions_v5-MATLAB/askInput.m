% Ask input from the user with the inputdlg function. A single cell result is converted to a char
        function result = askInput(cl_questions, cl_defaults)
            if length(cl_questions) == length(cl_defaults)
                temp = inputdlg(cl_questions, 'Input', 1, cl_defaults);
                
                if length(temp) > 1
                    result = temp;
                else
                    % Let's return a char
                    result = str2double(temp);
                end
            else
                error('The number of questions should equal the number of defaults');
            end
        end
        
