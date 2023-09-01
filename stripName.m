function newOne = stripName(one)
%function newOne = stripName(one)
%the metabolite names are really tedious, so I am going to strip all the
%annoying pieces out so that the matching is easier. Do this as a function
%and then set it loose everyone
%KL 10/5/2022; the BC lists have yet another special character 3/1/2023
%note - the easist way to find new oddball characters is to get a string of
%char with the special character and do this:
%double(str) -- there will be one number per char, and you can use that
%below in the first if statement to delete the special character
%There are three steps needed to make these names manageable
%1. strip out all the apostrophes --> char(39)
%1b. strip out the hyphens --> char(45)
%1c. strip out the apostrophe --> char(8242)
%2. make names all lower
%3. deblank the names
if ischar(one)  
    %char will only be one metabolite
    %do nothing...already ready for next step
    one(one==char(39))=[];
    one(one==char(45))=[];
    one(one==char(8242))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(226))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(8364))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(178))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(32))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(40))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(41))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(43))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(44))=[]; %3/1/2023 - this is in the BC list (eyeroll)
    one(one==char(63))=[]; %3/1/2023 - this is in the BC list (eyeroll)

    newOne = lower(deblank((one)));    
else
    %seems to be a cell - could be one metabolite, or a cell array of them
    %do the stripping and then give back a cell
    if length(one)==1
        one = one{:};
        one(one==char(39))=[];
        one(one==char(45))=[];
        newOne = {lower(deblank((one)))};
    else
        %have multiple names in a cell, each row is a char
        newOne = cell(size(one));
        for a = 1:size(one,1)
            for j = 1:size(one,2) 
                t = one{a,j};
                t(t==char(39))=[];
                t(t==char(45))=[];
                t(t==char(8242))=[];
                t(t==char(226))=[];
                t(t==char(8364))=[];
                t(t==char(178))=[];
                t(t==char(32))=[];
                t(t==char(40))=[];
                t(t==char(41))=[];
                t(t==char(43))=[];
                t(t==char(44))=[];
                t(t==char(63))=[];
                newOne{a,j} = lower(deblank(t));
                clear t
            end
        end
        clear a
    end
end
