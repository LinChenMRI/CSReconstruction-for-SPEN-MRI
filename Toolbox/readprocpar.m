function [ProcparStruct] = readprocpar(ProcparFilePath)
   
% [ProcparStruct] = ReadProcparFile(ProcparFilePath)
% 
% Reads a procpar file and stores it inside a Matlab structure
%
% Note 1: Not fully tested. Worked on a simple(?) procpar.
% Note 2: Not all consistancy tests are performed, e.g. that number of
%  values, or possible values, match the number set in the file (first
%  number in 2nd line, or in last line.
% Note 3: Not written nicley, should probably be rewritten.
%
% Follows definitions in "VrmrJ User Programming" chapter 5, page 273 
%  (255 in PDF).


   % Define Line type Constants
   NONE = 0 ;
   FIRST = 1 ;
   SECOND = 2 ;
   MID = 3 ;
   LAST = 4 ;
   
   % Define perl search expressions
   StringStartWithLetter = '[a-zA-Z]\S*' ;
   WhiteSpace = '\s*' ;
   InitialWhiteSpace = '^\s*' ;
   AnyString = '.*' ;
   StringInQuotes = '"[^"]*"' ;
   TokenizedStringsIQuotes = '"([^"]*)"' ;
   UintString = '\d*' ; % "unsigned integer", i.e. digits only.
   BlankLine = '^\s*$' ;
   
   % Variable types (as in procpar definitions)
   TYPE_UNKNOWN = 0 ;
   TYPE_REAL = 1 ; 
   TYPE_STR = 2 ;
   
   
   fid = fopen([ProcparFilePath,'\procpar']) ;
   
   LastLineType = NONE ; % NONE, FIRST, MID, LAST
   TextLine = fgetl(fid) ;
   
   
   while (TextLine ~= -1)
      
      if (~isempty(regexp(TextLine, BlankLine)) )
         % do nothing - skip blank line
         
      elseif ( (LastLineType == NONE || LastLineType == LAST) && ...
               ~isempty(regexp(TextLine,[InitialWhiteSpace, ...
                                         StringStartWithLetter])) )
                        
%          [ProcparStruct, VarName]= GetFirstLine(ProcparStruct, TextLine) ;
         
         TokensCell = regexp(TextLine,[InitialWhiteSpace, ...
                                       '(' StringStartWithLetter ')', ...
                                       '(' AnyString ')'], 'tokens' ) ;
         
         VarName = TokensCell{1}{1} ;
         VarParams = sscanf(TokensCell{1}{2}, '%f') ;
         
         ProcparStruct.(VarName).name = VarName ;
         ProcparStruct.(VarName).subtype = VarParams(1) ;
         ProcparStruct.(VarName).basictype = VarParams(2) ;
         ProcparStruct.(VarName).maxvalue = VarParams(3) ;
         ProcparStruct.(VarName).minvalue = VarParams(4) ;
         ProcparStruct.(VarName).stepsize = VarParams(5) ;
         ProcparStruct.(VarName).Ggroup = VarParams(6) ;
         ProcparStruct.(VarName).Dgroup = VarParams(7) ;
         ProcparStruct.(VarName).protection = VarParams(8) ;
         ProcparStruct.(VarName).active = VarParams(9) ;
         ProcparStruct.(VarName).intptr = VarParams(10) ;
         
         LastLineType = FIRST ; 
         
      elseif ( LastLineType == FIRST && ...
               ~isempty(regexp(TextLine,[InitialWhiteSpace, ...
                                         UintString])) )
                             
         switch ProcparStruct.(VarName).basictype
            case TYPE_REAL % Real valued variable
               ValuesTmp = sscanf(TextLine, '%f') ;
               ProcparStruct.(VarName).Line2Prefix = ValuesTmp(1) ;
               if length(ValuesTmp) > 1
                  ProcparStruct.(VarName).Values = ValuesTmp(2:end) ;
               else
                  ProcparStruct.(VarName).Values = [] ;
               end
            case TYPE_STR % String valued variable
               
               TokensCell = regexp(TextLine, ...
                                   [InitialWhiteSpace, ...
                                    '(' UintString ')', ...
                                    WhiteSpace, ...
                                    '(' StringInQuotes ')'], 'tokens' ) ;
               Line2Prefix = num2str(TokensCell{1}{1}) ;
               if length(TokensCell{1}) > 1
                  DoubleQuotedString = TokensCell{1}{2} ;
                  % discard double quotes
                  ValuesCounter = 1 ;
                  Values = DoubleQuotedString(2:(end-1)) ;
                  ProcparStruct.(VarName).Values{ValuesCounter} = Values ;
               else
                  Values = [] ;
                  ValuesCounter = 0 ;
                  ProcparStruct.(VarName).Values = {} ;
               end
               ProcparStruct.(VarName).Line2Prefix = str2double(Line2Prefix) ;
            otherwise
               error('Unexpected variable type')
         end
         
         LastLineType = SECOND ;
         
      elseif ( (LastLineType == SECOND || LastLineType == MID) && ...
               ProcparStruct.(VarName).basictype == TYPE_STR && ...
               ProcparStruct.(VarName).Line2Prefix > 1 && ...
               ~isempty(regexp(TextLine,[InitialWhiteSpace, ...
                                         StringInQuotes])) )
                             
         TokensCell = regexp(TextLine, ...
                             [InitialWhiteSpace, ...
                              '(' StringInQuotes ')'], 'tokens' ) ;
               
         DoubleQuotedString = TokensCell{1}{1} ;
         ValuesCounter = ValuesCounter + 1 ;
         Values = DoubleQuotedString(2:(end-1)) ;
         ProcparStruct.(VarName).Values{ValuesCounter} = Values ;
                             
         LastLineType = MID ;
            
      elseif ( (LastLineType == SECOND || LastLineType == MID) && ...
               ~isempty(regexp(TextLine,[InitialWhiteSpace, ...
                                         UintString])) )
                             
                             
         TokensCell = regexp(TextLine,[InitialWhiteSpace, ...
                                       '(' UintString ')', ...
                                       '(' AnyString ')'], 'tokens' ) ;
         NumPossibleValues = str2num(TokensCell{1}{1}) ;
         ProcparStruct.(VarName).NumPossibleValues = NumPossibleValues ;
         
         if NumPossibleValues == 0
            ProcparStruct.(VarName).PossibleValues = [] ;
         elseif (NumPossibleValues > 0 && ...
                 ProcparStruct.(VarName).basictype == TYPE_STR)
            ValuesTokenCells = regexp(TokensCell{1}{2}, ...
                                      TokenizedStringsIQuotes, 'tokens' ) ;
            ProcparStruct.(VarName).PossibleValues = ...
                                              horzcat(ValuesTokenCells{:}) ;
         elseif (NumPossibleValues > 0 && ...
                 ProcparStruct.(VarName).basictype == TYPE_REAL)
            ProcparStruct.(VarName).PossibleValues = ...
                                    sscanf(TokensCell{1}{2}, '%f') ;
                          
              
         else
            error('Unsupported case of possible values line') ;
         end
         
         LastLineType = LAST ;
      end
      
      
      
      TextLine = fgetl(fid) ;
   end
   
   fclose(fid) ;

return ;



